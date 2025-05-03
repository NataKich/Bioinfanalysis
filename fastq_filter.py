import Bio.SeqIO
import Bio.SeqUtils
import os
import click
import logging
import sys


logger = logging.getLogger()


def _in_bounds(value, bounds):
    if type(bounds) is int:
        bounds = [0, bounds]
    return bounds[0] <= value <= bounds[1]


def filter_fastq(input_fastq: str, output_fastq: str,
                 gc_bounds: int | tuple = (0, 100),
                 length_bounds: int | tuple = (0, 2**32),
                 quality_threshold: int = 0):
    """
    Function filter_fastq(input_fastq, output_fastq, gc_bounds,length_bounds,
    quality_threshold) returns new file with fastq sequnces according to
    filtering parameteres.
    Function create new file in the 'filtered' directory in same directory as
    input_fastq.
    Function get 5 arguments:
        input_fastq -absolut path (row line r'path') to input fastq file;
        output_fastq - name of output file with filtered fastq seqs.
        Write new file in 'filtered' directory in the same directory as
        input_file.
        gc_bounds - interval of quantity of GC in percenst, default
        gc_bound = (0, 100)from 0% to 100%, if gc_bound = (20, 80)-
        from 20% to 80% including. In case of one number, it counts
        that it is high limit, gc_bound = 45 - quantity of GC lower or
        equal 45%;
        length_bounds - interval of read length, default
        length_bounds = (0, 2**32),from 0 to 2 in power 32.
        In case of one number, it counts that it is high limit,
        length_bounds = 32 - from 0 to 32;
        quality_threshold - treshhold for quality of mean read,
        default 0(scale Phred33). Reads with quality lower then
        treshold are decline.
    """
    if not os.path.isfile(input_fastq):
        logger.error("ERROR: Can not find file %s", input_fastq)
        raise Exception(f'Can not find file {input_fastq}')     

    parent_dir = os.path.dirname(input_fastq)
    new_file_path = os.path.join(parent_dir, 'filtered', output_fastq)

    c = 0
    
    with open(new_file_path, "w") as output:
        for rec in Bio.SeqIO.parse(input_fastq, "fastq"):
            seq = rec.seq
            quality = rec.letter_annotations["phred_quality"]

            cond1 = _in_bounds(len(seq), length_bounds)
            cond2 = _in_bounds(Bio.SeqUtils.gc_fraction(seq) * 100, gc_bounds)
            cond3 = sum(quality) / len(quality) >= quality_threshold
            if cond1 and cond2 and cond3:
                Bio.SeqIO.write(rec, output, "fastq")
                c += 1

    logger.info("INFO from %s to %s filteren %s records", input_fastq, output_fastq, c)


class BiologicalSequence:
    def __len__(self):
        pass
    def __getitem__(self, index):
        pass
    def __str__(self):
        pass
    def check_alphabet(self):
        pass

class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq):
        self.seq = seq
    def complement(self):
        res = []
        for c in self.seq:
            res.append(self.complement_dict[c])
        return self.__class__("".join(res))
    def reverse(self):
        return self.__class__(self.seq[::-1])
    def reverse_complement(self):
        return self.complement().reverse()    
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, index):
        return self.seq[index]
    def __str__(self):
        return self.seq
    def check_alphabet(self):
        valid = set(self.complement_dict.keys())
        for c in self.seq:
            if c not in valid:
                return False
        return True


class DNASequence(NucleicAcidSequence):
    complement_dict = {'T': 'A', 'A': 'T',
                'a': 't', 't': 'a',
                'G': 'C', 'C': 'G',
                'g': 'c', 'c': 'g'}
    
    def transcribe(self):
        s = self.seq.replace("T", "U")
        s = s.replace("t", "u")
        return RNASequence(s)


class RNASequence(NucleicAcidSequence):
    complement_dict = {'G': 'C', 'C': 'G',
            'g': 'c', 'c': 'g',
            'A': 'U', 'U': 'A',
            'a': 'u', 'u': 'a'}


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, seq):
        self.seq = seq
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, index):
        return self.seq[index]
    def __str__(self):
        return "".join(self.seq)
    def check_alphabet(self):
        valid = set("ABCDEFGHIJKLMNOPQRSTUVWYZ")
        for c in self.seq:
            if c not in valid:
                return False
        return True



@click.group()
def main():
    pass


@click.command()
@click.option("-q", "--quality", default=0, type=click.IntRange(min=0, max_open=True))
@click.option("-gc-low", default=0, type=click.IntRange(min=0, max=100))
@click.option("-gc-high", default=100, type=click.IntRange(min=0, max=100))
@click.option("-len-low", default=0, type=click.IntRange(min=0, max=2**32))
@click.option("-len-high", default=2**32, type=click.IntRange(min=0, max=2**32))
@click.argument("infile")
@click.argument("outfile")
def filter(quality, gc_low, gc_high, len_low, len_high, infile, outfile):
    try:
        filter_fastq(infile, outfile, 
                 gc_bounds=(gc_low, gc_high), 
                 length_bounds=(len_low, len_high), 
                 quality_threshold=quality)
    except Exception as e:
        print(e, file=sys.stderr)


@click.command()
@click.argument("cmd", type=click.Choice(["compl", "reverse", "complreverse", "transcribe"]))
@click.argument("seq")
def dna(cmd, seq):
    x = DNASequence(seq)
    if not x.check_alphabet():
        print("invalid DNA sequence", file=sys.stderr)
        return 1
    if cmd == "compl":
        print(x.complement())
    elif cmd == "reverse":
        print(x.reverse())
    elif cmd == "complreverse":
        print(x.reverse_complement())
    elif cmd == "transcribe":
        print(x.transcribe())


@click.command()
@click.argument("cmd", type=click.Choice(["compl", "reverse", "complreverse"]))
@click.argument("seq")
def rna(cmd, seq):
    x = RNASequence(seq)
    if not x.check_alphabet():
        print("invalid RNA sequence", file=sys.stderr)
        return 1
    if cmd == "compl":
        print(x.complement())
    elif cmd == "reverse":
        print(x.reverse())
    elif cmd == "complreverse":
        print(x.reverse_complement())


main.add_command(filter)
main.add_command(dna)
main.add_command(rna)


if __name__ == "__main__":
    main()