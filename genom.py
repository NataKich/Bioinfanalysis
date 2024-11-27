from tools.dna_rna_tools import (valid, transcribe, reverse,
                                 complement, reverse_complement)
from tools.filter_parameters import (ch_length_bounds, ch_gc_bounds,
                                     ch_quality_threshold, read_fastq,
                                     write_fastq)


def run_dna_rna_tools(*args) -> str | list:
    """
    Function run_dna_rna_tools() returns nuc acid sequence after actions
    on initial RNA or DNA nucs sequence ordered by the procedure
    wrote in argument.
    Function make check for exsistance of initial nucs sequences, by checking
    simultanously presence of T and U in seq.
    In case of negarive resalt func returns:
    "Please insert correct sequence"
    Possible variant of func colling:
    run_dna_rna_tools('ATG','aT','reverse')
        'ATG','aT',... - optional quantity of arguments, consisting of
        nuc sequences in str format.
        'reverse' - always last parameter in function, with istruction
        for required procedure.
        Procedure list:
            'reverse' - return reversed sequence.
            'transcribe' - return transcribed sequence.
            'complement' - return comlemented sequence.
            'reverse_complement' - returns reversed and complemented seq.
    Return value: string with one sequence of nucs or list with multiple
    seqs of nucs, depending on quantity of initial seqs.
    """
    seq: tuple = args[: -1]
    if valid(seq) is False:
        return "Please insert correct sequence"
    elif args[-1] == 'transcribe':
        return transcribe(seq)
    elif args[-1] == 'reverse':
        return reverse(seq)
    elif args[-1] == 'complement':
        return complement(seq)
    elif args[-1] == 'reverse_complement':
        return reverse_complement(seq)


def filter_fastq(input_fastq: str, output_fastq: str,
                 gc_bounds: int | tuple = (0, 100),
                 length_bounds: int | tuple = (0, 2**32),
                 quality_threshold: int = 0) -> dict:
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
    import os
    parent_dir = os.path.dirname(input_fastq)
    new_file_path = os.path.join(parent_dir, 'filtered', output_fastq)
    '''Checks the name of new file is no the same as name of other files
    in directory. If it is the same, returns 'Please insert other name
    for output file, file is exsist'
    '''
    if os.path.isfile(new_file_path):
        return 'Please insert other name for output file, file is exsist'
    location_of_reading = 0
    while True:
        '''Call read func to read the sequence from input file, and take
        position of last reading in the file'''
        seqs, location_of_reading = read_fastq(input_fastq,
                                               location_of_reading)
        '''Check the end of all seqs, break the cicle if find '' in
        seqs dictionary'''
        if '' in seqs:
            break
        '''Checking complience for all filter parameters '''
        for seq in seqs:
            cond1 = ch_length_bounds(seqs[seq], length_bounds)
            cond2 = ch_gc_bounds(seqs[seq], gc_bounds)
            cond3 = ch_quality_threshold(seqs[seq], quality_threshold)
            if cond1 and cond2 and cond3:
                write_fastq(input_fastq, output_fastq, seqs)
    return 'Filtering is complete'
