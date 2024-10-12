Program Genom  and bio_files_processor
Developer:
Kichigina Natalia
Programm description :
Scrypt genom consist of two functions:

-Function filter_fastq(input_fastq, output_fastq, gc_bounds,length_bounds,
quality_threshold) returns new file with fastq sequnces according to filtering parameteres.  Function create new file in the 'filtered' directory in same directory as input_fastq.
Function get 5 arguments:
1) input_fastq -absolut path (row line r'path') to input fastq file;
2) output_fastq - name of output file with filtered fastq seqs.
        Write new file in 'filtered' directory in the same directory as
        input_file.
3) gc_bounds - interval of quantity of GC in percenst, default
        gc_bound = (0, 100)from 0% to 100%, if gc_bound = (20, 80)-
        from 20% to 80% including. In case of one number, it counts
        that it is high limit, gc_bound = 45 - quantity of GC lower or
        equal 45%;
4) length_bounds - interval of read length, default
        length_bounds = (0, 2**32),from 0 to 2 in power 32.
        In case of one number, it counts that it is high limit,
        length_bounds = 32 - from 0 to 32;
5) quality_threshold - treshhold for quality of mean read,
        default 0(scale Phred33). Reads with quality lower then
        treshold are decline.

Function run_dna_rna_tools() returns nuc acid sequence after actions on initial RNA or DNA nucs sequence ordered by the procedure wrote in argument. Function make check for exsistance of initial nucs sequences, by checking simultanously presence of T and U in seq.In case of negarive result func returns:
    "Please insert correct sequence"
Possible variant of func colling:
 run_dna_rna_tools('ATG','aT','reverse')
'ATG','aT',... - optional quantity of arguments, consisting of nuc sequences in str format.
        'reverse' - always last parameter in function, with instruction for required procedure.

Procedure list:
            'reverse' - return reversed sequence.
      'transcribe' - return transcribed sequence.
      'complement' - return comlemented sequence.
      'reverse_complement' - returns reversed and complemented seq.
Return value: string with one sequence of nucs or list with multiple  seqs of nucs, depending on quantity of initial seqs.

Script bio_files_processor consist to functions:
    Function convert_multiline_fasta_to_oneline(input_fasta: str,  output_fasta: str = 'oneline.fasta') take two arguments:
input_fasta -absolut path (row line r'path') to fasta file with  multiple lines of nucs
 output_fasta -is not obligate, name of output file with nucs writed in one line, as default file name is 'oneline_fasta.fasta'
    Write new file in the same directory as input_fasta

    Function parse_blast_output(input_file, output_file) take two arguments:
input_file -absolut path (row line r'path') to BLAST file
 output_file - name of output file with result protein list, placed in one column. Write new file in the same directory as input_file.
