def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str = 'oneline.fasta'):
    '''Function convert_multiline_fasta_to_oneline(input_fasta: str,
    output_fasta: str = 'oneline.fasta') take two arguments:
    input_fasta -absolut path (row line r'path') to fasta file with
    multiple lines of nucs
    output_fasta -is not obligate, name of output file with nucs writed in
    one line, as default file name is 'oneline_fasta.fasta'
    Write new file in the same directory as input_fasta
    '''
    import os
    input_file_path = os.path.dirname(input_fasta)
    new_file_path = os.path.join(input_file_path, output_fasta)
    '''Open and reading input_file'''
    with open(input_fasta, 'r') as input_file:
        counter = 0
        res_nuc = ''
        for line in input_file:
            if line[0] == '>':
                if counter >= 1:
                    with open(new_file_path, 'a') as output_file:
                        output_file.write(key + res_nuc + '\n')
                    res_nuc = ''
                key = line
            else:
                res_nuc += line.strip()
                counter += 1
        '''Open new output file and writing lines one by one'''
        with open(new_file_path, 'a') as output_file:
            output_file.write(key + res_nuc + '\n')


def parse_blast_output(input_file: str, output_file: str):
    '''Function parse_blast_output(input_file: str, output_file: str) take
    two arguments:
    input_file -absolut path (row line r'path') to BLAST file
    output_file - name of output file with result protein list,
    placed in one column. Write new file in the same directory as
    input_file
    '''
    import os
    input_file_path = os.path.dirname(input_file)
    new_file_path = os.path.join(input_file_path, output_file)
    protein_list = []
    flag = 0
    '''Open and reading fist column in the input file'''
    with open(input_file, 'r') as blast_file:
        for line in blast_file:
            ''''Take first line and add to  protein_list'''
            if flag == 1:
                protein_list.append(line[:65].strip())
                flag = 0
            if 'Description' in line:
                flag = 1
    '''Protein_list sorted in alphabetic order '''
    protein_list.sort(key=lambda x: x.lower())
    '''Open new output file and write there protein_list in one column'''
    with open(new_file_path, 'a') as file_with_proteins:
        for protein in protein_list:
            file_with_proteins.write(protein + '\n')
