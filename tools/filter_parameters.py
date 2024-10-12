import os


def ch_length_bounds(seq, length_bounds):
    """Function checks sequence length in read length interval"""
    if type(length_bounds) is int:
        length_bounds = [0, length_bounds]
        '''nuc sequence is first object in tuple 'seq' '''
    low_bound = length_bounds[0]
    high_bound = length_bounds[1]
    nuc_seq = seq[0]
    return low_bound <= len(nuc_seq) and len(nuc_seq) <= high_bound


def ch_gc_bounds(seq, gc_bounds):
    """Function checks of 'gc' quantity in nuc read in percents"""
    nuc_seq = seq[0]
    if type(gc_bounds) is int:
        gc_bounds = [0, gc_bounds]
        '''nuc sequence is first object in tuple 'seq'''
    GC_cont = (nuc_seq.count('G') + nuc_seq.count('C')) * 100 / len(nuc_seq)
    return GC_cont >= gc_bounds[0] and GC_cont <= gc_bounds[1]


def ch_quality_threshold(seq, quality_threshold):
    """Function checks that read has required quality """
    sum = 0
    for symb in seq[1]:
        sum += ord(symb)
        '''seq[1] - second object in tuple 'seq'''
        quality = seq[1]
    return sum/len(quality) >= quality_threshold+33


def read_fastq(input_fastq, last_position=0):
    '''Function read_fastq() receives absolute path (row line r'path')to file
    with fastq, and last position of reading in the file. Then function opens
    file and reads line from last_position, and return one line as dictionary
    seq[key] = (nuc, qual), and actual reading position input_file.tell()in
    the file.
    '''
    with open(input_fastq, 'r') as input_file:
        input_file.seek(last_position)
        key = input_file.readline().strip()
        nuc = input_file.readline().strip()
        input_file.readline()
        qual = input_file.readline().strip()
        seq = {}
        seq[key] = (nuc, qual)
        return seq, input_file.tell()


def write_fastq(write_file_path, file_name, seq):
    '''Function write_fastq() receives absolute path (row line r'path') where
    file should be create, name of new file, and dictionary seq with fastq.
    '''
    parent_dir = os.path.dirname(write_file_path)
    new_dir = 'filtered'
    file_path = os.path.join(parent_dir, new_dir)
    '''Creation of new file path with directory 'filtered' if it
    is not exsist'''
    if os.path.isdir(os.path.join(parent_dir, new_dir)):
        pass
    else:
        os.makedirs(file_path)
    for key in seq:
        '''Getting nucs and quality from tuple'''
        nucs = seq[key][0]
        quality = seq[key][1]
        with open(os.path.join(file_path, file_name), 'a') as filtered_fastq:
            filtered_fastq.write(key + '\n' + nucs +
                                 '\n' + key.replace('@', '+') +
                                 '\n' + quality + '\n')
