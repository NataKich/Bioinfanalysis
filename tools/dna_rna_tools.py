def valid(seq):
    '''Check sequence is real'''
    for nuc in seq[: -1]:
        if ('U' in nuc or 'u' in nuc) and ('T' in nuc or 't' in nuc):
            return False


def transcribe(seq):
    '''Nuc transcribe function'''
    trans = []
    for nuc in seq:
        ''''Replace uracil by thymine'''
        nuc = nuc.replace('T', 'U')
        nuc = nuc.replace('t', 'u')
        trans.append(nuc)
    if len(trans) == 1:
        '''Çheck the length of result sequence
        for output form'''
        return trans[0]
    else:
        return trans


def reverse(seq):
    '''Reverse sequence function'''
    rev = []
    for nuc in seq:
        '''Reverse sequence and add to list'''
        rev.append(nuc[len(nuc):: -1])
    if len(rev) == 1:
        '''Çheck the length of result sequence
        for output form'''
        return rev[0]
    else:
        return rev


def complement(seq):
    '''Complement function, find pair for
    corresponding nuc'''
    comp = []
    for nuc in seq:
        if 'U' in nuc or 'u' in nuc:
            rna_dict = {'G': 'C', 'C': 'G',
                        'g': 'c', 'c': 'g',
                        'A': 'U', 'U': 'A',
                        'a': 'u', 'u': 'a'}
            x = rna_dict
        else:
            dna_dict = {'T': 'A', 'A': 'T',
                        'a': 't', 't': 'a',
                        'G': 'C', 'C': 'G',
                        'g': 'c', 'c': 'g'}
            x = dna_dict
        res_seq = ''
        for j in nuc:
            '''Getting value from dictionary'''
            res_seq += x.get(j)
        comp.append(res_seq)
    if len(comp) == 1:
        '''Çheck the length of result sequence
        for output form'''
        return comp[0]
    else:
        return comp


def reverse_complement(seq):
    '''Reverse and Complement function, find pair for
    corresponding nuc and reversing sequence'''
    comp = []
    for nuc in seq:
        if 'U' in nuc or 'u' in nuc:
            rna_dict = {'G': 'C', 'C': 'G',
                        'g': 'c', 'c': 'g',
                        'A': 'U', 'U': 'A',
                        'a': 'u', 'u': 'a'}
            x = rna_dict
        else:
            dna_dict = {'T': 'A', 'A': 'T',
                        'a': 't', 't': 'a',
                        'G': 'C', 'C': 'G',
                        'g': 'c', 'c': 'g'}
            x = dna_dict
        res_seq = ''
        for j in nuc:
            '''Getting value from dictionary'''
            res_seq += x.get(j)
        comp.append(res_seq)
    rev = []
    for nuc in comp:
        rev.append(nuc[len(nuc):: -1])
    if len(rev) == 1:
        '''Çheck the length of result sequence
        for output form'''
        return rev[0]
    else:
        return rev
