def ch_length_bounds(seq, length_bounds):
    """Ф-ия проверки соответсвия интервалу длины рида"""
    if type(length_bounds) is int:
        length_bounds = [0, length_bounds]
    return length_bounds[0] <= len(seq[0]) and len(seq[0]) <= length_bounds[1]


def ch_gc_bounds(seq, gc_bounds):
    """Ф-ия проверки соответсвия интервалу состава GC в процентах"""
    if type(gc_bounds) is int:
        gc_bounds = [0, gc_bounds]
    GC_cont = (seq[0].count('G')+seq[0].count('C'))*100/len(seq[0])
    return GC_cont >= gc_bounds[0] and GC_cont <= gc_bounds[1]


def ch_quality_threshold(seq, quality_threshold):
    """Ф-ия проверки соответсвия пороговому значению качества среднего рида"""
    sum = 0
    for symb in seq[1]:
        sum += ord(symb)
    return sum/len(seq[1]) >= quality_threshold+33
