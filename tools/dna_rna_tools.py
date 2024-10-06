def valid(seq):  # ф-ия проверки существования посл-ти
    for d in seq[: -1]:
        if ('U' in d or 'u' in d) and ('T' in d or 't' in d):
            return False


def transcribe(seq):  # ф-ия транскрибирования
    trans = []
    for i in seq:
        i = i.replace('T', 'U')  # Замена 'T' на 'U'
        i = i.replace('t', 'u')  # Замена 't' на 'u'
        trans.append(i)
    if len(trans) == 1:
        return trans[0]  # Условие для вывода одной строки или списка
    else:
        return trans


def reverse(seq):  # ф-ия реверсирования посл-ти
    rev = []
    for i in seq:
        rev.append(i[len(i):: -1])  # Переворот посл-ти и добавление в список
    if len(rev) == 1:
        return rev[0]  # Условие для вывода одной строки или списка
    else:
        return rev


def complement(seq):  # ф-ия комплементирования
    comp = []
    for i in seq:
        if 'U' in i or 'u' in i:
            x = {'G': 'C', 'C': 'G',  # Словарь для РНК посл-ти
                 'g': 'c', 'c': 'g',
                 'A': 'U', 'U': 'A',
                 'a': 'u', 'u': 'a'}
        else:
            x = {'T': 'A', 'A': 'T',  # Словарь для ДНК посл-ти
                 'a': 't', 't': 'a',
                 'G': 'C', 'C': 'G',
                 'g': 'c', 'c': 'g'}
        res_seq = ''
        for j in i:
            res_seq += x.get(j)  # Поиск значений по словарю
        comp.append(res_seq)
    if len(comp) == 1:
        return comp[0]  # Условие для вывода одной строки или списка
    else:
        return comp


def reverse_complement(seq):  # ф-ия комплементирования и реверсирования
    comp = []
    for i in seq:
        if 'U' in i or 'u' in i:
            x = {'G': 'C', 'C': 'G',  # Словарь для РНК посл-ти
                 'g': 'c', 'c': 'g',
                 'A': 'U', 'U': 'A',
                 'a': 'u', 'u': 'a'}
        else:
            x = {'T': 'A', 'A': 'T',  # Словарь для ДНК посл-ти
                 'a': 't', 't': 'a',
                 'G': 'C', 'C': 'G',
                 'g': 'c', 'c': 'g'}
        res_seq = ''
        for j in i:
            res_seq += x.get(j)  # Поиск значений по словарю
        comp.append(res_seq)
    rev = []
    for i in comp:
        rev.append(i[len(i):: -1])
    if len(rev) == 1:
        return rev[0]  # Условие для вывода одной строки или списка
    else:
        return rev
