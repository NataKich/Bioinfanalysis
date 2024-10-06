# Genom 
 
## Описание программы :
Программа содержит две основные функции для работы с последовательностями нуклеиновых кислот.
Функция filter_fastq(seq, gc_bounds,length_bounds,quality_threshold) возвращает словарь состоящий из FASTQ-сиквенсов удовлетворящих параметрам   фильтрации. 
Функция принимает на вход 4 аргумента, seq - словарь состоящий из FASTQ-сиквенсов и параметры фильтрации:
        gc_bounds - интервал состава GC в процентах, по умолчанию
        gc_bound = (0, 100)от 0% до 100%, при вводе gc_bound = (20, 80)-
        от 20% до 80% включительно. При вводеодного числа, считается,
        что это верхняя граница gc_bound = 45 - состав GC меньше  или
        равен 45%.
        length_bounds - длина интервала для фильтрации, по умолчанию
        length_bounds = (0, 2**32),от 0 до 2 в степени 32.
        При вводе одного числа,считается, что это верхняя граница,
        length_bounds = 32 - от 0 до 32
        quality_threshold - пороговое значение качества среднего рида,
        по умолчанию равно 0(шкала Phred33). Риды со ступенью
        качества по всем нуклеотидам ниже порогового отбрасываются.
    Возвращаемое значение:словарь состоящий из FASTQ-сиквенсов
    удовлетворящих параметрамфильтрации, при наличии таковых или пустой
    при отсутсвии иквенсов удовлетворящих параметрам.
Функция run_dna_rna_tools() возвращает  последовательности нуклеиновых
 кислот после проведения указанной в аргументе процедуры над переданными в
 функцию изначальными последовательностями ДНК или РНК.
 Функция проводит проверку на возможность существованияи значальных
 последовательностей ДНК и РНК. При наличии Т и U в одно последовательности,
  функция прекращается и возвращает сообщение пользователю:
  "Please insert correct sequence"
  Вариант вызова функции:
    run_dna_rna_tools('ATG','aT','reverse')
        'ATG','aT',... - произвольное количество аргументов с
        последовательностями ДНК или РНК в формате str.
        'reverse' - всегда последний аргумент функции, с указанием процедуры
        которую необходимо выполнить в формате str.
        Список процедур:
            'reverse' - вернуть развёрнутую последовательность.
            'transcribe' - вернуть транскрибированную последовательность.
            'complement' - вернуть комплементарную последовательность.
            'reverse_complement' - вернуть обратную комплементарную.
            последовательность.
    Возвращаемое значение: строка с одной последовательностью нуклеиновых
    кислот или список последовательностей нуклеиновых кислот, в зависимости
    от количества переданных изначально последовательностей.

Разработчик:
Кичигина Наталья 