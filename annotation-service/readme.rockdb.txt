Упаковка данных

nohup java -cp annotation_v.0.6.1.3-full.jar \
    org.forome.annotation.makedatabase.main.MainMakeDatabase \
    -database /projects/annotations/database/annotation-hg38 \
    -assembly GRCh38 \
    -gerp19 /db/download/hg19.GERP_scores.tar.gz &

nohup java -cp annotation_v.0.6.1.13-full.jar \
    org.forome.annotation.makedatabase.main.MainMakeDatabase \
    -database /projects/annotations/database/annotation-hg37 \
    -assembly GRCh37 \
    -gerp19 /db/download/hg19.GERP_scores.tar.gz &

==================================================================
При упоковке по 100 вариантов:
139 объектов, всего 8,3 ГБ
===================================================================================
При упоковке по 200 вариантов:
135 объектов, всего 8,1 ГБ
********* conservation ***********
Колличество записей: 13241022

Использованные алгоритмы упаковки
SELECTIVE: Использовано: 0%, Средний размер пакета: 278 bytes
ORDERS: Использовано: 24%, Средний размер пакета: 801 bytes
ORDERS_WITH_DICTIONARY: Использовано: 73%, Средний размер пакета: 617 bytes
SELECTIVE_WITH_DICTIONARY: Использовано: 1%, Средний размер пакета: 259 bytes
===================================================================================
При упоковке по 200 вариантов:
114 объектов, всего 6,6 ГБ
********* conservation ***********
Колличество записей: 13241022
Чистый объем данных: 6412202922 bytes

Использованные алгоритмы упаковки
SELECTIVE_WITH_DICTIONARY: Использовано: 10589 (0,08%), Средний размер пакета: 31 bytes
ORDERS: Использовано: 1 (0,00%), Средний размер пакета: 801 bytes
SELECTIVE: Использовано: 5522 (0,04%), Средний размер пакета: 25 bytes
ORDERS_WITH_DICTIONARY_OVER_GZIP: Использовано: 1383388 (10,45%), Средний размер пакета: 86 bytes
ORDERS_OVER_GZIP: Использовано: 11841522 (89,43%), Средний размер пакета: 526 bytes