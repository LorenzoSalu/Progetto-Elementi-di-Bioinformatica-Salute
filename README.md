# Progetto-Elementi-di-Bioinformatica-Salute
Progetto del corso "Elementi di Bioinformatica" svolto da Lorenzo Salute (894461) durante l'a.a. 2024/25.

Dopo la compilazione del file "runner.cpp" è possibile eseguire il file eseguibile ottenuto per l'esecuzione automatica dello script python e l'ottenimento del report testuale.

Per garantire la corretta esecuzione dello script è necessaria l'installazione delle librerie pysam, Bio, HTSeq e subprocess.

------------------------------------

Tema di progetto:

Dato un file in formato GTF che fornisce i trascritti di geni umani su un certo cromosoma (in cui sia presente un solo trascritto per gene), un file in formato BAM che fornisce gli allineamenti di reads trascrittomici sequenziati da uno dei trascritti presenti nel file GTF e il file FASTA del cromosoma (da scaricare da Ensemble Genome Browser), produrre uno script/notebook che:

1) individui il trascritto da cui sono stati sequenziati i reads e il gene a cui appartiene tale trascritto;
2) produca un report che mostri:
  - il supporto di ogni introne del trascritto individuato (numero di reads che si allineano all’introne);
  - per ognuno dei reads allineati all’introne con il supporto maggiore, una visualizzazione testuale del loro allineamento al cromosoma

NOTA BENE: i reads spliced nel file BAM potrebbero non allinearsi esattamente ai siti di splicing dell’introne supportato nel file GTF. Occorre quindi prevedere una tolleranza di qualche base per fare corripondere un introne derivato da cigar string a un introne derivato dal file GTF.
