#!/usr/bin/env python
# coding: utf-8

# ### Importazione librerie utili

# In[2]:


#Importazione delle librerie utili
import pysam
from pysam import AlignmentFile
import pandas as pd
import numpy as np
import re
import Bio
from Bio import Align
from collections import Counter
import itertools
import subprocess


# ### Funzioni

# In[4]:


# Funzione per la pulizia e aggiustamento dei dati ricavati dal file gtf
def gtf_to_dataframe(gtf_file_location) :
    # Lettura file gtf e salvataggio dataframe
    df = pd.read_csv(gtf_file_location, sep = '\t', header = None)

    # Rimpiazzo dei nomi delle colonne
    replace_dict = {0 : 'reference', 1 : 'source', 2 : 'feature', 3 : 'start',
                4 : 'end', 5 : 'score', 6 : 'strand', 7 : 'frame', 8 : 'attributes'}
    df.rename(columns = replace_dict, inplace = True)
    
    # Rimozione colonne non utili
    df.drop(['source', 'score'], axis = 1, inplace = True)

    # Creazione colonne di id trascritto e nome gene + rimozione colonna attributi
    df['transcript_id'] = df['attributes'].apply(lambda x : re.search(r'transcript_id\s+"(.+?)";', x).group(1))
    df['gene'] = df['attributes'].apply(lambda x : re.search(r'gene_id\s+"(.+?)";', x).group(1))
    df.drop('attributes', axis = 1, inplace = True)

    # Creazione colonna che contiene la lunghezza
    df['length'] = df['end'] - df['start'] + 1

    # Reindicizzazione delle colonne del df
    df = df.reindex(columns = ['reference', 'feature', 'start', 'end',
                      'length', 'strand', 'frame', 'gene', 'transcript_id'])
    return df


# In[5]:


# Funzione per creare il file contenente i reads spliced
def get_spliced_reads_file(all_alignments, output_file):

    # Apro il file e lo svuoto ("reset")
    with open(output_file, "w") as f:
        pass 
    
    # Recuperare i reads spliced cercando la N nelle cigar string
    all_spliced_reads = list(set([alignment for alignment in all_alignments if 'N' in alignment.cigarstring]))

    # Costruzione del file
    with open(output_file, "w") as fastq_file:
        for read in all_spliced_reads:
            read_id = read.query_name
            read_sequence = read.query_sequence

            if read.query_qualities == None:
                # Viene salvata una quality string di default 
                read_quality_string = "I" * len(read_sequence)
            else:
                read_quality_string = read.query_qualities

            # Scrittura sul file
            fastq_file.write(f'@'+read_id+'\n'+read_sequence+'\n'+'+\n'+read_quality_string+'\n')

    # Ritorno la lista di tutti i read spliced
    return all_spliced_reads


# In[6]:


def get_transcripts(df):
    # Viene quindi trovato lo start e l'end per ogni trascritto e il suo gene di appartenenza
    transcript_dict = {}
    
    # Scorro sui geni
    for gene in set(df['gene']):
        
        # Maschera per prendere soltanto trascritti associati al gene e che sono esoni
        mask = (df.gene == gene) & (df.feature == 'exon')
        
        # Scorro sui trascritti per il gene
        for transcript_id in df[mask].groupby('transcript_id').groups:
    
            # Trovo tutti gli esoni associati al trascritto
            index_list = df[mask].groupby('transcript_id').groups[transcript_id]
            exon_list_sorted = sorted(zip(df.loc[index_list]['start'], df.loc[index_list]['end']))
            
            # Salvo la entry trascritto : (gene, start_trascritto, end_trascritto) in trascript_list
            transcript_dict[transcript_id] = (gene, exon_list_sorted[0][0], exon_list_sorted[-1][1])

    return transcript_dict


# In[7]:


def get_associated_transcripts(transcript_dict, all_alignments):
    # Dizionario read : trascritto da cui è stato sequenziato
    associated_transcript_dict = {}
    
    for alignment in all_alignments:
        for transcript_key in transcript_dict:
            if ((transcript_dict[transcript_key][1] <= alignment.reference_start <= transcript_dict[transcript_key][2]) &
                (transcript_dict[transcript_key][1] <= alignment.reference_end <= transcript_dict[transcript_key][2])):
                associated_transcript_dict[alignment.query_name] = (transcript_key, transcript_dict[transcript_key][0])

    return associated_transcript_dict


# In[8]:


def get_sequentied_transcript(associated_transcript_dict):
    # Risulta un unico elemento cioè il trascritto da cui sono stati sequenziati tutti i read
    from_transcript = list(set(associated_transcript_dict.values()))[0]

    return from_transcript[0]


# In[9]:


def get_gene_from_transcript(sequentied_transcript, transcript_dict):
    return transcript_dict[sequentied_transcript][0]


# In[10]:


def get_sequentied_transcript_subprocess(bam_file_location, gtf_file_location, output_file):
    # Comando unico per trovare da quali trascritti sono stati sequenziati i read
    cmd = [
        "htseq-count", "-f", "bam", "-r", "pos", "-s", "no",
        "-i", "transcript_id",
        bam_file_location, gtf_file_location
    ]

    # Scrittura file associazioni trascritto : numero read sequenziati da esso
    with open(output_file, "w") as out_file:
        subprocess.run(cmd, stdout=out_file)

    # Creo il dizionario per verificare i risultati del comando unico
    read_transcript_association_dict = {}

    # Leggo il file generato dal comando
    with open(output_file, "r") as file:
        for line in file:
            key, value = line.strip().split("\t")  # Divide sulla tabulazione
            read_transcript_association_dict[key] = int(value)  # Converte il valore in intero
    
    # Filtrare solo i record con value > 0
    read_transcript_association_dict = {key: value for key, value in read_transcript_association_dict.items() if value > 0}
    
    return read_transcript_association_dict


# In[11]:


def get_sequentied_transcript_file(sequentied_transcript, transcript_associated_gene, output_file):
    with open(output_file, "w") :
            pass
    
    # Scrive ogni coppia (chiave, valore) su una riga del file
    with open(output_file, "w") as f:
        f.write(f'SEQUENTIED TRANSCRIPT: {sequentied_transcript}\n')
        f.write(f'ASSOCIATED GENE: {transcript_associated_gene}')


# In[12]:


def get_best_intron(all_introns):
    # Trova il valore di copertura massimo
    max_cover = max(all_introns.values())
    
    # Trova tutti gli introni con copertura massima
    max_introns = [(interval, value) for (interval, value) in all_introns.items() if value == max_cover] 

    # Prendo come introne da "studiare" il primo della lista
    best_intron = max_introns[0]

    return best_intron


# In[13]:


def get_introns_support_file(all_introns, best_intron, output_file):
    with open(output_file, "w") :
            pass
    
    # Scrive ogni coppia (chiave, valore) su una riga del file
    with open(output_file, "w") as f:
        f.write(f'INTRONS: \n\n')
        
        for locus, count in all_introns.items():
            f.write(f'- {locus} : {count}\n')

        f.write(f'\nMax support intron -> {best_intron[0]} with support = {best_intron[1]} ')


# In[14]:


def get_alignment_file(intron_aligned_reads, reference, output_file):
    with open(output_file, "w") :
        pass
    
    with open(output_file, "w") as f:

        f.write(f'READS ALIGNED TO BEST INTRON: \n\n')
        
        for read in intron_aligned_reads:
            # Salvo info del read
            read_seq = read.seq
            read_cigar = read.cigarstring
            read_blocks = read.get_blocks() # Salva lista di (start, end) delle sottostringhe allineate e senza gap del read 
        
            # Lista di ogni carattere del read con: posizione sulla reference e operatore associato
            aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True, with_cigar=True)
        
            f.write(f'{read.query_name} -------------------------------------------------\n\n')
            f.write(f'READ INFO:\n\n')
            f.write(f'Read sequence -> {read_seq}\n')
            f.write(f'Cigar string -> {read_cigar}\n\n\n')
            f.write(f'ALIGNMENT: \n\n')
            
            # Scorro sui blocchi del read in modo da generare un allineamento testuale per ciascuno di essi
            for i, (block_start, block_end) in enumerate(read_blocks):
                query = ''
            
                # Scorro sulle aligned_pairs per costruire la query
                for (_, ref_position, char, _) in aligned_pairs:
                    # pair[1] contiene la posizione sulla reference, controllo che non sia None 
                    if (ref_position != None):
            
                        # Se il carattere fa parte del blocco considerato viene aggiunto alla query
                        if ((ref_position >= block_start) & (ref_position < block_end)):
                            query = query + char
            
                # Rendo la query uppercase dato che ho notato che ci sono dei caratteri minuscoli
                query = query.upper()
            
                # Recupero la reference relativa alle posizioni del blocco
                reference_seq = reference[block_start : block_end]
        
                # Eseguo l'allineamento
                result_alignments = aligner.align(reference_seq, query)
        
                # Recupero l'allineamento con score maggiore
                best_alignment = sorted(result_alignments)[0]
        
                # Stampo allineamento
                f.write(f'{best_alignment}\n')
        
                if i != len(read_blocks) - 1:
                    f.write(f'************************ ... INTRON ... ************************\n\n')
        
            f.write(f'\n\n\n')


# In[15]:


def write_report(sequentied_transcript_file, introns_file, alignment_file, report_file):
    with open(report_file, "w") as out_f:
        for file in [sequentied_transcript_file, introns_file, alignment_file]:
            with open(file, "r") as in_f:
                out_f.write(in_f.read() + "\n\n\n###################################################################################\n\n\n") 


# ### Parametri di input

# In[17]:


gtf_file_location = './input_files/annotation_one_tr_chr8.gtf'
bam_file_location = './input_files/sample-chr8.bam'
fasta_file_location = './input_files/Homo_sapiens.GRCh38.dna.chromosome.8.fa'


# ### File di output

# In[19]:


spliced_read_file = "./output_files/spliced_reads.fastq"
output_counts_transcript = "./output_files/output_counts_transcript.txt"
sequentied_transcript_file = "./output_files/sequentied_transcript_file.txt"
introns_file = "./output_files/introns_file.txt"
alignment_file = "./output_files/alignment_output.txt"
report_file = "./output_files/report.txt"


# ### Codice 

# In[21]:


df = gtf_to_dataframe(gtf_file_location)
#df


# In[22]:


#lettura del file BAM
pysam.index(bam_file_location)
bam_file = AlignmentFile(bam_file_location, 'rb')

#Lettura degli allineamenti
all_alignments = list(bam_file.fetch())


# In[23]:


# Creazione del file contenente i reads spliced e ritorna lista dei reads spliced
spliced_reads = get_spliced_reads_file(all_alignments, spliced_read_file)
#spliced_reads


# In[24]:


# Viene quindi trovato lo start e l'end per ogni trascritto e il suo gene di appartenenza
transcript_dict = get_transcripts(df)
#transcript_dict


# In[25]:


# Dizionario read : trascritto da cui è stato sequenziato
associated_transcript_dict = get_associated_transcripts(transcript_dict, all_alignments)
#associated_transcript_dict


# In[26]:


# Salvo il trascritto da cui sono stati sequenziati i read
sequentied_transcript = get_sequentied_transcript(associated_transcript_dict)
#sequentied_transcript


# In[27]:


# Salvo il gene associato al trascritto da cui sono stati sequenziati i read
transcript_associated_gene = get_gene_from_transcript(sequentied_transcript, transcript_dict)
#transcript_associated_gene


# In[28]:


read_transcript_association_dict = get_sequentied_transcript_subprocess(bam_file_location, gtf_file_location, output_counts_transcript)
#read_transcript_association_dict


# In[29]:


# Cerco gli introni del trascritto e numero di reads che supportano l'introne
all_introns = bam_file.find_introns(bam_file.fetch())
#all_introns


# In[30]:


best_intron = get_best_intron(all_introns)
#best_intron


# In[31]:


intron_start, intron_end = best_intron[0][0], best_intron[0][1]


# In[32]:


# Lettura del file FASTA indicizzato
fasta = pysam.FastaFile(fasta_file_location)

# Specifica il cromosoma e le coordinate (0-based)
chrom = fasta.references[0]

# Estrarre la sequenza
reference = fasta.fetch(chrom)


# In[33]:


# Salvo le posizioni di inizio e fine del trascritto sequenziato
transcript_position = (transcript_dict[sequentied_transcript][1], transcript_dict[sequentied_transcript][2])
transcript_start, transcript_end = transcript_position[0], transcript_position[1]
#print(transcript_start, transcript_end)


# In[34]:


# Cerco i reads allineati all'introne
intron_aligned_reads = [read for read in spliced_reads if ((read.reference_start <= intron_start) & (read.reference_end >= intron_end))]
#intron_aligned_reads


# In[35]:


aligner = Align.PairwiseAligner()


# In[36]:


get_sequentied_transcript_file(sequentied_transcript, transcript_associated_gene, sequentied_transcript_file)
get_introns_support_file(all_introns, best_intron, introns_file)
get_alignment_file(intron_aligned_reads, reference, alignment_file)

# Scrittura del report 
write_report(sequentied_transcript_file, introns_file, alignment_file, report_file)

