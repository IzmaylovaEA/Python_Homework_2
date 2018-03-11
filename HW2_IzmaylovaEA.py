
# coding: utf-8

# In[ ]:

from Bio import SeqIO
import re
import operator
from tabulate import tabulate

class Kmer:
    counter = 1
    
    def __init__(self,kmer_name):
        self.sequence = kmer_name
        self.position = []
    
    def increase(self):
        self.counter += 1
        
    # Функция для определения позиции k-мера.
    def coordinates(self,name,numbers):
        self.position.append([name,numbers])
    
    def __repr__(self):
        return str(self.counter)+', '+str(self.position)

# Функция для нахождения и вычисления координат наиболее частого k-мера.
def frequent_kmer(fasta, kmer_size):
    if kmer_size < 1 or type(kmer_size) == float:
        raise Exception('The second argument must be an integer greater than 0')
    else:
        kmer_dict = {}
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            # Название хромосомы.
            find = re.compile(r"^[^.]*")
            chrom_name = find.findall(record.name)
            seq = str(record.seq).upper()
            seq_lng = len(seq)
            for index in range(seq_lng-kmer_size+1):
                current_kmer = seq[index:(index+kmer_size)]
                if current_kmer in kmer_dict:
                    kmer_dict[current_kmer].increase()
                    if kmer_size == 1:
                        kmer_dict[current_kmer].coordinates(chrom_name,str(index+1))
                    else:
                        kmer_dict[current_kmer].coordinates(chrom_name,str(index+1)+'-'+str(index+kmer_size))
                else:
                    kmer_dict[current_kmer] = Kmer(current_kmer)
                    if kmer_size == 1:
                        kmer_dict[current_kmer].coordinates(chrom_name,str(index+1))
                    else:
                        kmer_dict[current_kmer].coordinates(chrom_name,str(index+1)+'-'+str(index+kmer_size))

    sorted_by_counter = sorted(kmer_dict.values(), key = operator.attrgetter('counter'), reverse = True)[0]
    print('The most common k-mer is:', sorted_by_counter.sequence,'\n')
    print(tabulate(sorted_by_counter.position, headers = ['Chromosome','Coordinates']))

frequent_kmer('/home/izmaylova/Downloads/seq_y_pestis.fasta', 23)

