import pandas as pd
import Bio as bio
from Bio import SeqIO
import glob

# generates list of fastq files to analyze
sources = glob.glob('*.fastq')
# reads the fastq files into a dictionary with the file names as keys
fastqdict = {}
for i in range(len(sources)):
    temp = list(SeqIO.parse(sources[i],"fastq"))
    fastqdict[sources[i]]= [str(temp[k].seq) for k in range(len(temp))]

# the referenced sequence to be searched for is entered into the following dictionary with
# #an appropriate key
scaffdict = {'HEK3':'CAGAGGACCGACTCGGTCCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTT GCTATTTCTAGCTCTAAAACTCACGTGCTCAGTCTGGGCCGGTG',  'EMX1':'ATCACGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACT TGCTATTTCTAGCTCTAAAACTTCTTCTTCTGCTCGGACTCGGTG', 'FANCF':'TTTCCGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAAC TTGCTATTTCTAGCTCTAAAACGGTGCTGCAGAAGGGATTCCGGTG', 'RNF2':'TCGTTGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTT GCTATTTCTAGCTCTAAAACCAGGTAATGACTAAGATGACGGTG'}
# matches and counts iterative slices of the reference string to the appropriate fastq files
# reference key must be contained in the name of the fastq file
# generated values represent cumulative counts for a minimum degree of sgRNA integration
# i.e. a given value x means x reads contain y or more bases of the scaffold
resultdict = dict.fromkeys(sources)
for key in fastqdict:
    for scaffold in scaffdict:
        if scaffold in str(key):
            resultlist = []
            for j in range(len(scaffdict[scaffold])):
                extent = scaffdict[scaffold][0:(j + 1)]
                counter = 0
            for i in range(len(fastqdict[key])):
                if extent in fastqdict[key][i]:
                    counter = counter + 1
            resultlist.append(counter)
        resultdict[key] = resultlist

# writes the results into a dataframe indexed from 1
resultdf = pd.DataFrame.from_dict(resultdict)
resultdf = resultdf.reindex(sorted(resultdf.columns), axis=1)
resultdf.index = range(1,len(resultdf)+1)

# converts the cumulative count values into specific counts
# i.e. a given value x means x reads contain exactly y bases of the scaffold
resultdf2=resultdf.copy()
for entry in resultdf:
    for i in range(1,len(resultdf[entry])+1):
        try:
            resultdf2[entry][i] = resultdf[entry][i]-resultdf[entry][i+1]
        except:
            resultdf2[entry][i] = resultdf[entry][i]

# converts the specific counts values into frequencies
resultdf3=resultdf2.copy()
for entry in resultdf3:
    resultdf3[entry]=resultdf2[entry].div(resultdf[entry][1])*100
# reads the results into excel files
resultdf.to_excel('cumulativecounts.xlsx')
resultdf2.to_excel('specificcounts.xlsx')
resultdf3.to_excel('specificfrequencies.xlsx')


