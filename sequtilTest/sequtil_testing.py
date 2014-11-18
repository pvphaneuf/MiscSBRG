
# coding: utf-8

# In[1]:

import pysam


                samfile = pysam.Samfile("test.sam", mode="w", referencelengths=[100], referencenames=["testchr"])
read = pysam.AlignedRead()
read.seq = "AAAAAAAAAA"
read.pos = 2
                
# In[2]:

sequence = (
"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
"TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
"ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
"CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
"GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
"AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
"AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT"
"GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTT"
"GCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGC"
)


# In[3]:

len(sequence)


# In[4]:

def reverse_compliment(seq):
    matches = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(matches[i] for i in reversed(seq))


# In[5]:

with open("test_SE_plus.fastq", "w") as outfile:
    for i in range(10):
        outfile.write("@TEST_%d\n" % i)
        outfile.write(sequence[i + 20:i + 55] + "\n")
        outfile.write("+\n")
        outfile.write("~" * 35 + "\n")


# In[6]:

with open("test_SE_minus.fastq", "w") as outfile:
    rev = reverse_compliment(sequence)
    for i in range(10):
        outfile.write("@TEST_%d\n" % i)
        outfile.write(rev[i + 20:i + 55] + "\n")
        outfile.write("+\n")
        outfile.write("~" * 35 + "\n")


# In[7]:

get_ipython().system(u'bowtie2 "/home/aebrahim/sequencing/indexes/e_coli" test_SE_plus.fastq -S test_SE_plus.sam')
get_ipython().system(u'bowtie2 "/home/aebrahim/sequencing/indexes/e_coli" test_SE_minus.fastq -S test_SE_minus.sam')
get_ipython().system(u'bowtie2 -X 1000 "/home/aebrahim/sequencing/indexes/e_coli" -1 test_SE_plus.fastq -2 test_SE_minus.fastq -S test_PE.sam')
get_ipython().system(u'bowtie2 -X 1000 "/home/aebrahim/sequencing/indexes/e_coli" -2 test_SE_plus.fastq -1 test_SE_minus.fastq -S test_PE_rev.sam')


# In[8]:

from sequtil import makegff


# In[9]:

makegff.write_samfile_to_gff("test_SE_plus.sam", "test_SE_plus.gff")
makegff.write_samfile_to_gff("test_SE_plus.sam", "test_SE_plus_5.gff", five_prime=True)
makegff.write_samfile_to_gff("test_SE_minus.sam", "test_SE_minus.gff")
makegff.write_samfile_to_gff("test_SE_minus.sam", "test_SE_minus_5.gff", five_prime=True)
makegff.write_samfile_to_gff("test_PE.sam", "test_PE.gff")
makegff.write_samfile_to_gff("test_PE.sam", "test_PE_5.gff", five_prime=True)
makegff.write_samfile_to_gff("test_PE_rev.sam", "test_PE_rev.gff")
makegff.write_samfile_to_gff("test_PE_rev.sam", "test_PE_rev_5.gff", five_prime=True)


# In[10]:

get_ipython().system(u'cat test_SE_plus.gff')


# In[11]:

get_ipython().system(u'cat test_SE_minus.gff')


# In[12]:

get_ipython().system(u'cat test_PE.gff')


# In[13]:

get_ipython().system(u'cat test_PE_rev.gff')


# In[14]:

get_ipython().system(u'cat test_PE_5.gff')


# In[14]:



