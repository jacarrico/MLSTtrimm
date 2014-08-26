import sys
from Bio import SeqIO

try:
        fn = sys.argv[1]
        thr = float(sys.argv[2])
        sz = int(sys.argv[3])
except IndexError:
    print "Usage: MLST_trimpattern.py <fasta_filename> <threshold (0-1)> <pattern size>"
    print "Reads a fasta file and outputs the trimming pattern sequences in IUPAC code by the selected size and threshold."
    print "jcarrico@fm.ul.pt, jmiranda@medicina.ulisboa.pt  30/05/2014"
    raise SystemExit



#fn = 'xpt_.fas'
#sz = 20

#thr=0.05



def NucToIUPAC(nuclist):
    IUPAC = {
        ('A',): 'A',
        ('A', 'C'): 'M',
        ('A', 'C', 'G'): 'V',
        ('A', 'C', 'G', 'T'): 'N',
        ('A', 'C', 'T'): 'H',
        ('A', 'G', 'T'): 'D',
        ('A', 'T'): 'W',
        ('C',): 'C',
        ('C', 'G', 'T'): 'B',
        ('G',): 'G',
        ('A', 'G'): 'R',
        ('C', 'G'): 'S',
        ('G', 'T'): 'K',
        ('T',): 'T',
        ('C', 'T'): 'Y',
        ('U',): 'U'}
    return IUPAC[nuclist] 

def patToTrimm(sz, pat, alelle_count, thr):
    
    TrimmPattern=''

    for i in range(sz):
        tmp=[]
        for j in ('A','C','G','T'):
            if float(pat[(i,j)])/alelle_count>thr:
                tmp.append(j)
        TrimmPattern+=NucToIUPAC(tuple(tmp))
    return TrimmPattern

def pat(fn, sz, pos):
    h = open(fn,'rU')
    pat={}
    
    for i in range(sz):
        for j in ('A','T','G','C'):
            pat[(i,j)]=0

    alelle_count=0

    for allele in SeqIO.parse(h, 'fasta'):
        alelle_count+=1
        if pos == "start":
            ss = allele.seq.tostring()[:sz]
        elif pos == "end":
            ss = allele.seq.tostring()[-sz:]
        ct=0
        for n in ss:
            pat[(ct, n)] += 1
            ct+=1
    h.close()
    return pat, alelle_count

pos = ["start", "end"]

pat_start, alelle_count = pat(fn, sz, pos[0])
pat_end, alelle_count = pat(fn, sz, pos[1])

TrimmPattern_start = patToTrimm(sz,pat_start,alelle_count,thr)
TrimmPattern_end = patToTrimm(sz,pat_end,alelle_count,thr)

print "TrimmPattern start:", TrimmPattern_start
print "TrimmPattern end:", TrimmPattern_end
