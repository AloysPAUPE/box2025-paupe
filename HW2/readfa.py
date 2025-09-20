#! /usr/bin/env python3
import matplotlib.pyplot as plt
size=0
length=0

def readfq(fp):  # this is a generator function
    # From https://github.com/lh3/readfq/blob/master/readfq.py
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

def size_length(filename):
    #Returns the number of sequences and the cumulative length of the sequences
    #in filename
    f=open(filename)
    g=readfq(f)
    temp_size=0
    temp_length=0
    for elem in g:
        temp_size+=1
        temp_length+=len(elem[1])
    return (temp_size,temp_length)

def reverse(w):
    #Returns the reverse (mirror) of w
    n=len(w)
    res=""
    for i in range(1,n+1):
        l=w[n-i]
        if l=='A':
            res=res+'T'
        elif l=='T':
            res=res+'A'
        elif l=='C':
            res=res+'G'
        elif l=='G':
            res=res+'C'
    return res

def F (w,k):
    #Returns the set of unique canonical k-mers in w
    subwords_set=set()
    n=len(w)
    for i in range(n-k+1):
        temp=w[i:i+k]
        if reverse(temp) not in subwords_set:
            subwords_set.add(w[i:i+k])
    return (subwords_set)

def canonical_kmers(filename,k):
    #Returns the set of unique canonical k-mers in filename
    f=open(filename)
    g=readfq(f)
    kmers_set=set()
    for elem in g:
        kmers_set=kmers_set.union(F(elem[1],k))
    return kmers_set

def canonical_20mers(filename):
    #Returns the set of unique canonical 20-mers in filename
    return(canonical_kmers(filename,20))

def jaccard_index(s1,s2):
    #Returns the jaccard index of s1 and s2
    return (len(s1.intersection(s2))/len(s1.union(s2)))

def pairwise_jaccard_indices():
    #Computes and prints the pairwise jaccard indices of the fasta files
    canonical_sets=[canonical_20mers('file'+str(i)+'.fa') for i in range(1,7)]
    for i in range(6):
        s1=canonical_sets[i]
        for j in range(i+1,6):
            s2=canonical_sets[j]
            print("the jaccard index of "+str(i+1)+" and "+str(j+1)+" is "+str(jaccard_index(s1,s2)))
            
def file_lengths():
    #Prints the lengths of all fasta files
    for i in range(1,7):
        (l1,l2)=size_length('file'+str(i)+'.fa')
        print(l2)

def draw_hist(filename,k):
    #Draws the histogram of the distribution of i-mers for i<=k
    plt.hist(x=[i for i in range(1,k+1)],bins=k, weights=[len(canonical_kmers(filename,i)) for i in range(1,k+1)])
    plt.show()
