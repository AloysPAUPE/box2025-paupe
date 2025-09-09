import string
import matplotlib.pyplot as plt

def Fibonacci (k):
    #Returns the k-th Fibonacci word
    if k==1:
        return "A"
    w=Fibonacci(k-1)
    n=len(w)
    strings=["A"]*n
    for i in range(n):
        if w[i]=='A':
            strings[i]="AB"
    return "".join(strings)

d_1="".join(open("dna_1.txt",'r').read().split('\n'))
#The DNA sequence from the genome
d_1_length=121393

import numpy as np
    
def F (w,k):
    #Returns the amount of unique subwords of length k in w
    subwords_set=set()
    n=len(w)
    for i in range(n-k+1):
        subwords_set.add(w[i:i+k])
    return len(subwords_set)

def random_dna_string (n):
    #Returns the equivalent of a random DNA string of length n
    return "".join([str(np.random.randint(4)) for k in range(n)])

def display ():
    x_points = [k for k in range(1,31)]
    y_points = [F(d_1,k) for k in range(1,31)]
    plt.subplot(1,3,1)
    plt.plot(x_points,y_points)
    plt.xlabel("Length of subwords - k")
    plt.ylabel("Amount of unique subwords - F(k)")
    #Plots F(k) for the genome
    d_2=random_dna_string(len(d_1))
    y_points = [F(d_1,k) for k in range(1,31)]
    plt.subplot(1,3,2)
    plt.plot(x_points,y_points)
    plt.xlabel("Length of subwords - k")
    plt.ylabel("Amount of unique subwords - F(k)")
    #Plots F(k) for a random DNA string
    w=Fibonacci(26)
    #len(Fibonacci(26))=d_1_length
    y_points = [F(w,k) for k in range(1,31)]
    plt.subplot(1,3,3)
    plt.plot(x_points,y_points)
    plt.xlabel("Length of subwords - k")
    plt.ylabel("Amount of unique subwords - F(k)")
    #Plots F(k) for the Fibonacci word of length d_1_length
    plt.tight_layout()
    plt.show()

#In the end, the DNA string from the genome behaves similarly to a randomly
#generated DNA string, whereas for the fibonacci word, F(k)=k+1

