import sys
import random

# get parameters

if len(sys.argv) != 4:
    print ()
    print ("USAGE: python",sys.argv[0],"<size_sort> <size_test> <percentage>")
    print ("  - size_sort: number of reads in sort_reads.txt")
    print ("  - size_test: number of reads in test_reads.txt")
    print ("  - percentage: percentage of reads present in test_reads.tx")
    print ()
    sys.exit()

# get parameters

sort_size = int(sys.argv[1])
test_size = int(sys.argv[2])
percentage = int(sys.argv[3])

# check parameters

if percentage < 0 and percentage > 100:
    print ("% of test reads must be in range 0..100")
    sys.exit()

NT = ['A','C','G','T']

# generate sorted list of reads
LS = []
for i in range(sort_size):
    rd = ""
    for k in range(100):
        x = random.randint(0,3)
        rd += NT[x]
    LS.append(rd)
LS.sort()

# compute the number of existing reads
nb_exist = int((test_size * percentage)/100)
# set the indice in the list of test reads where the existing reads will be inserted
T = []
for i in range(test_size):
    T.append(0)
i = 0
while i < nb_exist:
    x = random.randint(0,test_size-1)
    if T[x] == 0:
        T[x] = 1
        i += 1

# generate list of test reads (and the list of existing reads)
# this is a mix of random reads and reads from the sorted list 
LT = []
LE = []
for i in range(test_size):
    if T[i] == 0:
        rd = ""
        for k in range(100):
            x = random.randint(0,3)
            rd += NT[x]
        LT.append(rd)
    else:
        j = random.randint(0,sort_size-1)
        LT.append(LS[j])
        LE.append(LS[j])
LE.sort()

# write results

# sorted reads
ff = open("sort_reads.txt","w")
for i in range(len(LS)):
    ff.write(LS[i]+"\n")
ff.close()

# test reads
ff = open("test_reads.txt","w")
for i in range(len(LT)):
    ff.write(LT[i]+"\n")
ff.close()

# existing read --> put in spy reads (will be used to check the output of TP2)
ff = open("spy_reads.txt","w")
for i in range(len(LE)):
    ff.write(LE[i]+"\n")
ff.close()



    
