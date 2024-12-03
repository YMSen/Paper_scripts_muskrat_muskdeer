import sys
with open(sys.argv[1],'r') as r:
    b=r.readlines()

for j in range(len(b)):
    jj=b[j].strip().split()
    print jj[0]+' '+str(j+1)+' '+jj[1]
