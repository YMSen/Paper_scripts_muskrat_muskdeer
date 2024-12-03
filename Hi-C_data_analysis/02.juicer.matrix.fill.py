import sys
import numpy as np
import gzip

inp=sys.argv[1]
fai=sys.argv[2]
chr=sys.argv[3]
res=int(sys.argv[4])
oup=sys.argv[5]

FAI={}
ff=open(fai,"r")
for line in ff:
	line=line.strip()
	lis=line.split("\t")
	FAI[lis[0]]=int(lis[1])

if inp.endswith(".gz"):
	MAT=gzip.open(inp,'r')
else:
	MAT=open(inp,'r')
line=MAT.readline()
line=line.strip()
lis=line.split("\t")
le=len(lis)
nl=FAI[chr]/int(res)+1
add=nl-le
MAT.close()

st=[]
ed=[]
for i in range(0,nl):
	st.append(i*res)
	ed.append((i+1)*res)
ed[-1]=FAI[chr]

OUP=open(oup,'w')
if add>0:
	cc1=np.repeat("0",add)
	ho="\t".join(cc1)+"\n"
	cc2=np.repeat("0",nl)
	hoho="\t".join(cc2)+"\n"
	#MAT=open(inp,'r')
	if inp.endswith(".gz"):
        	MAT=gzip.open(inp,'r')
	else:
        	MAT=open(inp,'r')
	ind=0
	for line in MAT:
		line=line.strip()
		line=line.replace("NaN","0")
		OUP.writelines(chr+"\t"+str(st[ind])+"\t"+str(ed[ind])+"\t"+line+"\t"+ho)
		ind+=1
	for i in range(add):
		OUP.writelines(chr+"\t"+str(st[ind])+"\t"+str(ed[ind])+"\t"+hoho)
		ind+=1
	MAT.close()
elif add==0:
	#MAT=open(inp,"r")
	if inp.endswith(".gz"):
        	MAT=gzip.open(inp,'r')
	else:
        	MAT=open(inp,'r')
	ind=0
	for line in MAT:
		line=line.strip()
		line=line.replace("NaN","0")
		OUP.writelines(chr+"\t"+str(st[ind])+"\t"+str(ed[ind])+"\t"+line+"\n")
		ind+=1
	MAT.close()
else:
	print("The juicer matrix is not in corresct format!!!!!!!!!")
