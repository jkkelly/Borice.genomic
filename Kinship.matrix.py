import sys
from math import pow, log


prefix ="Maternal"

out1 = open(prefix+".kinship.txt","w") 
out1a = open(prefix+".kinmat.txt","w") 

snps={}
loci=[]
genos={}
SNPs={}
in1 = open("AF_Moms.csv","rU") 
for line_idx, line in enumerate(in1):
        cols = line.replace('\n', '').split(',') 
# "","scaffold_1_4870","scaffold_1_20610",
# "1",NA,0.1960571102,0.1769567697,NA,0.925359322,0.8688017889,

	if line_idx==0:
		for k in range(1,len(cols)):
			loci.append(cols[k])
			SNPs[k]=[0,0.0,0.0]
		NoSNPs=len(cols)-1

	else:
		genos[cols[0]]=[]
		for k in range(1,len(cols)):
			if cols[k]=="NA":
				genos[cols[0]].append(-9)
			else:
				# print cols[k]
				SNPs[k][0]+=1
				SNPs[k][1]+=float(cols[k])
				SNPs[k][2]+=( float(cols[k])*float(cols[k]) )
				genos[cols[0]].append(float(cols[k]))

plants = genos.keys()
print "Number of plants ",len(plants)

for k in SNPs:
	SNPs[k][1]=SNPs[k][1]/float(SNPs[k][0]) # mean score
	SNPs[k][2]=SNPs[k][2]/float(SNPs[k][0]) - (SNPs[k][1]*SNPs[k][1])


Relatedness={}
for j in range(len(plants)):
	Relatedness[plants[j]]={}
	for k in range(j,len(plants)):
		Relatedness[plants[j]][plants[k]]=[0,0.0]
		for x in SNPs:
			if genos[plants[j]][x-1] >=0 and genos[plants[k]][x-1] >=0:
				Relatedness[plants[j]][plants[k]][0]+=1
				Relatedness[plants[j]][plants[k]][1]+=( (genos[plants[j]][x-1]-SNPs[x][1])*(genos[plants[k]][x-1]-SNPs[x][1])/SNPs[x][2] )		

		if j==k:
			print plants[j],plants[k],Relatedness[plants[j]][plants[k]][0],Relatedness[plants[j]][plants[k]][1]/Relatedness[plants[j]][plants[k]][0]	
	
		out1.write(plants[j]+'\t'+plants[k]+'\t'+str(Relatedness[plants[j]][plants[k]][0])+'\t'+str(Relatedness[plants[j]][plants[k]][1]/Relatedness[plants[j]][plants[k]][0])+'\n')

	out1a.write(plants[j])
	for k in range(len(plants)):
		if k<j:
			out1a.write('\t'+str(Relatedness[plants[k]][plants[j]][1]/Relatedness[plants[k]][plants[j]][0]))
		else:
			out1a.write('\t'+str(Relatedness[plants[j]][plants[k]][1]/Relatedness[plants[j]][plants[k]][0]))

	out1a.write('\n')

