# dt2 :: deals with input datatype 2
# user defined parameters read from Control.txt

import copy
from math import pow, log, exp
from random import random
from random import randint
from numpy.random import beta
import scipy.stats 


inx=open("Control2.txt", "rU")
for line_idx, line in enumerate(inx):
        cols = line.split('\t') 
	u1=cols[0].split("=")
	if u1[0]=="InputDataType":
		InputDataType=int(u1[1])
		if InputDataType != 2:
			print "Using wrong program version for this input data type "
			break
	elif u1[0]=="FILEPREFIX":
		FILEPREFIX=u1[1]		# main input file is [FILEPREFIX].genotypes.txt 
	elif u1[0]=="rep":
		rep=u1[1]			# appended to front of outputs to distinguish this run of program (one chain)
	elif u1[0]=="NoSubpops":
		NoSubpops=int(u1[1])		# number of subpopulations (details in popfile)
	elif u1[0]=="NoFamilies":
		NoFamilies=int(u1[1])		# number of families
	elif u1[0]=="popfile":
		popfile=u1[1]			# file containing location of each maternal family
	elif u1[0]=="FullOutput":
		FullOutput=int(u1[1])		# set to 1 to output posterior probs for all genotypes
	elif u1[0]=="IHPriorModel":
		IHPriorModel=int(u1[1])		# 0=IH priors determined by t; 1= IH prior with uniform F
	elif u1[0]=="ChainLength":
		ChainLength=int(u1[1])		# run duration
	elif u1[0]=="afstep":
		afstep=float(u1[1])		# allele frequency step window size
	elif u1[0]=="tstep":
		tstep=float(u1[1])		# t step window size
	elif u1[0]=="BetaGstep":
		BetaGstep=float(u1[1])		# step size for BetaGu (relevant only for InputDataType=2)
	elif u1[0]=="AFboundary":
		AFboundary=float(u1[1])		# minimum for Minor Allele Frequency
	elif u1[0]=="thinningfreq":
		thinningfreq=int(u1[1])		# records taken from chain every thinningfreq steps
	elif u1[0]=="burnin":
		burnin=int(u1[1])		# length of burn-in
	elif u1[0]=="meanFst":
		meanFst=float(u1[1])		# Prior for degree of subpopulation differentiation
	elif u1[0]=="vfacFst":
		vfacFst=float(u1[1])		# variance factor for Fst Prior
		
inx.close()


FpopPrior = [meanFst,vfacFst*meanFst*(1-meanFst)] # mean, variance of beta prior
alphaFpop =FpopPrior[0]*(FpopPrior[0]*(1-FpopPrior[0])/FpopPrior[1] - 1)
betaFpop =(1-FpopPrior[0])*(FpopPrior[0]*(1-FpopPrior[0])/FpopPrior[1] - 1)
print "F pop prior ",alphaFpop, betaFpop


# Mendelism [mom genotype][dad allele][off]
OProb=[ [[1.0,0,0],[0,1.0,0]], [[0.5,0.5,0],[0,0.5,0.5]], [[0,1.0,0],[0,0,1.0]] ] # outcrossed progeny
SProb=[ [1.0,0,0], [0.25,0.5,0.25], [0,0,1.0] ] # selfed progeny
Fmom = [0,0.5,0.75,0.875,0.9375,0.96875,0.984375,1.0] # F given IH


# Genotyping probilities, "RA" values are updated for each snp 
GProb={"RR":{},"RA":{},"AA":{}}
GProb["RR"]["R"]=1.0 # [true genotype][observed data]
GProb["RR"]["H"]=0.0
GProb["RR"]["A"]=0.0
GProb["RR"]["N"]=1.0
GProb["RA"]["R"]=0.5
GProb["RA"]["H"]=0.5
GProb["RA"]["A"]=0.5
GProb["RA"]["N"]=1.0
GProb["AA"]["R"]=0.0
GProb["AA"]["H"]=0.0
GProb["AA"]["A"]=1.0
GProb["AA"]["N"]=1.0

########## Likelihood calculators


def fam1(Data,famID,snpID,delt): # family j probability for one SNP


	pop=SubPop[famID]
	q = Q[pop][snpID]
	f = Fmom[IH[famID]] 
	PriorMom={}
	PriorMom[0]= q*q*(1-f) + f*q
	PriorMom[1]= 2*q*(1-q)*(1-f)
	PriorMom[2]= 1.0-PriorMom[0]-PriorMom[1]

	dg=0
	Famsize=len(Data)
	# print j,Famsize
	# FruitNumber=PlantTypes[Famsize-1]
	for k in range(Famsize): # check for dataless fams
		if Data[k][0] != "N":
			dg = 1
	if dg==0:
		return 1.0 # no data for family

	else: # calculate mom prob
		nx=Data[0][1]-1
		if nx>0:
			BetaG = 1.0 - 0.5*exp(BetaGu[0]*float(CallCX[famID][0]-MaxCount))
			GProb["RA"]["H"]= 1.0 - BetaG**float(nx)
			GProb["RA"]["R"]= 0.5*(1.0-GProb["RA"]["H"])
			GProb["RA"]["A"]= 0.5*(1.0-GProb["RA"]["H"]) 
		else:
			GProb["RA"]["H"]= 0.0
			GProb["RA"]["R"]= 0.5
			GProb["RA"]["A"]= 0.5 


		if dg==1 and Famsize==1: # only a mom that has data (she produced no genotyped progeny)
			return PriorMom[0]*GProb["RR"][Data[0][0]] + PriorMom[1]*GProb["RA"][Data[0][0]] + PriorMom[2]*GProb["AA"][Data[0][0]] 

		elif dg==1 and Famsize>1: 

			FamLikelihood=0.0
			for MG in range(3):
				nx=Data[0][1]-1
				if nx>0:
					BetaG = 1.0 - 0.5*exp(BetaGu[0]*float(CallCX[famID][0]-MaxCount))
					GProb["RA"]["H"]= 1.0 - BetaG**float(nx)
					GProb["RA"]["R"]= 0.5*(1.0-GProb["RA"]["H"])
					GProb["RA"]["A"]= 0.5*(1.0-GProb["RA"]["H"]) 
				else:
					GProb["RA"]["H"]= 0.0
					GProb["RA"]["R"]= 0.5
					GProb["RA"]["A"]= 0.5 

				if MG==0:
					 Fk=GProb["RR"][Data[0][0]]
				elif MG==1:
					 Fk=GProb["RA"][Data[0][0]]
				elif MG==2:
					 Fk=GProb["AA"][Data[0][0]]

				# print "mom g ",MG,Fk

				for offspring in range(1,Famsize):
					nx=Data[offspring][1]-1
					if nx>0:
						BetaG = 1.0 - 0.5*exp(BetaGu[0]*float(CallCX[famID][offspring]-MaxCount))
						
						GProb["RA"]["H"]= 1.0 - BetaG**float(nx)
						GProb["RA"]["R"]= 0.5*(1.0-GProb["RA"]["H"])
						GProb["RA"]["A"]= 0.5*(1.0-GProb["RA"]["H"]) 
					else:
						GProb["RA"]["H"]= 0.0
						GProb["RA"]["R"]= 0.5
						GProb["RA"]["A"]= 0.5 

					
					if delt[offspring-1]==0:
						L0 = OProb[MG][0][0]*GProb["RR"][Data[offspring][0]] + OProb[MG][0][1]*GProb["RA"][Data[offspring][0]] + OProb[MG][0][2]*GProb["AA"][Data[offspring][0]] # dad gives R
						L1 = OProb[MG][1][0]*GProb["RR"][Data[offspring][0]] + OProb[MG][1][1]*GProb["RA"][Data[offspring][0]] + OProb[MG][1][2]*GProb["AA"][Data[offspring][0]] # dad gives A
						Fk*= (q*L0+(1-q)*L1)

						# print "zz ",L0,L1
					elif delt[offspring-1]==1:
						Fk*= ( SProb[MG][0]*GProb["RR"][Data[offspring][0]] + SProb[MG][1]*GProb["RA"][Data[offspring][0]] + SProb[MG][2]*GProb["AA"][Data[offspring][0]] )
							
				# print "given mom g ",MG,Fk
					
				FamLikelihood+=PriorMom[MG]*Fk
#		print "hx ",famID,FamLikelihood,Data,CallCX[famID]
		return FamLikelihood


def allfams1(snpID): # whole dataset probability for one SNP:: 

	ln_l = 0.0
	for j in range(1,NoFamilies+1):
		famL=fam1(FullData[snpID][j],j,snpID,delta[j])
		if famL>0.0:
			ln_l+=log(famL)
		else:
			print "allfams1, Impossible fam ",j,snpID,FullData[snpID][j]
			print "delta ",delta[j]
			return float('-inf')
	return ln_l


def allsnps1(famID): # fall snps for one family
	ln_l = 0.0
	for snpID in Q[0]:
		famL=fam1(FullData[snpID][famID],famID,snpID,delta[famID])
		if famL>0.0:
			ln_l+=log(famL)
		else:
			# print "allfsnps1, Impossible fam ",famID,snpID,FullData[snpID][famID]
			# print "delta ",delta[famID]
			return float('-inf')
	return ln_l


def FLL():
	BigLL=0.0 # LL of full data given all latents
	for k in range(NoSnps):
		BigLL+=allfams1(SNPList[k])
	return BigLL


########## UPDATES TO parameters

def UpdateBetaGu():

	take=[0,0]
	BetaC= BetaGu[0]
	BetaGp= BetaGu[0] + (random()-0.5)*BetaGstep
	if BetaGp < 0.1 and BetaGp > 0.0:
		cx=[0,0]

		LL0 =FLL()
		BetaGu[0]=BetaGp
		LL1 =FLL()

		if LL1>=LL0: 
			take[1]+=1
		else:
			if random()< exp(LL1-LL0):
				take[1]+=1
			else:
				take[0]+=1
				BetaGu[0]=BetaC # revert

	return take

def Updatet(): # currently uniform prior

	take=[0,0]
	tp= t[0] + (random()-0.5)*tstep
	if tp <= 1.0 and tp >= 0.0:
		cx=[0,0]
		for j in range(1,NoFamilies+1):
			for offspring in range(1,len(FullData[NamedSNP][j])):
				cx[ delta[j][offspring-1] ]+=1

		LL0 =float(cx[0])*log(t[0]) + float(cx[1])*log(1.0-t[0])
		LL1 =float(cx[0])*log(tp) + float(cx[1])*log(1.0-tp)
		if LL1>=LL0: 
			t[0]=tp
			take[1]+=1
		else:
			if random()< exp(LL1-LL0):
				t[0]=tp
				take[1]+=1
			else:
				take[0]+=1
	return take


def UpdateQA(): # currently uniform prior; this is "ancestral Q"
	take=[0,0]
	for snpID in Q[0]:
		q=Q[0][snpID]
		qp=Q[0][snpID]+ (random()-0.5)*afstep

		if qp > AFboundary and qp < 1.0-AFboundary:
			LL0=0.0
			LL1=0.0
			for spop in range(1,NoSubpops+1):
				LL0+=log( scipy.stats.beta.pdf(Q[spop][snpID],q*(1.0-Fpop[spop])/Fpop[spop],(1-q)*(1.0-Fpop[spop])/Fpop[spop]) )
				LL1+=log( scipy.stats.beta.pdf(Q[spop][snpID],qp*(1.0-Fpop[spop])/Fpop[spop],(1-qp)*(1.0-Fpop[spop])/Fpop[spop]) )

			if LL1>=LL0: 
				Q[0][snpID]=qp
				take[1]+=1
			else:
				if random()< exp(LL1-LL0):
					Q[0][snpID]=qp
					take[1]+=1
				else:
					take[0]+=1

	return take


def UpdateQsubpops(): # conditional on "ancestral Q"
	take=[0,0]
	for snpID in Q[0]:
		for spop in range(1,NoSubpops+1):
			q=Q[spop][snpID]
# 			qp= beta( Q[0][snpID]*(1.0-Fpop[spop])/Fpop[spop] , (1-Q[0][snpID])*(1.0-Fpop[spop])/Fpop[spop] ) # sample new value from the prior
			qp=Q[spop][snpID]+ (random()-0.5)*afstep # use prior ratio below
			if qp > AFboundary and qp < 1.0-AFboundary:
				PriorRatio= scipy.stats.beta.pdf(qp,Q[0][snpID]*(1.0-Fpop[spop])/Fpop[spop],(1-Q[0][snpID])*(1.0-Fpop[spop])/Fpop[spop])
				PriorRatio= PriorRatio/scipy.stats.beta.pdf(q,Q[0][snpID]*(1.0-Fpop[spop])/Fpop[spop],(1-Q[0][snpID])*(1.0-Fpop[spop])/Fpop[spop])
				LL0=0.0
				LL1=0.0

				for j in range(1,NoFamilies+1):
					if SubPop[j]==spop:
						LL0+=log( fam1(FullData[snpID][j],j,snpID,delta[j]) )
						Q[spop][snpID]=qp # temp change
						LL1+=log( fam1(FullData[snpID][j],j,snpID,delta[j]) )
						Q[spop][snpID]=q # revert 
				if LL1-LL0+log(PriorRatio)>0:
					Q[spop][snpID]=qp
					take[1]+=1
				elif random()< PriorRatio*exp(LL1-LL0):
					Q[spop][snpID]=qp
					take[1]+=1
				else:
					take[0]+=1

	return take

########## UPDATES TO LATENT VARIABLES

def UpdateFpop(): # 
	take=[0,0]
	for spop in range(1,NoSubpops+1):

		f=Fpop[spop]
		fp=Fpop[spop]+ 0.02*(random()-0.5) # use prior ratio below
		if fp>0.0001 and fp<(1-0.0001):
			PriorRatio= scipy.stats.beta.pdf(fp,alphaFpop, betaFpop)/scipy.stats.beta.pdf(f,alphaFpop, betaFpop)

			LL0=0.0
			LL1=0.0
			good=0
			for snpID in Q[0]:
				try:
					LL0+=log( scipy.stats.beta.pdf(Q[spop][snpID],Q[0][snpID]*(1.0-f)/f,  (1-Q[0][snpID])*(1.0-f)/f ) )
					LL1+=log( scipy.stats.beta.pdf(Q[spop][snpID],Q[0][snpID]*(1.0-fp)/fp, (1-Q[0][snpID])*(1.0-fp)/fp ) )
				except ValueError:
					take[0]+=1
					good=1
			if good==0:
				if LL1-LL0+log(PriorRatio)>0:
					Fpop[spop]=fp
					take[1]+=1
				else:
					if random() < PriorRatio*exp(LL1-LL0):
						Fpop[spop]=fp
						take[1]+=1
					else:
						take[0]+=1
	return take



def UpdateIH(): # sample proposed IH from prior (given current t)
	take=[0,0]
	cdf=[0.0 for k in range(8)]

	if IHPriorModel==1:
		cdf[0]=0.5
		for k in range(1,7):
			pk = 1.0/float(2.0+3.0*k+k*k)		
			cdf[k]=pk+cdf[k-1]

	else: # IH sampled based on current t
		cdf[0]=t[0]
		for k in range(1,7):
			cdf[k]= cdf[k-1] + t[0]*(1-t[0])**float(k)

	cdf[7]=1.0


	for j in range(1,NoFamilies+1):
		ih = IH[j] 
		val =random()
		for k in range(8):
			if val<cdf[k]:
				ihp = k
				break

		if ihp != IH[j]:

			LL0=allsnps1(j)
			IH[j]=ihp # temp change
			LL1=allsnps1(j)
			IH[j]=ih # revert

			if LL1-LL0>0:
				IH[j]=ihp
				take[1]+=1
			elif random()< exp(LL1-LL0):
				IH[j]=ihp
				take[1]+=1
			else:
				take[0]+=1

	return take


def UpdateDelta(): 
	take=[0,0]
	for j in range(1,NoFamilies+1):
		for offspring in range(1,len(FullData[NamedSNP][j])):
			deltp= copy.copy(delta[j])
			deltp[offspring-1]=1-delta[j][offspring-1]  # flat prior on t implies out/self equally likely (prior)

			PriorRatio=t[0]/(1.0-t[0]) # step to out
			if deltp[offspring-1]==1:
				PriorRatio=1.0/PriorRatio

			LL0=0.0
			LL1=0.0
			escape=0
			for snpID in Q[0]:
				L0=fam1(FullData[snpID][j],j,snpID,delta[j])
				L1=fam1(FullData[snpID][j],j,snpID,deltp)
				if L0==0.0:
					print "Should not be here, UpdateDelta()"
					break
				elif L1==0.0:
					take[0]+=1 # do not step
					escape = 1
					break
				else:
					LL0+=log(L0)
					LL1+=log(L1)

			if escape==0:

				if LL1-LL0+log(PriorRatio)>0: 
					delta[j][offspring-1]=deltp[offspring-1]
					take[1]+=1
				else:
					if random()< PriorRatio*exp(LL1-LL0):
						delta[j][offspring-1]=deltp[offspring-1]
						take[1]+=1
					else:
						take[0]+=1
	return take


def GenoPP(): # Posterior probabilities for genotypes
	for j in range(1,NoFamilies+1):
		famID=j
		for snpID in Q[0]:

			# Maternal post prob for snpID

			pop=SubPop[j]
			q = Q[pop][snpID]
			f = Fmom[IH[j]] 
			PriorMom={}
			PriorMom[0]= q*q*(1-f) + f*q
			PriorMom[1]= 2*q*(1-q)*(1-f)
			PriorMom[2]= 1.0-PriorMom[0]-PriorMom[1]
			PostMom={}

			Data=copy.copy(FullData[snpID][j])
			delt= copy.copy(delta[j])


			dg=0
			Famsize=len(Data)

			for k in range(Famsize): # check for dataless fams
				if Data[k][0] != "N":
					dg = 1
			if dg==0:
				pass # return 0 # no data for family

			elif dg==1 and Famsize==1: # only a mom that has data (she produced no genotyped progeny)
				nx=Data[0][1]-1
				BetaG = 1.0 - 0.5*exp(BetaGu[0]*float(CallCX[famID][0]-MaxCount))
				GProb["RA"]["H"]= 1.0 - BetaG**float(nx)
				GProb["RA"]["R"]= 0.5*(1.0-GProb["RA"]["H"])
				GProb["RA"]["A"]= 0.5*(1.0-GProb["RA"]["H"]) 				
				PostMom[0]=PriorMom[0]*GProb["RR"][Data[0][0]] 
				PostMom[1]=PriorMom[1]*GProb["RA"][Data[0][0]]
				PostMom[2]=PriorMom[2]*GProb["AA"][Data[0][0]]
				for ma in range(3):
					PostMom[ma]=PostMom[ma]/(PostMom[0]+PostMom[1]+PostMom[2])
					MatGPP[snpID][j][ma]+= PostMom[ma] 
					
					
			elif dg==1 and Famsize>1: 

				for MG in range(3):
					nx=Data[0][1]-1
					BetaG = 1.0 - 0.5*exp(BetaGu[0]*float(CallCX[famID][0]-MaxCount))
					GProb["RA"]["H"]= 1.0 - BetaG**float(nx)
					GProb["RA"]["R"]= 0.5*(1.0-GProb["RA"]["H"])
					GProb["RA"]["A"]= 0.5*(1.0-GProb["RA"]["H"]) 

					if MG==0:
						 Fk=GProb["RR"][Data[0][0]]
					elif MG==1:
						 Fk=GProb["RA"][Data[0][0]]
					elif MG==2:
						 Fk=GProb["AA"][Data[0][0]]

					for offspring in range(1,Famsize):
						nx=Data[offspring][1]-1
						BetaG = 1.0 - 0.5*exp(BetaGu[0]*float(CallCX[famID][offspring]-MaxCount))
						GProb["RA"]["H"]= 1.0 - BetaG**float(nx)
						GProb["RA"]["R"]= 0.5*(1.0-GProb["RA"]["H"])
						GProb["RA"]["A"]= 0.5*(1.0-GProb["RA"]["H"]) 

						if delt[offspring-1]==0:
							L0 = OProb[MG][0][0]*GProb["RR"][Data[offspring][0]] + OProb[MG][0][1]*GProb["RA"][Data[offspring][0]] + OProb[MG][0][2]*GProb["AA"][Data[offspring][0]] # dad gives R
							L1 = OProb[MG][1][0]*GProb["RR"][Data[offspring][0]] + OProb[MG][1][1]*GProb["RA"][Data[offspring][0]] + OProb[MG][1][2]*GProb["AA"][Data[offspring][0]] # dad gives A
							Fk*= (q*L0+(1-q)*L1)

						elif delt[offspring-1]==1:
							Fk*= ( SProb[MG][0]*GProb["RR"][Data[offspring][0]] + SProb[MG][1]*GProb["RA"][Data[offspring][0]] + SProb[MG][2]*GProb["AA"][Data[offspring][0]] )							
					PostMom[MG]=PriorMom[MG]*Fk
				for ma in range(3):
					PostMom[ma]=PostMom[ma]/(PostMom[0]+PostMom[1]+PostMom[2])
					MatGPP[snpID][j][ma]+= PostMom[ma]

				# paternal allele Post prob for each offspring
				for offspring in range(1,Famsize):
					if delt[offspring-1]==0:
						PostDad=[0.0,0.0]
						nx=Data[offspring][1]-1
						BetaG = 1.0 - 0.5*exp(BetaGu[0]*float(CallCX[famID][offspring]-MaxCount))
						GProb["RA"]["H"]= 1.0 - BetaG**float(nx)
						GProb["RA"]["R"]= 0.5*(1.0-GProb["RA"]["H"])
						GProb["RA"]["A"]= 0.5*(1.0-GProb["RA"]["H"]) 

						for MG in range(3):
							PostDad[0]+=PostMom[MG]*(OProb[MG][0][0]*GProb["RR"][Data[offspring][0]] + OProb[MG][0][1]*GProb["RA"][Data[offspring][0]] + OProb[MG][0][2]*GProb["AA"][Data[offspring][0]]) # dad gives R
							PostDad[1]+=PostMom[MG]*(OProb[MG][1][0]*GProb["RR"][Data[offspring][0]] + OProb[MG][1][1]*GProb["RA"][Data[offspring][0]] + OProb[MG][1][2]*GProb["AA"][Data[offspring][0]])

						PatGPP[snpID][j][offspring-1][0]+=(PostDad[0]*q)/(PostDad[0]*q+PostDad[1]*(1-q))
						PatGPP[snpID][j][offspring-1][1]+=(PostDad[1]*(1-q))/(PostDad[0]*q+PostDad[1]*(1-q))

	return 0

 
###########################



# Main program

src  =open(FILEPREFIX+".genotypes.txt", "rU") # FILEPREFIX determined in Control.txt
srx  =open("CX."+FILEPREFIX+".genotypes.txt", "rU") 

out1a =open(rep+"."+FILEPREFIX+".Q.MAP.txt", "w")
out2 =open(rep+"."+FILEPREFIX+".Chain.delta.txt", "w")
out2a =open(rep+"."+FILEPREFIX+".Delta.pp.txt", "w")
out3 =open(rep+"."+FILEPREFIX+".Chain.IH.txt", "w")
out3a =open(rep+"."+FILEPREFIX+".IH.pp.txt", "w")
out4 =open(rep+"."+FILEPREFIX+".Chain.t.txt", "w")
out4a =open(rep+"."+FILEPREFIX+".t.pp.txt", "w")
out5 =open(rep+"."+FILEPREFIX+".Chain.Fpop.txt", "w")
out5a =open(rep+"."+FILEPREFIX+".Fpop.pp.txt", "w")
out7 =open(rep+"."+FILEPREFIX+".step.history.txt", "w")


if FullOutput == 1:
	out8 =open(FILEPREFIX+"pat.allele.pp.txt", "w")
	out9 =open(FILEPREFIX+"mat.G.pp.txt", "w")
	MatGPP={}
	PatGPP={}


# set initial values for delta[j], IH of moms

IH={}
SubPop={}
delta={}
for j in range(1,NoFamilies+1):
	IH[j]=0	# everybody is outbred
	delta[j]=[]

inx=open(popfile, "rU")
for line_idx, line in enumerate(inx):
        cols = line.replace('\n', '').split('\t') 
	if line_idx>0:
		SubPop[int(cols[0])]=int(cols[1])

Q={}
Fpop = {}
for j in range(NoSubpops+1): # Q[0] is the "Ancestral Population" 
	Q[j]={}
	if j>0:
		Fpop[j]=meanFst # initial value for among subpop fst

NoSnps=0
FullData={}
CallCX={}

t ={0:0.9} # initial value for population outcrossing rate
BetaGu = {0:0.001} # prob that a het looks like a het = 1 - beta^(readdepth - 1)

MaxCount=0
callcounts=[]
for line_idx, line in enumerate(srx):
        cols = line.replace('\n', '').split('\t') 
# 0	0
# 1	1275
	callcounts.append(int(cols[1]))
	if int(cols[1])>MaxCount:
		MaxCount=int(cols[1])


SNPList=[]
vcx=0
for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 

	if line_idx==0:
		NamedSNP = cols[0]+"_"+cols[1]
# SNP	0	
# 1	par	N	0
# 1	off	H	12
	if cols[0][0]==NamedSNP[0]:
		NoSnps+=1
		snpID=cols[0]+"_"+cols[1]	
		FullData[snpID]={}
		if FullOutput == 1:
			MatGPP[snpID]={}
			PatGPP[snpID]={}
		SNPList.append(snpID)
		for j in range(NoSubpops+1):
			Q[j][snpID]=0.5

		for j in range(1,NoFamilies+1):
			FullData[snpID][j]=[] 
			if FullOutput == 1:
				MatGPP[snpID][j]=[0.0,0.0,0.0] 
	else:
		if cols[2]=='N':
			FullData[snpID][int(cols[0])].append(['N', 0]) 
		elif cols[3]=='0':
			FullData[snpID][int(cols[0])].append(['N', 0]) 
		elif cols[2]=='H' and int(cols[3])<2:
			FullData[snpID][int(cols[0])].append(['N', 0]) 
		else:
			FullData[snpID][int(cols[0])].append([cols[2], int(cols[3])]) # first list element is mom

		if snpID==NamedSNP:
			if cols[1]=="off":
				delta[int(cols[0])].append(0) # everybody is outcrossed at beginning of run
			try:
				uu1 = CallCX[int(cols[0])]
			except KeyError:
				CallCX[int(cols[0])]=[]

			CallCX[int(cols[0])].append(callcounts[vcx])
			vcx+=1

print "No Snps ",NoSnps
if FullOutput == 1:
	for snpID in Q[0]:
		for j in range(1,NoFamilies+1):
			Famsize=len(FullData[snpID][j])
			if Famsize>1:
				PatGPP[snpID][j]={}
				for offspring in range(1,Famsize):
					PatGPP[snpID][j][offspring-1]=[0.0,0.0]


ChainDelta={}
ChainIH={}
Chaint={}
ChainQ={}
ChainFpop={}
ChainBetaGu={}

for j in range(NoSubpops+1):
	ChainQ[j]={}
	if j>0:
		ChainFpop[j]={}
		for k in range(101):
			ChainFpop[j][k]=0
	for s in Q[0]: # cycle through all snps
		ChainQ[j][s]={}
		for k in range(101):
			ChainQ[j][s][k]=0

for j in range(101):
	Chaint[j]=0
	ChainBetaGu[j]=0

	
for j in range(1,NoFamilies+1):
	ChainIH[j]=[0 for x in range(8)]
	for offspring in range(1,len(FullData[NamedSNP][j])):
		ChainDelta[str(j)+"_"+str(offspring)]=[0,0]

movehistory=[[0,0] for zy in range(10)]
for steps in range(ChainLength):

	if steps % thinningfreq == 0:
		Chaint[int(t[0]*100+0.0049)]+=1		
		cLL=FLL()
                out4.write(str(steps)+'\t'+str(t[0])+'\t'+str(BetaGu[0])+'\t'+str(cLL)+'\n')

		for snpid in Q[0]:
			for j in range(NoSubpops+1):
				ChainQ[j][snpid][int(Q[j][snpid]*100+0.0049)]+=1
			
		for j in range(1,NoFamilies+1):
			out3.write(str(steps)+'\t'+str(j)+'\t'+str(IH[j])+'\n')
			if steps>=burnin:
				ChainIH[j][ IH[j] ]+=1
			for offspring in range(1,len(FullData[NamedSNP][j])):
				out2.write(str(steps)+'\t'+str(j)+'\t'+str(offspring)+'\t'+str(delta[j][offspring-1])+'\n')
				if steps>=burnin:
					ChainDelta[str(j)+"_"+str(offspring)][ delta[j][offspring-1] ]+=1
		out5.write(str(steps))
		for j in range(1,NoSubpops+1):
			out5.write('\t'+str(Fpop[j]))
			ChainFpop[j][int(Fpop[j]*100+0.0049)]+=1
		out5.write('\n')

		if FullOutput == 1:
			nothingnumber = GenoPP()

		print "Step ",steps," current LL ",cLL," t ",t[0]

#####	UPDATE parameters
	stay,go= UpdateFpop()
	movehistory[0][0]+=stay
	movehistory[0][1]+=go
	stay,go= Updatet()
	movehistory[1][0]+=stay
	movehistory[1][1]+=go
	stayD,goD = UpdateDelta()
	movehistory[2][0]+=stayD
	movehistory[2][1]+=goD
	stay,go= UpdateIH()
	movehistory[3][0]+=stay
	movehistory[3][1]+=go
	stay,go= UpdateQA()
	movehistory[4][0]+=stay
	movehistory[4][1]+=go
	stay,go= UpdateQsubpops()
	movehistory[5][0]+=stay
	movehistory[5][1]+=go
	stay,go= UpdateBetaGu()
	movehistory[6][0]+=stay
	movehistory[6][1]+=go


#output summaries
out2a.write('family\toffspring\tout\tself\n')
out3a.write('family\tIH=0\tIH=1\tIH=2\tIH=3\tIH=4\tIH=5\tIH=6\tIH=7+\n')
for j in range(1,NoFamilies+1):
	out3a.write(str(j))
	for z in range(8):
		out3a.write('\t'+str(ChainIH[j][z]))
	out3a.write('\n')
	for offspring in range(1,len(FullData[NamedSNP][j])):
		out2a.write(str(j)+'\t'+str(offspring)+'\t'+str(ChainDelta[str(j)+"_"+str(offspring)][0])+'\t'+str(ChainDelta[str(j)+"_"+str(offspring)][1])+'\n')

out4a.write('t value\tfrequency\n')
for j in range(101):
	out4a.write(str(float(j)/100 + 0.005)+'\t'+str(Chaint[j]) )
	out4a.write('\n')

out5a.write('Subpop\tvalue\tfrequency\n')
for k in range(1, NoSubpops+1):
	for j in range(101):
		out5a.write(str(k)+'\t'+str(float(j)/100 + 0.005)+'\t'+str(ChainFpop[k][j])+'\n')

out7.write("parameter	rejected	accepted	ARate\n")
stx=["Fpop","t","delta","IH","Qancestral","Qsubpop","BetaGu"]
for j in range(7):
	out7.write(stx[j]+'\t'+str(movehistory[j][0])+'\t'+str(movehistory[j][1])+'\t'+str(float(movehistory[j][1])/float(movehistory[j][1]+movehistory[j][0]))+'\n') 

for snpid in Q[0]:
	out1a.write(snpid)
	for j in range(NoSubpops+1):
		mx=0.0
		for k in range(101):
			if ChainQ[j][snpid][k]>mx:
				mx = ChainQ[j][snpid][k]
				val = k
		out1a.write('\t'+str(float(val)/100 + 0.005) )

	out1a.write('\n')


if FullOutput == 1:
	for snpID in Q[0]:
		for j in range(1,NoFamilies+1):
			rx=MatGPP[snpID][j][0]+MatGPP[snpID][j][1]+MatGPP[snpID][j][2]
			if rx>0.0:
				out9.write(snpID+'\t'+str(j)+'\t'+str(MatGPP[snpID][j][0]/rx)+'\t'+str(MatGPP[snpID][j][1]/rx)+'\t'+str(MatGPP[snpID][j][2]/rx)+'\n')

			for offspring in range(1,len(FullData[snpID][j])):
				rx=PatGPP[snpID][j][offspring-1][0]+PatGPP[snpID][j][offspring-1][1]
				if rx>0.0:
					out8.write(snpID+'\t'+str(j)+'\t'+str(offspring)+'\t'+str(PatGPP[snpID][j][offspring-1][0]/rx)+'\t'+str(PatGPP[snpID][j][offspring-1][1]/rx)+'\n')



