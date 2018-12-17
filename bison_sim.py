import msprime
import numpy as np
import random

#functions

def print_list(sim):
	for st in sim:
		print (st)

def het_calc(chrom1,chrom2,length):
	count=0
	for pos in range(0,length-1):
		if chrom1[pos]!=chrom2[pos]:
			count=count+1
	return count


#setting parameters

Ne=10
gen_time=8
T=20
mu=10e-3
rho=0
sim_length=5

population_configurations = [
	msprime.PopulationConfiguration(
		initial_size=Ne),
	msprime.PopulationConfiguration(
		initial_size=Ne),
	msprime.PopulationConfiguration(
		initial_size=Ne),
	msprime.PopulationConfiguration(
		initial_size=Ne),
	msprime.PopulationConfiguration(
		initial_size=Ne),
	]

samples = [
	msprime.Sample(population=0,time=0),
	msprime.Sample(population=0,time=0),

	msprime.Sample(population=1,time=0),
	msprime.Sample(population=1,time=0),

	msprime.Sample(population=2,time=0),
	msprime.Sample(population=2,time=0),

	msprime.Sample(population=3,time=0),
	msprime.Sample(population=3,time=0),

	msprime.Sample(population=4,time=0),
	msprime.Sample(population=4,time=0),
]

demographic_events=[
	msprime.MassMigration(
		time=T,source=0,destination=4,proportion=1),
	msprime.MassMigration(
		time=T,source=1,destination=4,proportion=1),
	msprime.MassMigration(
		time=T,source=2,destination=4,proportion=1),
	msprime.MassMigration(
		time=T,source=3,destination=4,proportion=1),
]

def bison_sim():
	tree_sequence=msprime.simulate(Ne=Ne,length=sim_length,recombination_rate=rho,mutation_rate=mu,population_configurations=population_configurations,
		demographic_events=demographic_events,samples=samples)
#	print(tree_sequence.genotype_matrix())
	haplotypes=tree_sequence.haplotypes()
	haplotype_list=[]
	for i in haplotypes:
				haplotype_list.append(i)
#	haplotype_list.append(tree_sequence.get_num_mutations())
	polymorph=tree_sequence.get_num_mutations()
	sim=haplotype_list
	print_list(sim)
	print(" ")
#	print(sim[0][0])
#	print(sim[1])
	print(polymorph)
	print(het_calc(sim[0],sim[1],length=polymorph)/sim_length)
	print(" ")




bison_sim()


########

def abc():
	#for seed in range (1,10000):
	count=0
	for _ in range (0,1000000):
		seed=random.randint(1,(2**32)-1)
		#t2=random.uniform(11700,13700)
		t2=12800
		t1=random.uniform(903,t2)
	#	print(seed,t1,t2)
		tree_sequence=msprime.simulate(Ne=N_DNG,length=1e7,recombination_rate=0.33e-8,samples=ancient_samples,population_configurations=population_configurations,demographic_events=double_split_time(t1,t2),
            	mutation_rate=0.66e-8, random_seed=seed)
		haplotypes=tree_sequence.haplotypes()
		haplotype_list=[]
		for i in haplotypes:
			haplotype_list.append(i)
		haplotype_list.append(tree_sequence.get_num_mutations())
		sim=haplotype_list
	#	print_list(sim)
		het_list=het_B_finder(sim)
		derived_list=derived_finder(sim,het_list)
	#	print("het list is",het_list)
		fab_k=fab_calc(sim,derived_list)
	#	print(fab_k,d_calc(fab_k,fab_obs),t1,t2,len(derived_list[1]),seed)
		count=count+1
		print(fab_k,t1,t2,len(derived_list[1]),seed,count)


def test_abc():
	t1=random.uniform(1,14900)
	t2=random.uniform(t1,14900)
	print(t1,t2)
	tree_sequence=msprime.simulate(Ne=N_DNG,length=1000,recombination_rate=5e-7,population_configurations=population_configurations,demographic_events=split_times(t1),
            	mutation_rate=1e-8, random_seed=1,samples=ancient_samples)
	haplotypes=tree_sequence.haplotypes()
	haplotype_list=[]
	for i in haplotypes:
		haplotype_list.append(i)
	haplotype_list.append(tree_sequence.get_num_mutations())
	sim=haplotype_list
	print_list(sim)
	het_list=het_B_finder(sim)
	derived_list=derived_finder(sim,het_list)
	print("het list is",het_list)
	fab_k=fab_calc(sim,derived_list)
	print(fab_k)

#test_abc()





#fab functions go here

#sim=simulate()

def print_list(sim):
	for st in sim:
		print (st)
		print (" ")

#print_list()

def het_B_finder(sim):
	het_B_list=[]
	for base in range(0,int(sim[8])):
		if sim[2][base]!=sim[3][base]:
			het_B_list.append(base)
	return het_B_list
	
#het_list=het_B_finder()
#print("het list is", het_list)

#allele in B that is not in outgroup is the derived allele

def derived_finder(sim,het_list):
	ancestral_allele_list=[]
	ancestral_allele_position_list=[]
	derived_allele_list=[]
	for base in het_list:
		#check if outgroup is homozygous first
		if sim[4][base]==sim[5][base]:
			ancestral_allele_list.append(sim[4][base])
			ancestral_allele_position_list.append(base)
	for base in range(0,len(ancestral_allele_list)):
		if ancestral_allele_list[base]=="0":
			derived_allele_list.append("1")
		elif ancestral_allele_list[base]=="1":
			derived_allele_list.append("0")
		else:
			print ("ERROR IN DERIVED_FINDER FUNCTION")
	derived_allele_list=[derived_allele_list,ancestral_allele_position_list]
	return derived_allele_list

#derived_list=derived_finder()

#print ("derived list is",derived_list)

"""
def fab_calc():
	count=0


	for base in range(0,len(derived_list[1])):
		if sim[0][derived_list[1][base]]==derived_list[0][base]:
			count=count+1
		if sim[1][derived_list[1][base]]==derived_list[0][base]:
			count=count+1
	fab=count/len(derived_list[1])*0.5


	##testing..
allele_1=[]
allele_2=[]
	for base in range(0,len(derived_list[1])):
		allele_1.append(sim[0][derived_list[1][base]])
		allele_2.append(sim[1][derived_list[1][base]])
	allele=[allele_1,allele_2]
	for st in allele:
		print (st)
		print (" ")
	print(derived_list[0])

	for base in range(0,len(derived_list[1])):
		if allele_1[0][base]==derived_list[0][base]:
			count=count+1
		if allele_2[0][base]==derived_list[0][base]:
			count=count+1
	fab=count/len(derived_list[1])*0.5
	return fab
"""

def fab_calc(sim,derived_list):
	if len(derived_list[1])==0:
		return None
	count=0
	allele_1=[]
	allele_2=[]
	for base in range(0,len(derived_list[1])):
		allele_1.append(sim[0][derived_list[1][base]])
		allele_2.append(sim[1][derived_list[1][base]])
	allele=[allele_1,allele_2]

#	for st in allele:
#		print (st)
#		print (" ")
#	print(derived_list[0])
#	print(len(derived_list[1]))
#	print(len(allele[0]))
#	print(len(allele[1]))


	for base in range(0,len(derived_list[1])):
		if allele[0][base]==derived_list[0][base]:
			count=count+1
	#		print("count is ", count)
		if allele[1][base]==derived_list[0][base]:
			count=count+1
	#		print("count is ", count)
#	print("Final count is ", count)
	count=count*1.0
#	print("Did it work??", count)
	fab=count/len(derived_list[1])*0.5
	return fab

#print("fab stat is", fab_calc())



#simulate function returns a string list of haplotypes
"""
def simulate():
	tree_sequence=msprime.simulate(Ne=N_DNG,length=1e8,recombination_rate=1e-9,population_configurations=population_configurations,demographic_events=demographic_events,
            mutation_rate=1e-8, random_seed=50)
	haplotypes=tree_sequence.haplotypes()
	haplotype_list=[]
	for i in haplotypes:
		haplotype_list.append(i)
	haplotype_list.append(tree_sequence.get_num_mutations())
	return haplotype_list
"""

def d_calc(fab_k,fab_obs):
	if fab_k==None:
		return 9000
	else:
		return (fab_obs-fab_k)**2


