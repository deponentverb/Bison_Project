import msprime
import numpy as np
import random
import sys

#functions

def print_list(sim):
	for st in sim:
		print (st)

def het_calc(chrom1,chrom2):
	count=0
	length=len(chrom1)
	if len(chrom1)!= len(chrom2):
		print ("error")
	for pos in range(0,length-1):
		if chrom1[pos]!=chrom2[pos]:
			count=count+1
	return count

#obtain the average pairwise het between haplotypes a to b
def pair_wise_het (start,end,chrom_list):
	count=0
	for i in range (start,end):
		count=count+het_calc(chrom_list[start],chrom_list[start+i],)




#setting parameters
"""
Ne=100000
gen_time=8
mu=2e-8
rho=1e-8
sim_length=1.5e8
demes=5
num_samples=20
num_sim=100
bottleneck_size=20
bottleneck_length=50


"""
#testing parameters

Ne=10
gen_time=8
bottleneck_length=50
mu=1e-2
rho=0
sim_length=2
demes=5
num_samples=2
num_sim=10
bottleneck_size=20
#"""

def pop_config(bottleneck_size,Ne):
	population_configurations = [
		msprime.PopulationConfiguration(
			initial_size=bottleneck_size),
		msprime.PopulationConfiguration(
			initial_size=bottleneck_size),
		msprime.PopulationConfiguration(
			initial_size=bottleneck_size),
		msprime.PopulationConfiguration(
			initial_size=bottleneck_size),
		msprime.PopulationConfiguration(
			initial_size=Ne),
		]
	return population_configurations



samples = list()
for i in range(demes):
	for j in range(num_samples):
		samples.append(msprime.Sample(population=i,time=0))

def demo_events(bottleneck_length):
	demographic_events=[
		msprime.MassMigration(
			time=bottleneck_length,source=0,destination=4,proportion=1),
		msprime.MassMigration(
			time=bottleneck_length,source=1,destination=4,proportion=1),
		msprime.MassMigration(
			time=bottleneck_length,source=2,destination=4,proportion=1),
		msprime.MassMigration(
			time=bottleneck_length,source=3,destination=4,proportion=1),
	]
	return demographic_events

def bison_sim(bottleneck_size,bottleneck_length,Ne):
	print("class","het", "num_polymorph")
	population_configurations=pop_config(bottleneck_size,Ne)
	demographic_events=demo_events(bottleneck_length)
	for i in range (num_sim):
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
	#	print_list(sim)
	#	print(" ")
	#	print(sim[0][0])
	#	print(sim[1])
	#	print(polymorph)
	#	print(len(sim[1]))
		random_chrom1=random.randint(1,demes*(num_samples-1)-1)
		random_chrom2=random.randint(1,demes*(num_samples-1)-1)
		while(random_chrom2==random_chrom1):
			random_chrom2=random.randint(1,demes*num_samples-num_samples-1)
		inbred_het=het_calc(sim[0],sim[1])/sim_length
		outbred_het=het_calc(sim[random_chrom1],sim[random_chrom2])/sim_length
		print("inbred",inbred_het,polymorph)
		print("outbred",outbred_het,polymorph)
	#	print(" ")




bison_sim(bottleneck_size=bottleneck_size,bottleneck_length=bottleneck_length,Ne=Ne)


########

