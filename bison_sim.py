import msprime
import numpy as np
import random
import sys

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

Ne=10000
gen_time=8
T=20
mu=2e-8
rho=1e-8
sim_length=1e7
demes=5
num_samples=20


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



samples = list()
for i in range(demes):
	for j in range(num_samples):
		samples.append(msprime.Sample(population=i,time=0))

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
	print("inbred_het","outbred_het", "num_polymorph")
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
	random_chrom1=random.randint(1,demes*num_samples-1)
	random_chrom2=random.randint(1,demes*num_samples-1)
	inbred_het=het_calc(sim[0],sim[1],length=polymorph)/sim_length
	outbred_het=het_calc(sim[random_chrom1],sim[random_chrom2],length=polymorph)/sim_length
	print(inbred_het,outbred_het,polymorph)
#	print(" ")




bison_sim()


########

