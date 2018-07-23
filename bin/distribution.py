import os
import argparse
import itertools
import multiprocessing
import numpy as np
#from functools import partial

parser = argparse.ArgumentParser(description='Calculate all diversity score via average pairwise distance')
parser.add_argument('location',type=str, help='location of files to calculate score')

args = parser.parse_args()
#filee = open('parallel_cluster_distribution','w')
#numbers = []

manager = multiprocessing.Manager()
numbers = manager.list()

def score(pair):
	global numbers
        number_different = 0
	#print("GOING IN")
	a = np.array(list(pair[0]))
	b = np.array(list(pair[1]))
	#print(len(pair[0]))
	#print(len(pair[1]))
	number_different = np.sum(a != b)
	#print("COMPARED")
	numbers.append(number_different)
	#for i,j in zip(pair[0], pair[1]):
        #        if(i != j):
        #                number_different += 1
        #numbers.append(number_different)
	#print(number_different)
	return number_different

def process(filename):
	global numbers
        file = open(filename,'r').readlines()
        file = [i.strip().rstrip() for i in file if '>' not in i]
        if(len(file)):
                #numbers = []
		#manager = multiprocessing.Manager()
		#numbers = manager.list()
		#print("ABOUT TO DO:{}\t{}\t{}".format(filename, len(file), len(file[0])))
		pairs = itertools.combinations(file,2)
		pool = multiprocessing.Pool()
		pool.map(score, pairs )
		pool.close()
		pool.join()
		#print(numbers)
		#print(len(numbers))
		result = "{}\t{}\t{}\t{}".format(filename, sum(numbers)/len(numbers), len(file), len(file[0]))
		#filee.write(result+"\n")
		print(result)
		return sum(numbers)/len(numbers)
process(args.location)
#l_dirs = [process(args.location+"/"+i+"/concnenate.fasta") for i in os.listdir(args.location) if 'lineage' in i]
#pool = multiprocessing.Pool()
#pool.map(process, l_dirs)
#pool.close()
#pool.join()


#l_dirs = [process("/home/ye12/analysis_pipeline/"+i+"/concnenate.fasta") for i in os.listdir("/home/ye12/analysis_pipeline/") if 'lineage' in i]
#pool = multiprocessing.Pool()
#pool.map(process, l_dirs)
#pool.close()
#pool.join()

#process(args.location)
