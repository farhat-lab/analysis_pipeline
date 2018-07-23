"""

Takes in threshold coverage and %of strain before you say a region is low coverage
First: List of strains of interest
Second: Coverage Threshold 
Third: Sequence threshold (At least this threshold % of the # of strains must have the coverage threshold or else we reject)
Fourth: Your reference genome name

Note this script assumes your path for the depth data will be at /n/data1/hms/dbmi/farhat/rollingDB/genomic_data/<strain_name>/depth/<strain_name>.depth.gz

"""
import pandas as pd
import sys
import csv
import os


arguments = sys.argv
list_of_strains  = arguments[1]
coverage_threshold = float(arguments[2])
threshold = arguments[3]
reference_name = arguments[4]

sequence = {}

def update(position, coverage):
	if(coverage < coverage_threshold):
		if(sequence[position] != 'F'):
			failure = sequence[position] - 1
			if(failure <= threshold):
				sequence[position] = 'F'
			else:
				sequence[position] = failure

def output():
	file = open('low_coverage.bed','w')

	start = None
	end = None
	end_sequence = False

	for key in sequence:
		if(sequence[key] == 'F'):
			if(not start):
				start = key
				end = key
			else:
				end = key
			end_sequence = True
		else:
			if(end_sequence):
				file.write('{}\t{}\t{}\n'.format(reference_name, start, end))
				start = None
				end = None
				end_sequence = False

#iterate through all the depth files in 
#directory = os.fsencode(location)
strains = [i.strip() for i in open(list_of_strains,'r').readlines()]

#path, dirs, files = next(os.walk(location))
file_count = len(strains)

threshold = file_count*(1-float(threshold))
import os.path

for file in strains:
    filename = "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/{}/depth/{}.depth.gz".format(file,file)
    if os.path.isfile(filename): 
    	df = pd.read_csv(filename, compression='gzip', header=0, sep='\t', quotechar='"')
    	columns = list(df.columns)
    	first =  columns[1]
    	second = columns[2]
    	current_sequence = pd.Series(df[second].values,index=df[first]).to_dict()
    	current_sequence[first] = second

    	for key in current_sequence:
    		if(key not in sequence):
    			sequence[key] = file_count
    		update(key, float(current_sequence[key]))
output()
