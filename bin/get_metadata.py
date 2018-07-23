import sys
import pandas as pd 

strains_of_interest = [i.rstrip().strip() for i in open(sys.argv[2],'r').readlines()]
drug = sys.argv[3]
df = pd.read_csv(sys.argv[1], sep='\t')
df.index = df['strain']

#Get one metadata that tells you whether your strain is resistant or not
metadata = df.loc[strains_of_interest][drug.upper()+'_R']
metadata.columns = ['strains',drug.upper()+'_R']
metadata.to_csv("drug_metadata.txt", sep="\t")

#Get another metadata that tells what date your strain was sampled
#metadata = df.loc[strains_of_interest]['date']
#metadata.columns = ['strains', 'date']
#metadata.loc['MT_H37Rv'] = 1998
#metadata.to_csv("date_metadata.txt", sep="\t")


