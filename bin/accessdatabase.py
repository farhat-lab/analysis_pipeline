"""
Written by Yasha Ektefaie
June 2018

python3 accessdatabase.py <geography.py output> <country> <lineage> <drug>

Takes in the outputs of geography.py (the main strain info one), country, lineage, and drug of interest and then outputs the strains
associated with your query or says that exists--note that depending on whether you insert specific arguments or not will give you different responses
so you can put no lineage or drug of interest and get all countries...etc.

Note for drugs it checks if both suscetible and resistance for that drug exists or else it won't output anything

"""

import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep='\t')
country = sys.argv[2]

if(len(sys.argv) >  3):
	lineage = sys.argv[3]
else:
	lineage = None

if(len(sys.argv) >  4):
	drug = sys.argv[4].upper()
else:
	drug= None

result = None
error = "Uh oh, no results were found"

if(lineage and drug and country):
	result = df.loc[(df['country'] == country) & (df[lineage] == 1) & ((df[drug+"_R"] == 1) | (df[drug+"_S"] == 1))]['strain']
elif(country and lineage):
	result = df.loc[(df['country'] == country) & (df[lineage] == 1)]['strain']
elif(country and drug):
	result = df.loc[(df['country'] == country) & ((df[drug+"_R"] == 1) | (df[drug+"_S"] == 1))]['strain']
elif(lineage and drug):
	result = df.loc[(df[lineage] == 1) & ((df[drug+"_R"] == 1) | (df[drug+"_S"] == 1))]['strain']
elif(lineage):
	result = df.loc[(df[lineage] == 1)]['strain']
elif(country):
	result = df.loc[(df['country'] == country)]['strain']
elif(drug):
	result = df.loc[(df[drug+"_R"] == 1) | (df[drug+"_S"] == 1)]['strain']

if(drug):
	suscetible = False
	resistant = False
	suscetible_number = 0
	resistant_number = 0

	for strain in list(result):
		answer = df.loc[df['strain'] == strain]
		if(answer[drug+"_R"].item() == 1):
			resistant = True
			resistant_number += 1
		elif(answer[drug+"_S"].item() == 1):
			suscetible = True
			suscetible_number += 1
	total = suscetible_number+resistant_number
	if(not(suscetible and resistant)):
		result = None
		error = "Sorry the drug you were interested in did not have both resistant and suscetible strains in our database"
	elif(suscetible_number/float(total) < 0.2):
		result = None
		error = "Sorry the drug you were interested in had less than 20% suscetible strains"
if(result is not None and not result.empty):
	#print(list(result))
	for i in list(result):
		print(i)
else:
	print(error)
