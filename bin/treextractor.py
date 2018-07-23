from ete2 import Tree
from nexus import NexusReader
import sys
import re
import csv

arguments = sys.argv
n = NexusReader(arguments[1])

#First get that resistance data--figure out which strain is resistant to what
drug_resistance = csv.DictReader(open(arguments[2]),delimiter='\t', fieldnames = ( "strain", "drug"))

#Then We gather a mapping from name to number
number = 1
number_name = {}
for i in n.taxa.taxa:
	number_name[i] = number
	number += 1
keys = list(number_name.keys())
keys.sort(reverse=True)

strain_to_resistance = {}

for row in drug_resistance:
	strain_to_resistance[row['strain']] = row["drug"]

#Next get the date strains were taken data to figure out the oldest strain that the tree is based upon
#Maha said not to do this so let's default to 2018
# reference_date = csv.DictReader(open(arguments[3]),delimiter='\t', fieldnames = ( "strain", "date"))

# youngest_date = 0
# for row in reference_date:
# 	if(row['date'] != 'missing' and float(row['date']) > youngest_date and row['strain'] in number_name):
# 		youngest_date = float(row['date'])

# if(not youngest_date):
# 	print("UH OH THERE IS NO DATES ON FILE FOR YOUR STRAINS!")
reference_date = 2018

#We parse the tree to pull the info that we actually want
basic = n.trees[0]
basic = basic[13:]
finaltree = "".join(re.compile("\[.*?]").split(basic))

# basic = basic.split(',')
# finaltree = ''
# for i in l
# finaltree = [for i in l]

# for i in basic:
# 	if('[&height=' in i):
# 		print(i)
# 		i = i.replace('[&height=',':')
# 		if(')' not in i):
# 			i = i +','
# 		finaltree += i
# 	elif('height=' in i):
# 		print(i)
# 		i = i.replace('height=',':')
# 		if(')' not in i):
# 			i = i +','
# 		finaltree += i
# 	elif('(' in i):
# 		print(i)
# 		finaltree = finaltree+'('+','
# 	elif(')' in i):
# 		print(i)
# 		finaltree = finaltree+')'
# finaltree = finaltree.replace(',)',')')
# finaltree = finaltree[:len(finaltree) - 1]
# finaltree += ';'

print(finaltree)
t = Tree(finaltree)


oldest = 0
oldest_strain = None

print(strain_to_resistance)

for strain in strain_to_resistance:
	if(int(strain_to_resistance[strain]) and strain in number_name):
		distance = t.search_nodes(name=str(number_name[strain]))[0].dist
		if(distance > oldest):
			oldest = distance
			oldest_strain = strain

print("At long last the jounry is over")
print(oldest_strain)
print(oldest)
print(float(reference_date)-oldest)




# #Step 2--We gather a mapping from number to height
# heights = re.findall('[0-9]\[&height=[0-9]+.[0-9E-]+', basic)
# name_to_height = {}

# for height in heights:
# 	data = height.split("[&height=") 
# 	name_to_height[data[0]] = data[1]


# results = re.findall('([(]+[0-9]+)|([)]+)|([:][0-9]+.[0-9]+,[0-9])|([:][0-9]+.[0-9]+,)',basic)
# resultat = ''

# for i in results:
# 	for x in i:
# 		if(":" in x):
# 			find = re.findall(',[0-9]+', x)
# 			if(len(find) == 0):
# 				resultat += ","
# 			else:
# 				find = find[0]
# 				for key in keys:
# 					if(str(key) in find):
# 						height = name_to_height[str(key)]
# 						find = find.replace(str(key), number_name[key]+":"+height)
# 						break
# 				resultat += find
# 		else:
# 			for key in keys:
# 					if(str(key) in x):
# 						height = name_to_height[str(key)]
# 						x= x.replace(str(key), number_name[key]+":"+height)
# 						break
# 			resultat += x
# resultat += ';'

#print(resultat)
# t = Tree(resultat)
# print(resultat)

