import sys

badfile = sys.argv[1]
referencefile = sys.argv[2]

numbers_to_strand  = {}

for line in open(referencefile,'r'):
	data = line.rstrip().split('\t')
	key = data[0]+"\t"+data[1]
	strand = data[2]
	numbers_to_strand[key] = strand

file = open(badfile+'_fixed','w')

for line in open(badfile,'r'):
	data = line.rstrip().split('\t')
	key = data[1]+"\t"+data[2]
	if(key in numbers_to_strand):
		strand = numbers_to_strand[key]
		print("WOO BOY WE HIT")
	else:
		if(int(data[1]) > int(data[2])):
			strand = "-"
		else:
			strand = "+"
	if(data[1] > data[2]):
		key = data[2]+"\t"+data[1]
	file.write(data[0]+'\t'+key+'\t'+strand+'\n')


