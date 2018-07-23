import sys

fasta = {}
#https://stackoverflow.com/questions/29333077/reading-a-fasta-file-format-into-python-dictionary

with open(sys.argv[1]) as file_one:
    for line in file_one:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_sequence_name = line[1:]
            if active_sequence_name not in fasta:
                fasta[active_sequence_name] = []
            continue
        sequence = line
        fasta[active_sequence_name].append(sequence)

tree = open(sys.argv[2]).readlines()
tree = tree[0]


from bs4 import BeautifulSoup

file = open(sys.argv[3])
e = BeautifulSoup(file, "xml")

length_of_alignment = []

for strain in fasta:
    new_tag = e.new_tag("taxon")
    new_tag['id'] = strain
    e.data.taxa.append(new_tag)

    new_tag = e.new_tag("sequence")
    sub_tag = e.new_tag("taxon")
    sub_tag['idref']=strain

    new_tag.append(sub_tag)
    new_tag.append(fasta[strain][0])
    length_of_alignment.append(len(fasta[strain][0]))

    e.data.alignment.append(new_tag)


e.rescaledTree.newick.append(tree)

# print(length_of_alignment)
length = length_of_alignment[0]
A = 33*(4411532 - length)/2
T = .33*(4411532 - length)/2
G = .66*(4411532 - length)/2
C = .66*(4411532 - length)/2


value = str(A)+" "+str(C)+" "+str(G)+" "+str(T)
e.mergePatterns.constantPatterns.counts.parameter['value']=value
# print("boom done")

file = open('input.xml','w')
file.write(e.prettify())





