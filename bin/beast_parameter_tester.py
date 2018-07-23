import pandas as pd
import argparse
import sys
import re

parser = argparse.ArgumentParser()
parser.add_argument("input_file",type=str, help="Path to the input.xml for the BEAST File")
parser.add_argument("first_test", type=str, help="Path to the first test file")
args = parser.parse_args()

with open(args.first_test, 'r') as f:
	first_ESS_value = [line.strip().rstrip() for line in f]
	#first_ESS_value = first_ESS_value[1:]
	first_ESS_Value = [i for i in first_ESS_value if i]
	first_ESS_value = [float(i) for i in first_ESS_value if not re.search("[^a-zA-Z0-9.]+",i)]

minimum_first = min(first_ESS_value)
while(minimum_first == 1):
	first_ESS_value.remove(1)
	minimum_first = min(first_ESS_value)

if(minimum_first < 200):
	if(minimum_first < 40):
		print(40)
	elif(minimum_first > 100):
		print(100)
	else:
		print(int(minimum_first))
else:
	print("True")
