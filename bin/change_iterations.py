import argparse
from bs4 import BeautifulSoup

parser = argparse.ArgumentParser()
parser.add_argument("input_file",type=str, help="Path to the input.xml for the BEAST File")
parser.add_argument("number_iterations", type=int, help="The number of iterations you want mcmc to run")
args = parser.parse_args()

e = BeautifulSoup(open(args.input_file),"xml")
e.mcmc['chainLength']=str(args.number_iterations)

file = open('input.xml','w')
file.write(e.prettify())



