#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 28-00:00                         # Runtime in D-HH:MM format
#SBATCH -p long                           # Partition to run in
#SBATCH --mem=80G                          # Memory total in MB (for all cores)
#SBATCH -o Bangladeshlineage1AMIKACIN.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e Bangladeshlineage1AMIKACIN.err                 # File to which STDERR will be written, including job ID

#Welcome to the analysis pipeline
#This pipeline was created by Yasha Ektefaie in June 2018
#yasha_ektefaie@berkeley.edu if you find any issues
#ENJOY
#Read the README if you wanna know what this all does


#The country variable is which country you want to look at the strains of
#The lineage variable is which lineage within that country you want to look at the strains of
#The drug variable is which drug within that lineage you want to look at the strain of

#The coverage_threshold variable is used for what coverage do you consider as too little
#The percent threshold is at least what percentage of strains should have the coverage below the threshold at a position before you identify the region as a low coverage region
#The name_of_lineage is where you put the name of the gene that you want on the very left column of the bed file
#The path to python is because some of my scripts use panda and the default installation on o2 doesn't have panda so give me a path of python that has pandas
#The path to bedtools merge is needed for the bed file processing part
#The path for the fastml is needed for the final fastml part of the pipeline
export LD_LIBRARY_PATH=$HOME/bin/beagle-lib/lib:$LD_LIBRARY_PATH

country="$1" 
lineage="$2"
drug="$3"
truedrug="$4"
coverage_threshold=10
percentage_threshold=0.05
name_of_lineage="NC_0002.1"
path_to_python="../../anaconda2/bin/python2.7"
path_to_bedtools_merge="../../bin/bedtools2/bin/mergeBed"
path_to_fastml="../../bin/FastML.v3.1/www/fastml/FastML_Wrapper.pl"

filename="$country$lineage$drug"
#rm -r "$filename"
mkdir "$filename"
cd "$filename"

#################################
#				#
#    Producing SNP Alignment	#
#				#
#################################

echo "Accessing Database"
$path_to_python ../bin/accessdatabase.py ../bin/strain_info.tsv "$country" "$lineage" "$truedrug" > strains.txt
echo "Creating Metadata"
$path_to_python ../bin/get_metadata.py ../bin/strain_info.tsv strains.txt "$truedrug"
echo "Accessed now creating depth file"
$path_to_python ../bin/yasha_depth_organizer.py strains.txt "$coverage_threshold" "$percentage_threshold" "$name_of_lineage"
echo "Now fixing the bed files and merging with the low coverage areas"
$path_to_python ../bin/fix_bed.py low_coverage.bed ../bin/h37rv_genome_summary.txt
cat ../bin/drr.bed >> low_coverage.bed_fixed
cat low_coverage.bed_fixed | sort -k1,1 -k2n > almost_final.bed 
cat almost_final.bed|awk '{print $1,$2,$3,$5,$6,$4}'|tr ' ' '\t' > ugh.bed
cat ugh.bed | sort -k1,1 -k2n > penultimate.bed
$path_to_bedtools_merge -i penultimate.bed -s > final.bed
rm almost_final.bed
rm ugh.bed
rm penultimate.bed
rm low_coverage.bed_fixed
echo "Using Maha Script to concatenate with exclusion the SNPS"
perl ../bin/modified_snpConcatenater_w_exclusion_frompilonvcf_2.0.pl final.bed ../bin/blank.txt SNP WHOLE > concnenate.fasta

#################################
#                               #
#    Applying Diversity Filter  #
#                               #
#################################


echo "NOW RUNNING DIVERSITY FILTER: Depending on whether the sequence is sufficiently diverse or not we will or won't run"

mean=100
std=50
threshold=`expr $mean - $std`
result=`$path_to_python bin/distribution.py concnenate.fasta`
if [ $result -lt $threshold ]
then
	echo "I'm sorry but this failed the genetic diversity filter exiting now"
	exit 0
fi
#echo "Running FastMl on the concanenate.fasta"
#mkdir fastml
#perl $path_to_fastml --MSA_File /home/ye12/analysis_pipeline/$filename/concnenate.fasta --seqType NUC --outDir /home/ye12/analysis_pipeline/$filename/fastml/
#rm final.bed
#rm concnenate.fasta
#rm strains.txt

#################################
#                               #
# Prepping BEAST Input File	#
#                               #
#################################

megacc="../../bin/megacc"
BEAST="../../bin/BEASTv1.10.0/lib/beast.jar"
LOG="../../bin/BEASTv1.10.0/bin/loganalyser"
BEAGLE="/home/ye12/bin/beagle-lib/lib/"
raxml="/home/ye12/bin/standard-RAxML-master/raxmlHPC"

#First we construct a simple NJ tree to use for the clock test
echo "Constructing simple NJ Tree"
mkdir MEGAresult

$raxml -s concnenate.fasta -m GTRGAMMA  -n raxml -p 35425 -o MT_H37Rv
wait 

echo "Moving to test clock rate"
#Then we test clock see if you want a relaxed or a not relaxed clock rate

#cd MEGAresult
treefile="$(ls RAxML_bestTree.raxml)"
#cd ..

echo "Testing Clock Rate"
$megacc -a ../bin/ML_clock_test_coding.mao -d concnenate.fasta -t "$treefile" -g ../bin/group.txt -o MEGAresult
wait 

#Next we find the file name that was outputted and check whether the null was accepted or rejected
echo "Finished Testing Moving on"
echo "Building Xml file for Beast"

#cd MEGAresult
filename="$(find ./MEGAresult/ -maxdepth 1 -name "*_clockTest*" -print)"
#cd ..

#If accepted we make the xml file of interest one way if not we make it another way
if [ grep -q "The null hypothesis of equal evolutionary rates throughout the tree was rejected at a 5% significance level" "$filename" ];
then
        echo "Relaxed Clock Rate"
        $path_to_python ../bin/fasta_to_beast.py concnenate.fasta "$treefile" ../bin/avika_uncorrelated_relaxed_rate.xml
else
        echo "Strict Clock Rate"
        $path_to_python ../bin/fasta_to_beast.py concnenate.fasta "$treefile" ../bin/avika_strict_clock_rate.xml
fi

wait 

#################################
#                               #
#        Running BEAST		#
#                               #
#################################

#export LD_LIBRARY_PATH=$HOME/bin/beagle-lib/lib:$LD_LIBRARY_PATH
echo "Running BEAST with parameter setup"
#java -Djava.library.path="$BEAGLE" -jar "$BEAST" input.xml
iteration=5000000
converge=False

while [ $converge == False ]
do
        $path_to_python ../bin/change_iterations.py input.xml $iteration
        java -Djava.library.path="$BEAGLE" -jar "$BEAST" -overwrite input.xml
        $LOG saveme.log > savemysoul
        sed -n -e '19,36 p' savemysoul | awk '{ print $7 }' > beasttest1
        cat beasttest1
	output=`$path_to_python ../bin/beast_parameter_tester.py input.xml beasttest1`
        if [ $output == True ]
        then
                echo "We found convergence"
                converge=$output
                echo "Just set converge to true heading out"
        else
                echo "Ey sorry did didn't find convergence"
                (( iteration=$iteration*(200 / $output) ))
                echo "Now we will do more iterations:"
                echo $iteration
                echo "This was the ESS value"
                echo $output
        fi
	cp saveme.log beastoutput.log
done
wait
 
#################################
#                               #
#    Interpreting BEAST Results #
#     Treeannotator/Bootstrap   #
################################# 

echo "Running Treeannotator"
burn=0
echo $iteration
(( burn=$iteration*20/100 ))
echo $burn
../../bin/BEASTv1.10.0/bin/treeannotator -burnin $burn saveme.trees final.tree
wait
echo "RUNNING BOOTSTRAP ANALYSIS"
raxml="/home/ye12/bin/standard-RAxML-master/raxmlHPC"
$raxml -f a -N 10000 -p 123 -x 123 -o 'MT_H37Rv' -s concnenate.fasta -n TEST -m GTRGAMMA
sumtrees.py -d0 -p -t final.tree --force-rooted  -o final_bootstrap.tree RAxML_bootstrap.TEST
echo "MOVING NOW"
cd ..
mv "$country$lineage$drug" output
