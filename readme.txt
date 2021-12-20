#Amna Farooq 2020
avg_read_coverage.py 


This script is used for
1) Cleaning and sorting input methylation data (removes unidentified chromosomes, removes unnecessary data except the required columns)
2) Extracting data of single chromosome from input methylation files.
3) Split input methylation file chromosome wise.
4) Calculates average read count between the range 1-5, 6-10, 11-15, 16-20 for quality check of input data.

NOTE: Input file must have following order of columns : chr, chr_start, chr_end, methylation percentage, total reads.
-All tasks can be performed alone and collectively. However, task 2 and 3 can not be performed collectively.


For requirements see : requirements.txt


usage: avg_read_coverage.py [-h] --In_path IN_PATH --Out_path OUT_PATH --Org
                            ORG [--Chr CHR] [--E E] [--Avg_read AVG_READ]

This script is used for 1) Cleaning and sorting input methylation data
(removes unidentified chromosmes) 2) Extracting Chr wise data from input
methylation files. 3) Calculates average read count between the range 1-5,
6-10, 11-15, 16-20. NOTE: Input file must have following order of columns :
chr, chr_start, chr_end, meth_perc, t_reads


required arguments:
  --In_path IN_PATH    Enter the path for input file (default: None)
  --Out_path OUT_PATH  Enter the path for output files (default: out)
  --Org ORG            Enter hum for human, mou for mouse. Default is hum.
                       (default: hum)

optional arguments:
  --Chr CHR            Enter the chromosome number to extract data for
                       specific chromosome or 'all' to split the input
                       chromosome wise or 'none' to treat input file as it is.
                       Default = none (default: none)
  --E E                Extract sorted version of primary input file [y/n].
                       Default = n (default: n)
  --Avg_read AVG_READ  If you want to know average read coverage in the input
                       file, type 'y'. (default: None)


  -h, --help           show this help message and exit





			*****Example******

****For cleaning and sorting input :
Input: data-sample.txt

Output: out/chr/chr_data-sample.txt

Command:
 	python avg_read_coverage.py --In_path data-sample.txt --Out_path out  --Org hum



****For cleaning sorting input and extracting data for one chromosome :

Input: data-sample.txt

Output: out/chr17/chr17_data-sample.txt

Command:

python avg_read_coverage.py --In_path data-sample.txt --Out_path out  --Org hum --Chr 17




****For cleaning sorting input and splitting data chromosome wise :

Input: data-sample.txt

Output: out/   (respective chromosome folders)

Command:
python avg_read_coverage.py --In_path data-sample.txt --Out_path out  --Org hum --Chr all

*Output for this option is not included in demo as it would produce many files.



****For cleaning sorting input and calculating average read count :

Input: data-sample.txt

Output: out/chr/chr_data-sample.txt
	ARC_data-sample.txt.png

Command:
	python avg_read_coverage.py --In_path data-sample.txt --Out_path out  --Org hum --Avg_read y




