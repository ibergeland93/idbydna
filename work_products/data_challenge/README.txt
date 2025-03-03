### IDbyDNA DATA CHALLENGE ###

Step 1: Download FastA file from google drive --> https://drive.google.com/file/d/1zYt0XHzPOI37klCEpy0pjQIBpaf5wrz4/view?usp=sharing

Step 2: Create Python script that computes the following

A) Total number of 25mers
B) Number of distinct mers
C) K-mer with the highest count

Step 3: Check work with Jellyfish

############################ Step 1 ##################################

# Download FastA file 

(base) Indias-MacBook-Pro-780:~ ibergeland$ mkdir -p Desktop/idbydna
(base) Indias-MacBook-Pro-780:~ ibergeland$ cd Desktop/idbydna/

# Create README file
(base) Indias-MacBook-Pro-780:idbydna ibergeland$ nano README.txt


########################### Step 2 ######################################

# Create python script using Jupyter Notebook

#update data rate limit 
jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000000


## See idbydna_data_challenge.ipynb



########################## Step 3 #######################################

#Verify work using Jellyfish


##Install Jellyfish

(base) Indias-MacBook-Pro-780:~ ibergeland$ pip install jellyfish

#Collecting jellyfish
#  Downloading https://files.pythonhosted.org/packages/3f/80/bcacc7affb47be7279d7d35225e1a932416ed051b315a7f9df20acf04cbe/jellyfish-0.7.2.tar.gz (133kB)
#    100% |████████████████████████████████| 143kB 7.5MB/s 
#Building wheels for collected packages: jellyfish
#  Running setup.py bdist_wheel for jellyfish ... done
#  Stored in directory: /Users/ibergeland/Library/Caches/pip/wheels/e8/fe/99/d8fa8f2ef7b82a625b0b77a84d319b0b50693659823c4effb4
#Successfully built jellyfish
#Installing collected packages: jellyfish
#Successfully installed jellyfish-0.7.2


(base) Indias-MacBook-Pro-780:Jellyfish-2.3.0 ibergeland$ ./configure --prefix=$HOME
(base) Indias-MacBook-Pro-780:Jellyfish-2.3.0 ibergeland$ make
(base) Indias-MacBook-Pro-780:Jellyfish-2.3.0 ibergeland$ sudo make install


# Configure & Verify Jellyfish

(base) Indias-MacBook-Pro-780:idbydna ibergeland$ . ~/.bash_profile
(base) Indias-MacBook-Pro-780:idbydna ibergeland$ which jellyfish
#/Users/ibergeland/bin/jellyfish

(base) Indias-MacBook-Pro-780:jellyfish-2.3.0 ibergeland$ nano ~/.bash_profile
(base) Indias-MacBook-Pro-780:jellyfish-2.3.0 ibergeland$ . ~/.bash_profile
(base) Indias-MacBook-Pro-780:jellyfish-2.3.0 ibergeland$ echo $PATH
#/opt/local/bin:/opt/local/sbin:/Users/ibergeland/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.7/bin:/Users/ibergeland/anaconda3/bin:/Users/ibergeland/anaconda3/condabin:/opt/local/bin:/opt/local/sbin:/Users/ibergeland/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.7/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Users/ibergeland/Library/Python/2.7/bin:/opt/X11/bin:/Users/ibergeland/bin

(base) Indias-MacBook-Pro-780:jellyfish-2.3.0 ibergeland$ which jellyfish
#/Users/ibergeland/bin/jellyfish

(base) Indias-MacBook-Pro-780:jellyfish-2.3.0 ibergeland$ jellyfish -h
#Usage: jellyfish <cmd> [options] arg...
#Where <cmd> is one of: count, bc, info, stats, histo, dump, merge, query, cite, mem, jf.
#Options:
#  --version        Display version
#  --help           Display this message


# Create .jf file for further analysis 

(base) Indias-MacBook-Pro-780:idbydna ibergeland$ jellyfish count -m 25 -s 100M -t 10 -C SRR1748776.fa
(base) Indias-MacBook-Pro-780:idbydna ibergeland$ ls

#README.txt			SRR1748776.fa			idbydna_data_challenge.ipynb	mer_counts.jf


# output all the counts for all the k-mers in the file
(base) Indias-MacBook-Pro-780:idbydna ibergeland$ jellyfish dump mer_counts.jf > mer_counts_dumps.fa

# query the counts of a particular 25mer
## There are inconsistensies between some of the outputs. The first two examples do not align with python output but the second two answers do.

## Inconsistent
(base) Indias-MacBook-Pro-780:idbydna ibergeland$ jellyfish query mer_counts.jf CGGAAGAGCGGTTCAGCAGGAATGC
#CGGAAGAGCGGTTCAGCAGGAATGC 72003 --> jellyfish output , python output ---> 62777

(base) Indias-MacBook-Pro-780:idbydna ibergeland$ jellyfish query mer_counts.jf GGAAGAGCGGTTCAGCAGGAATGCC
#GGAAGAGCGGTTCAGCAGGAATGCC 71556 --> jellyfish output , python output -> 62513

## Consistent
(base) Indias-MacBook-Pro-780:idbydna ibergeland$ jellyfish query mer_counts.jf GTTCAGCAGGAATGCCGAGACCGGA
GTTCAGCAGGAATGCCGAGACCGGA 53116 --> jellyfish output , python output -> 53116

(base) Indias-MacBook-Pro-780:idbydna ibergeland$ jellyfish query mer_counts.jf CGAGACCGGATAGCGATCTCGTATG
CATACGAGATCGCTATCCGGTCTCG 48062 ---> jellyfish output , python output --> 48062
