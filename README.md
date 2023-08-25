# Introduction
CICADA (circRNA coding capability and product detection algorithm) software is a powerful tool designed to accurately distinguish the coding potential and products of circular RNAs (circRNAs). Recent studies have indicated the potential of circRNAs for translation. Certain circRNA-encoded proteins have been detected and found to have significant roles in the regulation of cancer cell behavior. However, the identification of circRNA-encoded proteins remains a significant challenge. Here, we present a pipeline aimed at assisting researchers in the exploration of unexplored translatable circRNAs.

### CICADA install
#### 2.1 Dependencies
```
conda > v3.8.5
Python3
R >= v3.5.0 
perl > v5.8
PhyloCSF
```

2.2 Conda can be downloaded as part of the Anaconda or the Miniconda plattforms. We recommend to install miniconda3. Using Linux you can get it with:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```


This script is written for Python3.

2.3 R and perl can be installed by conda :
```
conda install -c bioconda perl  
conda install r-base=3.5  
conda install r-randomForest  
conda install r-LncFinder  
```


2.4 Other python packages can be installed by pip :  
```
pip install sklearn pandas rpy2
```

2.5 Refer to the link below for installation of PhyloCSF:  

[PhyloCSF]https://github.com/cpockrandt/PhyloCSFpp


### USAGE  

```
python3 qyx_fold_test3_circ_9_aliyun.py [-h] [-f] [-o] [-t]

options:
-h or --help:            show the help and exit*************************
-f or --file :             input files
(Required) Input files transcript nucleotide sequences
-o or --out             output files path
(Required) The path where you want to output
-t or --top              top rank protein products 
(Optional) The top rank protein products you want  

```

### Example
You can use CICADA subroutines like our example:  
```
Python3 circ_9_aliyun.py -f example/example.fa -g -o example/test
```







