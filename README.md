# Introduction


CICADA (circRNA coding capability and product detection algorithm) software is a powerful tool designed to accurately distinguish the coding potential and products of circular RNAs (circRNAs). Recent studies have indicated the potential of circRNAs for translation. Certain circRNA-encoded proteins have been detected and found to have significant roles in the regulation of cancer cell behavior. However, the identification of circRNA-encoded proteins remains a significant challenge. Here, we present a pipeline aimed at assisting researchers in the exploration of unexplored translatable circRNAs.

### 1 CICADA install

***Note that we also provide a website for review to test our program in case the installation of CICADA is not smooth***
***[CICADA]http://121.196.54.69:8002/CICADA/home***

#### 1.1 Clone CICADA to local from github
```
git clone https://github.com/SunLab-biotool/CICADA.git
```

#### 1.2 Dependencies
```
conda > v3.8.5
Python3
R >= v3.5.0 
perl > v5.8
PhyloCSF
```

#### 1.3 Conda can be downloaded as part of the Anaconda or the Miniconda plattforms. We recommend to install miniconda3. Using Linux you can get it with:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```


This script is written for Python3.

#### 1.4 R and perl can be installed by conda :
```
conda install -c bioconda perl  
conda install r-base=3.5  
conda install r-randomForest  
conda install r-LncFinder  
```


#### 1.5 Other python packages can be installed by pip :  
```
pip install sklearn pandas rpy2
```

#### 1.6 Refer to the link below for installation of PhyloCSF:  

[PhyloCSF]https://github.com/cpockrandt/PhyloCSFpp


### 2 Usage  

```
python3 CICADA.py [-h] [-f] [-o] [-t]

options:
-h or --help:            show the help and exit
-f or --file :             input files
(Required) input files transcript nucleotide sequences
-o or --out             output files path
(Required) the path where you want to output
-t or --top              top rank protein products 
(Optional) the top rank protein products you want  

```

### 3 Example
You can use CICADA subroutines like our example:  
```
cd CICADA
python3 CICADA.py -f example/example.fa -t all -o example/example_out
```
