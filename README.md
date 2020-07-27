# PTMC
PTMC: Sequence based prediction of transmembrane protein crystallization propensity



# Installation
* Install Python 3.5 in Linux and Windows.
* Because the program is written in Python 3.5, python 3.5 with the pip tool must be installed first. 
* BBPpred uses the following dependencies: numpy, pandas, collection, biopython and sklearn You can install these packages first, by the following commands:
```
pip install numpy
pip install biopython==1.73 
pip install sklearn
pip install xgboost==
```

* SCRATCH (version SCRATCH-1D release 1.2) (http://scratch.proteomics.ics.uci.edu, Downloads: http://download.igb.uci.edu/#sspro)
```
wget http://download.igb.uci.edu/SCRATCH-1D_1.2.tar.gz
tar -xvzf SCRATCH-1D_1.2.tar.gz
cd SCRATCH-1D_1.2
perl install.pl
cd ..
```
* Run SCRATCH-1D on the provided test dataset:
```
    ../bin/run_SCRATCH-1D_predictors.sh test.fasta test.out 4
```
(if your computer has less than 4 cores, replace 4 by 1 in the command line above)

* The 4 output files:
```
- test.out.ss
- test.out.ss8
- test.out.acc
- test.out.acc20
```

# Running PTMC
open cmd in Windows or terminal in Linux, then cd to the PTMC-master/code folder which contains predict.py 

To predict general fasta sequences using our model, run: 

`python PTMC.py [custom predicting data in test.fasta format] [custom predicting data in test.acc20 format] [ predicting results in csv format]`


