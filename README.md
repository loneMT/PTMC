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
pip install xgboost
```

DISOPRED (version 3.16) (http://bioinfadmin.cs.ucl.ac.uk/downloads/DISOPRED/)
```
Run `wget http://bioinfadmin.cs.ucl.ac.uk/downloads/DISOPRED/DISOPRED3.16.tar.gz`
Run `tar -xvzf DISOPRED3.16.tar.gz`
Run `cd DISOPRED/src/`
Run `make clean; make; make install`
```
    
