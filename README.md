# Mr Mosaic
#### A Genomic Mosaic Structural Variant Caller

## Authors

Dan King (Creator, Developer)

Alejandro Sifrim (Developer)

Tomas Fitzgerald (Developer)

Matthew Hurles (Group Leader)

We are affiliated with the [Wellcome Trust Sanger Institute](http://www.sanger.ac.uk/science/groups/hurles-group), Cambridge, United Kingdom

## Abstract

## What does it do?

## Installation
MrMosaic requires a C-compiler such as gcc, samtools-1.3 and a couple of R packages.

MrMosaic requires the following R packages:

* gada
* yaml
* argparse
* devtools

To install gada run the following command in the R console:
```R
install.packages("gada", repos="http://R-Forge.R-project.org")
```
GADA requires a FORTRAN compiler to build correctly.

Some newer OSX versions complain about not finding the gfortran-4.8 compiler. You can compile this by using [homebrew](http://brew.sh) to install gcc:
```
brew install gcc
```
And then pointing R to the correct directories in `~/.R/MakeVars`:
```
F77=gfortran
FC=gfortran
FLIBS="-L/Users/as33/homebrew/lib/gcc/5 -lgfortran -lquadmath -lm"
```

Once the requirements are installed, MrMosaic can be installed as follows:
```
git clone https://github.com/asifrim/mrmosaic
cd mrmosaic
make
```
Alternatively (if you don't want to install git, shame...) one can download the tarball and install MrMosaic from there:
```
wget https://github.com/asifrim/mrmosaic/archive/master.zip
unzip master.zip
cd mrmosaic-master
make
```




 


## Usage
