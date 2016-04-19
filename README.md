# Mr Mosaic
#### A Genomic Mosaic Structural Variant Caller for Massively Parallel Sequencing Data

## Authors

Dan King (Creator, Developer)

Alejandro Sifrim (Developer)

Tomas Fitzgerald (Developer)

Matthew Hurles (Group Leader)

We are affiliated with the [Wellcome Trust Sanger Institute](http://www.sanger.ac.uk/science/groups/hurles-group), Cambridge, United Kingdom

## Abstract

Structural mosaic abnormalities are large post-zygotic mutations present in a subset of cells and have been implicated in developmental disorders. While such mutations are routinely assessed in clinical diagnostics using cytogenetic or microarray testing, an adequate method using targeted or whole-genome sequencing data is lacking. Here, we present a method to detect structural mosaic abnormalities using deviations in allele fraction and read coverage from next generation sequencing data. 

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

MrMosaic works in 2 steps:


1. Extract depth and B-allele fraction from BAM files
2. Call mosaic SV events by measuring deviations in copy number and allele frequency and running
a segmentation algorithm to determine consecutive stretches of datapoints deviating from the null.

####Step 1: Extracting information from BAM/CRAM files

The following command takes either an indexed bam or cram file and a file containing the
 captured regions as inputs and outputs the necessary metrics (e.g. in this example piped into a file with the .mrm extension).

```
cd mrmosaic
./bin/cnv_baf_data -R ./resources/cnv_baf_bait_region_file.txt sample.bam > sample.mrm
```

####Step 2: Calling mosaic structural variants
In this step MrMosaic uses the deviations in coverage (by using a group of other samples as a reference) and deviations in B-allele proportion of from the expected heterozygous 0.5
to call structural variants.

This script takes the output of the previous step into account and a configuration file which needs to be edited to point to the correct directories:
* `resource_dir`: directory containing the ecdf_GADAamp.merged.gada.10.R3.Robj and ecdf_nprobes.merged.gada.10.R3.Robj files
* `reference_dir`: directory containing a set of cnv_baf_data output files which MrMosaic can use as a reference for the mean coverage for the bait. Ideally these are 10+ samples sequenced
using the same capturing platform/mean sequencing depth. 
* `output_dir`: directory to output the resulting SV calls to.\

```
cd mrmosaic/mrmosaic/scripts/
./run_mrmosaic.sh --input sample.mrm --config /nfs/users/nfs_a/as33/Projects/MrMosaic/config.yaml
```


##Output
Two files are outputted after running MrMosaic:

1. A file containing all the unfiltered candidate SV calls
2. A file containing candidate SV calls based on their loglikelihood (<= -8) of being a false positive given the number of baits the event spans and the amplitude of the total signal.
This is based on the observation that in our simulation study we observed that most false positives consisted of events spanning few bait regions and of low signal intensity.

Both files are formatted identically and consist of the following columns:

* `Chrom`: chromosome where the event is located
* `IniPos`: putative start position of event
* `EndPos`: putative end position of event
* `LRR`: mean log2ratio for the event
* `bdev`: mean b-allele frequence deviation for the event
* `clon`: putative clonality of the event given the putative event type (e.g. gain, loss, LOH)
* `type`: putative event type (e.g. gain, loss, LOH)
* `nProbe`: number of baits spanning the event
* `GADAamp`: mean amplitude of the segment called by GADA
* `p_nprobes`: probability of event being a false positive given the number of baits spanned by the event (based on simulations). 
* `p_GADAamp`: probability of event being a false positive given the amplitude of the signal (based on simulations).
* `ll`: log-likelihood of event being a false positive given the number of probes and the amplitude. 





