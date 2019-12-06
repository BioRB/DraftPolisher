# DraftPolisher

**DraftPolisher** is a tool for fast polishing of draft sequences.

## Description

DraftPolisher can produce an improved consensus sequence for a draft  linear or circular genome assembly. 

## Prerequisites

DraftPolisher is designed for the polishing of draft sequences. Two different releases of the tool are present: DraftPolisher.py and DraftPolisher_cov.py. The cov (i.e. coverage version)  takes into account also the coverage of the sequences from the reference database used for the polishing.The "cov" release was designed based on the SPAdes assembling output file format, where the coverage value is reported in the last part of the sequences IDs, so if using this version of the tool take care of the reference sequences file format (see SPAdes output format for more details). A general prerequisite for all the versions of the tool is the formatting of query and subject fasta files. It is mandatory to put "QRY" as query ID for the draft sequence and "SBJ" as subject ID for the reference genome. In the folder "Test" are present data for testing the tool.

## Installation
This tool uses MUSCLE Sequence alignment tool to perform the alignment thus the installation is required (We used MUSCLE v3.8.1551 to test the tool).
The last version can be installed as follows:

```
conda install muscle
```
Or you can download it from this site and install it manually: http://drive5.com/muscle/downloads.htm

No other installation steps are required. You need just to run the python script of the program. 
Python 3.6.0 (https://www.python.org/downloads/release/python-360/) or higher is required but also older version of python3 should work.

The following python packages are required:
- os
- re
- subprocess
- sys
- itertools
- pandas
- numpy
- Bio
- argparse

For an easier installation we suggest to use conda like this:

1- create a conda environment:

```conda create --name envname```

2- activate the environment:

```conda activate envname```

3- install the following packadges:

```conda config --add channels defaults```

```conda config --add channels bioconda```

```conda config --add channels conda-forge```

```conda install pandas```

```conda install -c anaconda numpy```

```conda install -c conda-forge biopython```

```conda install -c anaconda argparse```



## Parameters

  * #### Mandatory
| Name  | Example value | Description     |
|-------|---------------|-----------------|
| -q --input1| query.fa | draft sequence in Fasta format |
| -s --input2| subject.fa | reference genome in Fasta format |
| -f --input3| reference_contigs.fa | a Fasta file containing all the reads or sequences to use as reference database for the polishing of the draft sequence|

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| -h   | Display help |

## Usage 

```
python DraftPolisher.py -q query.fa -s subject.fa -f reads.fa
```

```
python DraftPolisher_cov.py -q query.fa -s subject.fa -f reads.fa
```


## Main Output

| Type      | Description     |
  |-----------|---------------|
  | consensus.fa    | the polished sequence in Fasta format |

## Other outputs

Several files are generated during the process and they can be used to evaluate the error rate, error correction rate and in general to have different checkpoints of the process. At the end of the process will be asked if you want to discard these process files or not.

## Contributions

| Name      | Email | Description     |
|-----------|---------------|-----------------|
  | Rosario Nicola Brancaccio | rosariobrancaccio@yahoo.it | Developer to contact for support |
  | Massimo Tommasino | tommasinom@iarc.fr
  | Tarik Gheit | gheitt@iartc.fr
  
## Versioning

Versions: 1.0, 1.0.cov


## Authors

**Rosario Nicola Brancaccio** - (https://github.com/BioRB/)

## License
[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)


## Acknowledgments
I would like to thanks Vincenzo Minissale and Sergey Senkin for their help in coding.

## References
