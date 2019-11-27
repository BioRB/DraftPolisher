# DraftPolisher

**DraftPolisher** is a tool for fast polishing of draft circular genomes.

## Description

DraftPolisher can produce an improved consensus sequence for a draft circular genome assembly. 

## Prerequisites

DraftPolisher is designed for the polishing of draft circular sequences. The number of nucleotide gaps may not exceed 26 bp in a window of 34 bp before or after any mismatch. The presence of large gaps will cause an error and block in the process. In this case, we suggest to perform a preliminary alignment and remove the bigger gaps before to precess it with this tool. The tool will only evaluate the mismatches and gaps present and will not affect the portion of the sequence where no mismatches or gaps have been identified. DraftPolisher_cov.py is an alternative release that takes into account also the coverage of the contigs produced in the upstream assembling step that you have presumably carried out before. DraftPolisher_cov.py was designed based on the SPAdes assembling output file format, where the coverage value is reported in the last part of the sequences IDs, so if you have in mind to use this version of the tool be careful about the formatting of your assembling output file (see SPAdes output format for more details). A general prerequisite for both the versions of the tool is the formatting of query and subject fasta files. It is mandatory to put "QRY" as query ID for the draft sequence and "SBJ" as subject ID for the reference genome. In the folder "Test" you can find data to test the tool.

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
  
## Versioning

Versions: 1.0, 1.0.cov


## Authors

**Rosario Nicola Brancaccio** - (https://github.com/BioRB/)

## License
[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)


## Acknowledgments
I would like to thanks Sergey Senkin  and Vincenzo Minissale for their help in coding.

## References
