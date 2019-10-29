# DraftPolisher

**DraftPolisher** is a tool for fast polishing of draft circular genomes.

## Description

DraftPolisher can produce an improved consensus sequence for a draft circular genome assembly. 

## Prerequisites

DraftPolisher is designed for the polishing of draft circular sequences. The number of mismatches may not exceed 26 bp in a window of 34 bp before or after any matching base.This tool is intended for the polishing of draft circular genomes (or any circular sequence) that are not so far from the reference sequence. The presence of large gaps will cause an error in the process. In this case, we suggest to perform a preliminal alignment and remove the bigger gaps before to precess it with this tool. The tool will only evaluate the mismatches and gaps present and will not affect the portion of the sequence where no mismatches or gaps have been identified. 

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
| -q --input1| query.fa | query sequence in Fasta format |
| -s --input2| subject.fa | subject sequence in Fasta format |
| -f --input3| reference_contigs.fa | a Fasta file containing all the reads or sequences that you want to use for the polishing of the draft sequence|

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

| Type      | Description     |
  |-----------|---------------|

## Contributions

| Name      | Email | Description     |
|-----------|---------------|-----------------|
  | Rosario Nicola Brancaccio | rosariobrancaccio@yahoo.it | Developer to contact for support |
  
## Versioning

Version 1.0

## Authors

**Rosario Nicola Brancaccio** - (https://github.com/BioRB/)

## License
[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)


## Acknowledgments
I would like to thanks Sergey Senkin  and Vincenzo Minissale for their help in coding.

## References
