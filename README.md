# SVanGogh
Pixelate SVs!

## Install

1. `git clone --recursive https://github.com/dantaki/SVanGogh`
2. `cd SVanGogh && pip install dist/SVanGogh-0.1.tar.gz` 

#### Requires:
* Pillow (PIL)
* pybedtools
* pysam
* scipy
* tqdm 

## Usage

```
$ svangogh --help

usage: svangogh [-h] -i I [-r R] [-b B] [-v V] [-t T] [-ci] [-w W] [-c C]
                [-f F] [-n N] [-m M] [-s S] [-hs HS] [-ws WS] [-V] [-P] [-o O]

svangogh     --paint SV breakpoints--

optional arguments:
  -h, --help    show this help message and exit

required arguments:
  -i I          BAM file

SV arguments:
  -r R          Breakpoint <chr:start-end>
  -b B          Breakpoint BED file, tab-delimited, <chr start end type>
  -v V          VCF file
  -t T          SV type <DEL|DUP|INV|INS>
  -ci           Search for clips within confidence intervals. Requires VCF. Overrides <-c>
  -w W          Flanking bp to search for supporting reads. [100]
  -c C          Maximum clipped distance to breakpoint. [50]

SV painting arguments:
  -f F          Flanking bp to paint. [20]
  -n N          Maximum number of reads to paint. [10]
  -m M          Maximum MAPQ. [Maximum in sample of reads]
  -s S          Scaling multiplier. Adjust the scaled image size. [5]
  -hs HS        Height scaling multiplier. Adjust the height of the scaled image size. [5]
  -ws WS        Width scaling multiplier. Adjust the width of the scaled image size [5]. 

Optional arguments:
  -V            Verbose.
  -P            Display a progress bar.
  -o O, -out O  output
```

## Output

svangogh creates three output files for each SV: 

* data file containing pixels with RGB values 
* PNG image of the SV, scaled to 800x300
* PNG image of the SV, unscaled

The prefix of the output file defined by the `<-o output>` argument

## Methodology

:construction::warning: **WIP** :warning::construction:

### Subsampling Reads 
To make images consistent for learning, the default maximum number of reads (rows) to pixelate is 10.  

Reads are then ordered, first printing out supporting reads. A more detailed list of ordering is found below:

#### Deletions, Duplications

1. Reads with alignments on the same strand. Minimum sum of the difference between both clipped positions to the median clipped positions and MAPQ (maximum 60). 
  * the read has two clipped positions
2. Alignments on the same strand. Minimum sum of the difference between a clipped position to the median and MAPQ.
  * the read has one clipped position
3. Alignments on the same strand. Minimum sum of he MAPQ.
  * the read has no clipped positions
4. Alignments on *different* strands. Minimum sum in 1.
  * the read has alignments on different strands, which for deletions and duplications are more likely to be erroreous than valid. Two clipped positions.
5. Similar to 2. but with alignments on different strands.
6. Similar to 3. but with alignments on different strands.

##### Inversions

1. Alignments on different strands. Minimum sum of the differences between both clipped positions to the median and MAPQ (maximum 60). 
  * For inversions, at least one alignment should map to an opposite strand. Two clipped positions
2. Alignments on different strands. Minimum sum of one clipped positions and MAPQ.
3. Alignments on the same strand. Minimum sum of the difference between a clipped position to the median and MAPQ
  * Sometimes the inverted sequence portion of the read is unmapped, but has an informative clipped alignment
4. Alignments on the same strand. Munimum sum of the MAPQ.
  * No clips, properly aligned reads


--- 

## Contact

This is an alpha release of svangogh. Use at your own risk:exclamation:

Author: Danny Antaki `dantaki@ucsd.edu`
