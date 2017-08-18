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

---------------------------------------------------------
8""""8  88   8                8""""8                   
8       88   8  eeeee  eeeee  8    "  eeeee  eeeee  e   e 
8eeeee  88  e8  8   8  8   8  8e      8  88  8   8  8   8 
    88  "8  8   8eee8  8e  8  88  ee  8   8  8e     8eee8 
e   88   8  8   88  8  88  8  88   8  8   8  88 "8  88  8 
8eee88   8ee8   88  8  88  8  88eee8  8eee8  88ee8  88  8

pixelate SVs         Author:Danny Antaki dantaki@ucsd.edu
---------------------------------------------------------

Usage: svangogh (-i BAM) [-r REGION] [-t SVTYPE] [-v VCF] [-b BED] 
                         [--ci] [-w WINDOW] [-c CLIP] [--min-inv MININV] [--min-indel INDEL] 
                         [-f FLANK] [--max-reads MAXREADS] [--max-mapq MAXMAPQ] [--min-sr MINSR] 
                         [-s SCALE] [--hs HSCALE] [--ws WSCALE] 
                         [-V --verbose] [-P] [-o OUT] [-h --help]

  -h --help               show this help message and exit
  --version               print the version number

Required Arguments:
  -i BAM                  BAM file

SV Arguments:  
  -r REGION               breakpoint [chr:start-end]
  -t SVTYPE               SV type. Required if -r is defined. [DEL|DUP|INV|INS]
  -v VCF                  VCF file
  -b BED                  BED file

SV Options:
  --ci                    search for clips within confidence intervals, requires -v, overrides -c
  -w WINDOW               flanking bp to search for supporting reads [default: 100]
  -c CLIP                 maximum distance of clipped position to breakpoint [default: 50]
  --min-inv MININV        minimum overlap of the alignment to inversion, flags supporting reads [default: 0.5]
  --min-indel INDEL       minimum indel size [default: 7]

Pixelating Options:  
  -f FLANK                flanking bp to paint [default: 20]
  --max-reads MAXREADS    maximum number of reads to paint [default: 10]
  --max-mapq MAXMAPQ      maximum MAPQ [Maximum in sample of reads]
  --min-sr MINSR          minimum number of supporting reads [default: 0]
  -s SCALE                scaling multiplier, adjust the scaled image size [default: 5]
  --hs HSCALE             height scaling multiplier, adjust scaled image height [default: 5]
  --ws WSCALE             width scaling multiplier, adjust scaled image width [default: 5]

Options:  
  -V --verbose            verbose mode
  -P                      display a progress bar
  -o OUT                  output prefix [default: svangogh]
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
