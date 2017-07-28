# SVanGogh
Pixelate SVs!

## Install

`pip install https://github.com/dantaki/SVanGogh/releases/tag/v0.1`

#### Requires:

* pysam
* pybedtools
* scipy

## Usage

```
$ svangogh --help

optional arguments:
  -h, --help    show this help message and exit

required arguments:
  -i I          BAM file

SV painting arguments:
  -r R          Breakpoint <chr:start-end>
  -b B          Breakpoint BED file, tab-delimited, <chr start end type>
  -v V          VCF file
  -t T          SV type <DEL|DUP|INV|INS>
  -c C          Maximum clipped distance to breakpoint. Default:50
  -f F          Flanking bp to paint. Default 20

Optional arguments:
  -ci           Search for clips within confidence intervals. Requires VCF. Overrides <-c>
  -w W          Flanking bp to search for supporting reads. Default 100
  -o O, -out O  output

```

## Output

svangogh creates three output files for each SV: 

* data file containing pixels with RGB values 
* PNG image of the SV, scaled to 800x300
* PNG image of the SV, unscaled

The prefix of the output file defined by the `<-o output>` argument

## Contact

This is an alpha release of svangogh. Use at your own risk:exclamation:

Author: Danny Antaki `dantaki@ucsd.edu`
