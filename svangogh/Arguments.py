"""

---------------------------------------------------------
8\"\"\"\"8  88   8                8\"\"\"\"8                   
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

"""
import os,sys,subprocess,pysam
from docopt import docopt
from SV import SV
def tryInt(v,opt):
	try: int(v)
	except ValueError: 
		sys.stderr.write('FATAL ERROR: {} is not an integer in [{}]\n'.format(v,opt))
		sys.exit(1)
def tryFloat(v,opt):
	try: float(v)
	except ValueError: 
		sys.stderr.write('FATAL ERROR: {} is not a number in [{}]\n'.format(v,opt))
		sys.exit(1)
class Arguments():
	def __init__(self):
	    arg=docopt(__doc__,version='0.1')
	    for x in ['-w','-c','--min-indel','-f','--max-reads','--min-sr']: tryInt(arg[x],x)
	    for x in ['--min-inv','-s','--hs','--ws']: tryFloat(arg[x],x)
	    self.ifh = arg['-i']
	    region = arg['-r']
	    bed = arg['-b']
	    vcf = arg['-v']
	    self.breakType,self.ci,self.windowFlank,self.maxClip,self.minInv,self.minLen = arg['-t'],arg['--ci'],int(arg['-w']),int(arg['-c']),float(arg['--min-inv']),int(arg['--min-indel'])
	    self.maxFlank,self.maxReads,self.maxMapq,self.minSR = int(arg['-f']),int(arg['--max-reads']),arg['--max-mapq'],int(arg['--min-sr'])
	    self.scaling,self.hscaling,self.wscaling = float(arg['-s']),float(arg['--hs']),float(arg['--ws'])
	    self.verbose,self.progress,self.ofh = arg['--verbose'],arg['-P'],arg['-o']
	    if self.verbose==0: self.verbose=False 
	    else: self.verbose=True
	    if self.maxMapq==None:
	      c, maxQ=0,0
	      for al in pysam.AlignmentFile(self.ifh).fetch(until_eof=True):
	        c+=1
	        if c>5000: break
	        if al.mapping_quality>maxQ: maxQ=al.mapping_quality
	      self.maxMapq=maxQ
	      if self.verbose==True: print "Maximum MAPQ in sample:{}".format(self.maxMapq)
	    if region == None and bed ==None and vcf==None:
	      sys.stderr.write('FATAL ERROR: Please supply either a region <-r chr:start-end> or a BED file <-b> or a VCF file <-v>\n')
	      sys.exit(1)
	    if region != None and self.breakType==None:
	      sys.stderr.write('FATAL ERROR: Please define the SV type <-t <DEL|DUP|INV|INS>\n')
	      sys.exit(1)
	    regions=[]
	    if region!=None:
	      c,p = region.split(':')
	      s,e = p.split('-')
	      sv = SV(c,s,e,self.breakType,(-1*self.maxClip,self.maxClip),(-1*self.maxClip,self.maxClip))
	      if sv.qc==1: regions.append(sv)
	    if bed !=None:
	      with open(bed,'r') as f:
	        for l in f:
	          a = l.rstrip('\n').split('\t')
	          sv = SV(a[0],a[1],a[2],a[3],(-1*self.maxClip,self.maxClip),(-1*self.maxClip,self.maxClip))
	          if sv.qc==1: regions.append(sv)     
	    if vcf != None:
	      with open(vcf,'r') as f:
	        for l in f:
	          if l.startswith('#'): continue
	          r = l.rstrip('\n').split('\t')
	          c,s,e,svtype=str(r[0]),int(r[1]),0,None
	          leftCI, rightCI  = (-1*self.maxClip, self.maxClip), (-1*self.maxClip,self.maxClip)
	          s-=1
	          for i in r[7].split(';'):
	            if 'SVTYPE=' in i: svtype=i.replace('SVTYPE=','')
	            if i.startswith('END='): e=int(i.replace('END=',''))
	            if self.ci==True and i.startswith('CIPOS='): leftCI=tuple(map(int,i.replace('CIPOS=','').split(',')))
	            if self.ci==True and i.startswith('CIEND='): rightCI=tuple(map(int,i.replace('CIEND=','').split(',')))
	          sv = SV(c,s,e,svtype,leftCI,rightCI)
	          if sv.qc==1: regions.append(sv)
	    self.regions=regions