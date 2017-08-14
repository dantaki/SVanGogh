#!/usr/env python
import os,argparse,sys,subprocess
from argparse import RawTextHelpFormatter
from SV import SV
import pysam
class Arguments():
	def __init__(self):
		splash='svangogh     --paint SV breakpoints--\n'
		parser = argparse.ArgumentParser(description=splash,formatter_class=RawTextHelpFormatter)
		reqArgs = parser.add_argument_group('required arguments')
		svArgs = parser.add_argument_group('SV arguments')
		pixArgs = parser.add_argument_group('SV painting arguments')
		optArgs = parser.add_argument_group('Optional arguments')
		reqArgs.add_argument('-i', help='BAM file',required=True,type=str)
		svArgs.add_argument('-r', help='Breakpoint <chr:start-end>',required=False,type=str,default=None)
		svArgs.add_argument('-b', help='Breakpoint BED file, tab-delimited, <chr start end type>',required=False,type=str,default=None)
		svArgs.add_argument('-v', help='VCF file',required=False,type=str,default=None)
		svArgs.add_argument('-t', help='SV type <DEL|DUP|INV|INS>',required=False,type=str) 
		svArgs.add_argument('-ci',help='Search for clips within confidence intervals. Requires VCF. Overrides <-c>',required=False,default=False,action='store_true')
		svArgs.add_argument('-w', help='Flanking bp to search for supporting reads. [100]',required=False,type=int,default=100)
		svArgs.add_argument('-c', help='Maximum clipped distance to breakpoint. [50]',required=False,type=int,default=50)
		pixArgs.add_argument('-f', help='Flanking bp to paint. [20]',required=False,type=int,default=20)
		pixArgs.add_argument('-n', help='Maximum number of reads to paint. [10]',required=False,type=int,default=10)
		pixArgs.add_argument('-m', help='Maximum MAPQ. [Maximum in sample of reads]',required=False,type=int,default=None)
		pixArgs.add_argument('-s', help='Scaling multiplier. Adjust the scaled image size. [5]',required=False,type=int,default=5)
		pixArgs.add_argument('-hs', help='Height scaling multiplier. Adjust the height of the scaled image size. [5]',required=False,type=int,default=None)
		pixArgs.add_argument('-ws', help='Width scaling multiplier. Adjust the width of the scaled image size [5]. ',required=False,type=int,default=None)
		optArgs.add_argument('-V', help='Verbose. Display a progress bar.',required=False,default=False,action='store_true')
		optArgs.add_argument('-o','-out', help='output',required=False,default="breakPainter",type=str)
		args = parser.parse_args()
		self.ifh = args.i
		region = args.r
		bed=args.b
		vcf=args.v
		self.maxMapq=args.m
		self.breakType=args.t
		self.maxClip = args.c 
		self.maxFlank = args.f
		self.maxReads=args.n
		self.scaling = args.s
		self.hscaling = args.hs
		self.wscaling = args.ws
		self.verbose= args.V
		self.ci, self.windowFlank, self.ofh = args.ci, args.w,args.o
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
		if region != None and bed==None and self.breakType==None and vcf==None:
			sys.stderr.write('FATAL ERROR: Please define the SV type <-t <DEL|DUP|INV|INS>\n')
			sys.exit(1)
		regions=[]
		if bed !=None:
			with open(bed,'r') as f:
				for l in f:
					a = l.rstrip('\n').split('\t')
					sv = SV(a[0],a[1],a[2],a[3],(-1*self.maxClip,self.maxClip),(-1*self.maxClip,self.maxClip))
					if sv.qc==1: regions.append(sv)	
		if region!=None:
			c,p = region.split(':')
			s,e = p.split('-')
			sv = SV(c,s,e,self.breakType,(-1*self.maxClip,self.maxClip),(-1*self.maxClip,self.maxClip))
			if sv.qc==1: regions.append(sv)		
		if vcf != None:
			with open(vcf,'r') as f:
				for l in f:
			       		if l.startswith('#'): continue
					r = l.rstrip('\n').split('\t')
					c=str(r[0])
					s=int(r[1])
					s-=1
					e=0
					svtype=None
					leftCI, rightCI  = (-1*self.maxClip, self.maxClip), (-1*self.maxClip,self.maxClip)
					for i in r[7].split(';'):
						if 'SVTYPE=' in i: svtype=i.replace('SVTYPE=','')
						if i.startswith('END='): e=int(i.replace('END=',''))
						if self.ci==True and i.startswith('CIPOS='): leftCI=tuple(map(int,i.replace('CIPOS=','').split(',')))
						if self.ci==True and i.startswith('CIEND='): rightCI=tuple(map(int,i.replace('CIEND=','').split(',')))
					sv = SV(c,s,e,svtype,leftCI,rightCI)
					if sv.qc==1: regions.append(sv)
		self.regions=regions	
