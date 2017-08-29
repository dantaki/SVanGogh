#!/usr/env python
import sys
class SV():
	def __init__(self,c=None,s=None,e=None,cl=None,leftCI=(None,None),rightCI=(None,None)):
		ins=['INS','SVA','LINE1','ALU','HERV']
		self.chrom=c
		self.start=int(s)
		self.end=int(e)
		self.svtype=cl
		if self.svtype in ins: self.svtype='INS'
		self.leftCI=(self.start+leftCI[0],self.start+leftCI[1])
		self.rightCI=(self.end+rightCI[0],self.end+rightCI[1])
		self.region='{}:{}-{}'.format(self.chrom,self.start,self.end)
		self.qc=0 #Fail
		if self.qc()==False: self.qc=1
	def qc(self):
		skip=False
		if self.svtype!='DEL' and self.svtype!='DUP' and self.svtype!='INV' and self.svtype!='INS':
			skip=True
			sys.stderr.write('WARNING:{} not an acceptable SV type. Accepted SV types: DEL, DUP, INV, INS (INS,ALU,HERV,LINE1,SVA).\nError found here: {} Skipping ...\n'.format(self.svtype,self.region))
		if (self.svtype=='DEL' or self.svtype=='DUP' or self.svtype=='INV') and self.end<=self.start:
			skip=True
			sys.stderr.write('WARNING:End position less than the start position.\nError found here: {} {} Skipping ...\n'.format(self.region,self.svtype))
		return skip