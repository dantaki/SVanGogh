#!/usr/env python
import sys
def qc(c,s,e,t):
	skip=False
	if t!='DEL' and t!='DUP' and t!='INV' and t!='INS':
		skip=True
		sys.stderr.write('WARNING:{} not an acceptable SV type. Accepted SV types: DEL, DUP, INV, INS.\nError found here: {}:{}-{} Skipping ...\n'.format(t,c,s,e))
	if (t=='DEL' or t=='DUP' or t=='INV') and e<=s:
		skip=True
		sys.stderr.write('WARNING:End position less than the start position.\nError found here: {}:{}-{} {} Skipping ...\n'.format(c,s,e,t))
	return skip
class SV():
	def __init__(self,c=None,s=None,e=None,cl=None,leftCI=(None,None),rightCI=(None,None)):
		self.chrom=c
		self.start=int(s)
		self.end=int(e)
		if self.end == 0: self.end+1
		self.svtype=cl
		self.leftCI=(self.start+leftCI[0],self.start+leftCI[1])
		self.rightCI=(self.end+rightCI[0],self.end+rightCI[1])
		self.region='{}:{}-{}'.format(self.chrom,self.start,self.end)
		self.qc=0 #Fail
		if qc(self.chrom,self.start,self.end,self.svtype)==False: self.qc=1
