#!/usr/env python
from pybedtools import BedTool as Bed
from Cigar import Cigar
import numpy as np
def overlap(s1,e1,s2,e2):
	s=sorted([s1,s2])
	e=sorted([e1,e2])
	ovr = e[0]-s[1]+1
	o=sorted([float(ovr)/(e2-s2+1),float(ovr)/(e1-s1+1)])
	return o[0]
class Alignment():
	def __init__(self,al=None):
		strand='+'
		if al.is_reverse: strand='-'
		self.pos=al.get_reference_positions()
		self.strand=strand
		self.mapq=int(al.mapping_quality)
		self.startClip=None
		self.endClip=None
		self.qStart=None
		self.qEnd=None
		self.insertion=None
	def queryPos(self,cig): 
		cig.qPos()
		self.qStart,self.qEnd = cig.qStart,cig.qEnd
	def setClips(self,cig,leftPos,rightPos,leftCI,rightCI): 
		if cig.leftClip==True and leftCI[0]<=leftPos<=leftCI[1]: self.startClip=leftPos
		elif cig.leftClip==True and rightCI[0]<=leftPos<=rightCI[1]:self.endClip=leftPos
		if cig.rightClip==True and leftCI[0]<=rightPos<=leftCI[1]: self.startClip=rightPos
		elif cig.rightClip==True and rightCI[0]<=rightPos<=rightCI[1]: self.endClip=rightPos
	def cigarSV(self,cig,left,breakStart,breakEnd,svtype,leftCI,rightCI,minLen):
		ind=0
		for (flg,leng) in cig:
			if (flg==1 or flg==2) and leng > (minLen-1):
				s1 = left+ind
				e1 = left+ind+leng-1
				if svtype=='DUP': 
					s2 = left+ind-leng
					e2 = s1-1
					ovr1 = overlap(breakStart,breakEnd,s1,e1)
					ovr2 = overlap(breakStart,breakEnd,s2,e2)
					if ovr1 < ovr2: s1,e1=s2,e2
				if leftCI[0]<=s1<=leftCI[1]: self.startClip=s1
				if rightCI[0]<=e1<=rightCI[1]: self.endClip=e1
				if svtype=='INS':
					if leftCI[0]<=s1<=leftCI[1]: self.startClip, self.insertion=s1,leng
					if rightCI[0]<=e1<=rightCI[1]: self.endClip,self.insertion=e1,leng
			if flg==0 or flg==2 or flg==3 or flg==7 or flg==8: ind+=leng	
