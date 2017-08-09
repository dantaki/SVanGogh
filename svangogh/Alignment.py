#!/usr/env python
from pybedtools import BedTool as Bed
from Cigar import Cigar
import numpy as np
def isOverlapping(s1,e1,s2,e2): 
	s=sorted([s1,s2])
	e=sorted([e1,e2])
	if s[1] < e[0]: return True
	else: return False 
def overlap(s1,e1,s2,e2):
	s=sorted([s1,s2])
	e=sorted([e1,e2])
	ovr = e[0]-s[1]+1
	o=sorted([float(ovr)/(e2-s2+1),float(ovr)/(e1-s1+1)])
	return o[0]
def getBreak(s1,e1,s2,e2,svtype):
	a = sorted([s1,e1,s2,e2])
	if svtype=='DEL' or svtype=='INV' or (svtype=='DUP' and isOverlapping(s1,e2,s2,e2)==True): return a[1],a[2]
	elif svtype=='DUP' and isOverlapping(s1,e2,s2,e2)==False: return a[0],a[3]
class Alignment():
	def __init__(self):
		self.pos=[None]
		self.strand=None
		self.mapq=None
		self.leftClip=None
		self.rightClip=None
		self.isRight=False #isRight is True when the alignment is on the right hand-side of the breakpoint
		self.qStart=None
		self.qEnd=None
	def refPos(self,cig,left,posit):
		pos=[]
		index=left
		for (flg,leng) in cig.cig:
			if flg==0 or flg==7 or flg==8:
				for p in range(leng):
					pos.append(index+p)
			if flg==0 or flg==2 or flg==3 or flg==7 or flg==8: index+=leng
		self.pos=posit
		#self.pos=pos
	def queryPos(self,cig): 
		cig.qPos()
		self.qStart,self.qEnd = cig.qStart,cig.qEnd
	def orient(self,strand,sidebool): 
		self.strand=strand
		if sidebool==True: self.isRight=True
	def quality(self,mapq):  self.mapq=int(mapq)
	def setClips(self,cig,leftPos,rightPos,leftCI,rightCI,svtype): 
		leftClipPos=None
		rightClipPos=None
		if svtype=='DEL' or svtype=='INV':
			if cig.leftClip!=None and rightCI[0]<=leftPos<=rightCI[1]: leftClipPos=leftPos
			if cig.rightClip!=None and leftCI[0]<=rightPos<=leftCI[1]: rightClipPos=rightPos
		if svtype=='DUP':
			if cig.leftClip!=None and leftCI[0]<=leftPos<=leftCI[1]: leftClipPos=leftPos
			if cig.rightClip!=None and rightCI[0]<=rightPos<=rightCI[1]: rightClipPos=rightPos
		if svtype=='INS':
			if cig.leftClip!=None and leftCI[0]<=leftPos<=leftCI[1]: leftClipPos=leftPos
			if cig.rightClip!=None and leftCI[0]<=rightPos<=leftCI[1]: rightClipPos=rightPos
		self.leftClip,self.rightClip=leftClipPos,rightClipPos
	def chimericClip(self,pChrom,pStart,pEnd,chromBreak,leftBreak,rightBreak,svtype,salns,leftCI,rightCI):
		if not salns[-1].startswith('SA'): del salns[-1]
		for saln in salns:
			salnList = saln.split(',')
			sRef,sLeft,sStrand,sCigar= salnList[0:4]
			sCig = Cigar(sCigar)
			sLeft=int(sLeft)-1
			sRight = int(sLeft)+sCig.aLen()
			if pChrom == sRef and self.strand == sStrand:
				brkS,brkE = getBreak(pStart,pEnd,sLeft,sRight,svtype)
				pLeft=True
				if pStart > sLeft: pLeft=False
				if Bed('{} {} {}'.format(chromBreak,leftBreak,rightBreak),from_string=True).intersect(Bed('{} {} {}'.format(sRef,brkS,brkE),from_string=True),f=0.90,F=0.90).count() >0:
					if svtype=='DUP' and isOverlapping(pStart,pEnd,sLeft,sRight)==False: 
						if pLeft==True and leftCI[0]<=brkS<=leftCI[1]: self.leftClip=brkS
						if pLeft==False and rightCI[0]<=brkE<=rightCI[1]: self.rightClip=brkE
					elif svtype=='DUP' and isOverlapping(pStart,pEnd,sLeft,sRight)==True: 
						if pLeft==False and leftCI[0]<=brkS<=leftCI[1]: self.leftClip=brkS
						if pLeft==True and rightCI[0]<=brkE<=rightCI[1]: self.rightClip=brkE
					elif svtype=='DEL' or svtype=='INV':
						if pLeft==True and leftCI[0]<=brkS<=leftCI[1]: self.rightClip=brkS
						if pLeft==False and rightCI[0]<=brkE<=rightCI[1]: self.leftClip=brkE
	def cigarSV(self,cig,left,breakStart,breakEnd,svtype,leftCI,rightCI):
		overlaps=[]
		ind=0
		for (flg,leng) in cig:
			if (flg==1 or flg==2) and leng > 100:
				s1 = left+ind
				e1 = left+ind+leng
				s2 = left+ind-leng
				e2 = s1
				ovr=None
				bs = s1
				be = e1
				if svtype=='DEL': ovr=overlap(breakStart,breakEnd,s1,e1)
				if svtype=='DUP': 
					ovr1 = overlap(breakStart,breakEnd,s1,e1)
					ovr2 = overlap(breakStart,breakEnd,s2,e2)
					if ovr1 > ovr2: ovr=ovr1
					else: 
						ovr=ovr2
						bs=s2
						be=e2
				if ovr < 0: continue
				if ovr >= 0.90:
					if leftCI[0]<=bs<=leftCI[1]: self.leftClip=bs
					if rightCI[0]<=be<=rightCI[1]: self.rightClip=be
			if flg==0 or flg==2 or flg==3 or flg==7 or flg==8: ind+=leng	
