#!/usr/env python
import numpy as np
class Read():
	def __init__(self):
		self.name=None
		self.alignments=[]
		self.left=0
		self.right=0
		self.insertion=None
		self.strandPix=None
		self.sameStrand=True
		self.startClip=None
		self.endClip=None
		self.mapq=None
		self.score=None
	def label(self,name): self.name=name
	def loadAlignments(self,aln): self.alignments.append(aln)
	def pixelPrep(self,MAX):
		self.readPosition()
		self.countStrands()
		self.medianMapq(MAX)
		self.scoreAlignments()
	def readPosition(self):
		self.left=self.alignments[0].pos[0]
		self.right=self.alignments[0].pos[-1]
		for Aln in self.alignments:
			if Aln.pos[0] < self.left: self.left=Aln.pos[0]
			if Aln.pos[-1] > self.right: self.right=Aln.pos[-1]
	def countStrands(self):
		fow,rev=0,0
		for Aln in self.alignments:
			if Aln.strand=='+': fow+=1
			if Aln.strand=='-': rev+=1
		self.strandPix='+'
		if rev > fow: self.strandPix='-'
		if fow>0 and rev>0: self.sameStrand=False
	def medianMapq(self,MAX):
		q=[]
		for Aln in self.alignments: q.append(Aln.mapq)
		mapq=np.median(q)
		if mapq>MAX:mapq=MAX
		self.mapq=abs(MAX-mapq)
	def scoreAlignments(self): 
		score=sum([x for x in [self.startClip,self.endClip] if x!= None])
		if self.startClip!=None or self.endClip!=None: self.score=score
