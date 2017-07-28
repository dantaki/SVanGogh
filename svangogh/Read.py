#!/usr/env python
import pysam
class Read():
	def __init__(self):
		self.name=None
		self.alignments=[]
		self.left=0
		self.right=0
		self.insertion=None
	def label(self,name): self.name=name
	def loadAlignments(self,aln): self.alignments.append(aln)
	def readPosition(self):
		self.left=self.alignments[0].pos[0]
		self.right=self.alignments[0].pos[-1]
		for Aln in self.alignments:
			if Aln.pos[0] < self.left: self.left=Aln.pos[0]
			if Aln.pos[-1] > self.right: self.right=Aln.pos[-1]
