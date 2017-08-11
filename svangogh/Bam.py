#!/usr/env python
import pysam
from Alignment import Alignment
from Cigar import Cigar
from Read import Read
import numpy as np
def unionize(reads,qname,qAln):
	if reads.get(qname)!= None:
		read = reads[qname]
		read.loadAlignments(qAln)
		reads[qname]=read
	else:
		read=Read()
		read.label(qname)
		read.loadAlignments(qAln)
		reads[qname]=read
	return reads
class Bam():
	def __init__(self,ifh):
		self.ifh=ifh
		self.bam = pysam.AlignmentFile(ifh,'rb')
		self.reads=None
		self.clips=None
		self.medLeftClip=None
		self.medRightClip=None
	def leftFlank(self,SV,Args):
		READS={}
		CLIPS=[]
		for al in self.bam.fetch(str(SV.chrom),SV.leftCI[0]-Args.windowFlank,SV.leftCI[1]+Args.windowFlank):
			if al.cigarstring==None or len(al.get_reference_positions())==0: continue
			qAln = Alignment()
			pStrand='+'
			if al.is_reverse: pStrand='-'
			qCig = Cigar(al.cigarstring)
			##### LOAD qALN #####
			qAln.refPos(qCig,al.reference_start,al.get_reference_positions())
			if SV.svtype=='INS': qAln.queryPos(qCig)
			qAln.orient(pStrand)
			qAln.quality(al.mapping_quality)
			qAln.setClips(qCig,al.reference_start,al.reference_end,SV.leftCI,SV.rightCI,SV.svtype)
			if SV.svtype=='DEL' or SV.svtype=='DUP':qAln.cigarSV(qCig.cig,al.reference_start,SV.start,SV.end,SV.svtype,SV.leftCI,SV.rightCI)
			#####################
			CLIPS.append((qAln.leftClip,qAln.rightClip))
			READS=unionize(READS,al.query_name,qAln)
		self.reads=READS
		self.clips=CLIPS
	def rightFlank(self,SV,Args):
		READS=self.reads
		CLIPS=self.clips
		for al in self.bam.fetch(str(SV.chrom),SV.rightCI[0]-Args.windowFlank,SV.rightCI[1]+Args.windowFlank):
			if al.cigarstring==None or len(al.get_reference_positions())==0: continue
			qAln = Alignment()
			pStrand='+'
			if al.is_reverse: pStrand='-'
			qCig = Cigar(al.cigarstring)
			##### LOAD qALN #####
			qAln.refPos(qCig,al.reference_start,al.get_reference_positions())
			qAln.orient(pStrand)
			qAln.quality(al.mapping_quality)
			qAln.setClips(qCig,al.reference_start,al.reference_end,SV.leftCI,SV.rightCI,SV.svtype)
			if SV.svtype=='DEL' or SV.svtype=='DUP':qAln.cigarSV(qCig.cig,al.reference_start,SV.start,SV.end,SV.svtype,SV.leftCI,SV.rightCI)
			#####################
			CLIPS.append((qAln.leftClip,qAln.rightClip))
			READS=unionize(READS,al.query_name,qAln)
		self.clips=CLIPS
		self.reads=READS
	def pixelPrep(self,SV,Args):
		self.maxMapq(Args)
		self.medianClip(SV)
		self.assignClips()
	def medianClip(self,SV):
		left = [x[0] for x in self.clips if x[0] != None]
		right= [x[1] for x in self.clips if x[1] != None]
		if len(left)>0: self.medLeftClip=int(np.median(left))
		else: self.medLeftClip=SV.start
		if len(right)>0: self.medRightClip=int(np.median(right))
		else: self.medRightClip=SV.end
	def assignClips(self):
		for name in self.reads:
			Read=self.reads[name]
			leftClip,rightClip=[],[]
			for Aln in Read.alignments:
				if Aln.leftClip!=None: leftClip.append(abs(Aln.leftClip-self.medLeftClip))
				if Aln.rightClip!=None: rightClip.append(abs(Aln.rightClip-self.medRightClip))
			if len(leftClip)>0: Read.leftClip=sorted(leftClip).pop(0)	
			if len(rightClip)>0: Read.rightClip=sorted(rightClip).pop(0)
	def maxMapq(self,Args):
		c=0;
		maxQ=0
		for al in self.bam.fetch(until_eof=True):
			c+=1
			if c>5000: break
			if al.mapping_quality>maxQ: maxQ=al.mapping_quality
		Args.maxMapq=maxQ	
