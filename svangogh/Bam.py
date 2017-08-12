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
		self.medStart=None
		self.medEnd=None
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
			qAln.setClips(qCig,al.reference_start,al.reference_end,SV.leftCI,SV.rightCI,)
			if SV.svtype=='DEL' or SV.svtype=='DUP':qAln.cigarSV(qCig.cig,al.reference_start,SV.start,SV.end,SV.svtype,SV.leftCI,SV.rightCI)
			#####################
			CLIPS.append((qAln.startClip,qAln.endClip))
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
			qAln.setClips(qCig,al.reference_start,al.reference_end,SV.leftCI,SV.rightCI)
			if SV.svtype=='DEL' or SV.svtype=='DUP':qAln.cigarSV(qCig.cig,al.reference_start,SV.start,SV.end,SV.svtype,SV.leftCI,SV.rightCI)
			#####################
			CLIPS.append((qAln.startClip,qAln.endClip))
			READS=unionize(READS,al.query_name,qAln)
		self.clips=CLIPS
		self.reads=READS
	def pixelPrep(self,SV):
		self.medianClip(SV)
		self.assignClips()
	def medianClip(self,SV):
		start = [x[0] for x in self.clips if x[0] != None]
		end= [x[1] for x in self.clips if x[1] != None]
		if len(start)>0: self.medStart=int(np.median(start))
		else: self.medStart=SV.start
		if len(end)>0: self.medEnd=int(np.median(end))
		else: self.medEnd=SV.end
	def assignClips(self):
		for name in self.reads:
			Read=self.reads[name]
			start,end=[],[]
			for Aln in Read.alignments:
				if Aln.startClip!=None: start.append(abs(Aln.startClip-self.medStart))
				if Aln.endClip!=None: end.append(abs(Aln.endClip-self.medEnd))
			if len(start)>0: Read.startClip=sorted(start).pop(0)	
			if len(end)>0: Read.endClip=sorted(end).pop(0)
