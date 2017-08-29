#!/usr/env python
import pysam,sys
from Alignment import Alignment
from Cigar import Cigar
from Read import Read
import numpy as np
from operator import itemgetter
def unionize(reads,qname,qAln):
	if reads.get(qname)!= None:
		read = reads[qname]
		read.loadAlignments(qAln,qAln.strand)
		reads[qname]=read
	else:
		read=Read()
		read.label(qname)
		read.loadAlignments(qAln,qAln.strand)
		reads[qname]=read
	return reads
class Bam():
	def __init__(self,SV=None,Args=None):
		self.ifh=Args.ifh
		self.bam = pysam.AlignmentFile(self.ifh,'rb')
		self.reads=None
		self.clips=None
		self.medStart=None
		self.medEnd=None
		self.startClips=None # (minClip,maxClip)
		self.endClips=None   # (minClip,maxClip)
		self.verbose=Args.verbose
		self.hasIns=False
		READS={}
		CLIPS=[]
		if self.verbose==True: print "processing left breakpoint"
		for al in self.bam.fetch(str(SV.chrom),SV.leftCI[0]-Args.windowFlank,SV.leftCI[1]+Args.windowFlank):
			if al.cigarstring==None or len(al.get_reference_positions())==0: continue
			qAln, qCig = Alignment(al), Cigar(al.cigarstring)
			if SV.svtype=='INS': qAln.queryPos(qCig)
			qAln.setClips(qCig,al.reference_start,al.reference_end,SV.svtype,SV.leftCI,SV.rightCI)
			if SV.svtype!='INV':qAln.cigarSV(qCig.cig,al.reference_start,SV,Args.minLen,Args.minOvr)
			CLIPS.append((qAln.startClip,qAln.endClip))
			READS=unionize(READS,al.query_name,qAln)
		if self.verbose==True: print "left breakpoint complete"
		if self.verbose==True: print "processing right breakpoint"
		for al in self.bam.fetch(str(SV.chrom),SV.rightCI[0]-Args.windowFlank,SV.rightCI[1]+Args.windowFlank):
			if al.cigarstring==None or len(al.get_reference_positions())==0: continue
			qAln, qCig = Alignment(al), Cigar(al.cigarstring)
			if SV.svtype=='INS': qAln.queryPos(qCig)
			qAln.setClips(qCig,al.reference_start,al.reference_end,SV.svtype,SV.leftCI,SV.rightCI)	
			if SV.svtype!='INV':qAln.cigarSV(qCig.cig,al.reference_start,SV,Args.minLen,Args.minOvr)
			CLIPS.append((qAln.startClip,qAln.endClip))
			READS=unionize(READS,al.query_name,qAln)
		if self.verbose==True: print "right breakpoint complete"
		self.clips=CLIPS
		self.reads=READS
		if SV.svtype=='INS': self.insertion(Args)
		self.pixelPrep(SV,Args)
	def insertion(self,Args=None):
		if self.verbose==True: print "processing insertions"
		for name in self.reads:
			qPos,qGaps=[],[]
			if len(self.reads[name].alignments) < 1: continue
			for Aln in self.reads[name].alignments: qPos.append((Aln.qStart,Aln.qEnd,Aln.strand,Aln.insertion))
			alns = sorted(list(set(qPos)),key=itemgetter(0,1))
			for i in range(len(alns)-1):
				qGap = alns[i+1][0]-alns[i][1]
				if alns[i][2]==alns[i+1][2]:
					if qGap >= Args.minLen: qGaps.append(qGap)
				if alns[i][3] != None: qGaps.append(alns[i][3])
				if alns[i+1][3] != None: qGaps.append(alns[i+1][3])
			if len(qGaps)>0: 
				self.hasIns,self.reads[name].insertion=True,max(qGaps)
		if self.verbose==True: print "insertion processing complete"
	def pixelPrep(self,SV=None,Args=None):
		if self.verbose==True: print "painting ..."
		self.medianClip(SV)
		self.assignClips()
		if self.medStart==self.medEnd: self.medEnd=self.medStart+1
		for name in self.reads: self.reads[name].pixelPrep(Args.maxMapq,self.medStart,self.medEnd,SV.svtype,Args.minInv)
	def medianClip(self,SV):
		start = [x[0] for x in self.clips if x[0] != None]
		end= [x[1] for x in self.clips if x[1] != None]
		if SV.svtype!='INS':
			if len(start)>0: self.medStart, self.startClips = int(np.median(start)), (min(start),max(start))
			else: self.medStart, self.startClips =SV.start , (SV.start,SV.start)
			if len(end)>0: self.medEnd, self.endClips =int(np.median(end)), (min(end),max(end))
			else: self.medEnd, self.endClips =SV.end, (SV.end,SV.end)
		else: 
			start=sorted(start)
			if len(start)>0: self.medStart,self.medEnd, self.startClips, self.endClips = start[0],start[-1],(min(start),max(start)),(min(start),max(start))
			else: self.medStart,self.medEnd, self.startClips,self.endClips=SV.start,SV.end,(SV.start,SV.start),(SV.end,SV.end)
		if self.endClips[0]<= self.startClips[1]:
			tmp = self.endClips[0]
			self.endClips[0],self.startClips[1]=self.startClips[1],tmp
			if self.endClips[0]==self.startClips[1]: self.endClips[0]+=1
	def assignClips(self):
		for name in self.reads:
			Read=self.reads[name]
			start,end=[],[]
			for Aln in Read.alignments:
				if Aln.startClip!=None: 
					start.append(abs(Aln.startClip-self.medStart))
					Read.clips+=1
				if Aln.endClip!=None: 
					Read.clips+=1
					end.append(abs(Aln.endClip-self.medEnd))
			if len(start)>0: Read.startClip=sorted(start).pop(0)	
			if len(end)>0: Read.endClip=sorted(end).pop(0)
			score=sum([x for x in [Read.startClip,Read.endClip] if x!= None])
			if Read.startClip!=None or Read.endClip!=None: Read.score=score
