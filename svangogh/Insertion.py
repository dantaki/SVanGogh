#!/usr/env python
from operator import itemgetter
class Insertion():
	def __init__(self):
		self.size=0
		self.leftClip=None
		self.rightClip=None
	def findInsertions(self,reads):
		leftClips=[]
		rightClips=[]
		size=[]
		for name in reads:
			qPos=[]
			qGaps=[]
			if len(reads[name].alignments) < 1: continue
			for Aln in reads[name].alignments:
				qPos.append((Aln.qStart,Aln.qEnd,Aln.strand,Aln.leftClip,Aln.rightClip))
			alns = sorted(list(set(qPos)),key=itemgetter(0,1))
			for i in range(len(alns)-1):
				qGap = alns[i+1][0]-alns[i][1]
				qGaps.append(qGap)
				if alns[i][2]==alns[i+1][2]:
					if qGap >= 20:
						size.append(qGap)
						if alns[i][3] != None: leftClips.append(alns[i][3])
						if alns[i+1][3] != None: leftClips.append(alns[i+1][3])
						if alns[i][4] != None: rightClips.append(alns[i][4])
						if alns[i+1][4] != None: rightClips.append(alns[i+1][4])
			if len(qGaps)>0: reads[name].insertion=max(qGaps)
		if len(size)>0: self.size = max(size)
		if len(leftClips)>0: self.leftClip=max(leftClips)
		if len(rightClips)>0: self.rightClip=min(rightClips)
