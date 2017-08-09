#!/usr/env python
import scipy.misc as smp
from operator import itemgetter
def pixelUnion(tmp,pix):
	unionPix=[]
	for i in xrange(len(tmp)):
		if tmp[i] == [0,0,0] and pix[i]==[0,0,0]:unionPix.append([0,0,0])
		if tmp[i] != [0,0,0] and pix[i]==[0,0,0]:unionPix.append(tmp[i])
		if tmp[i] == [0,0,0] and pix[i]!=[0,0,0]:unionPix.append(pix[i])
		if tmp[i] != [0,0,0] and pix[i]!=[0,0,0]:unionPix.append([max(tmp[i][0],pix[i][0]),max(tmp[i][1],pix[i][1]),max(tmp[i][2],pix[i][2])])
	return unionPix
def pixelUnionInsertion(tmp,pix,clip1,clip2):
	unionPix=[]
	for i in xrange(len(tmp)):
		if tmp[i] == [0,0,0] and pix[i]==[0,0,0]:unionPix.append([0,0,0])
		if tmp[i] != [0,0,0] and pix[i]==[0,0,0]:unionPix.append(tmp[i])
		if tmp[i] == [0,0,0] and pix[i]!=[0,0,0]:unionPix.append(pix[i])
		if tmp[i] != [0,0,0] and pix[i]!=[0,0,0]:
			if tmp[i][0]==clip1 and pix[i][0]==clip1: unionPix.append([clip2,max(tmp[i][1],pix[i][1]),max(tmp[i][2],pix[i][2])])
			else: unionPix.append([max(tmp[i][0],pix[i][0]),max(tmp[i][1],pix[i][1]),max(tmp[i][2],pix[i][2])])
	return unionPix
def appendOrder(i,o,t):
	for x in sorted(i, key=itemgetter(0)):
		if x[1] in t: o.append(x[1])
class Painter():
	def __init__(self,maxFlank,maxMapq):
		self.canvas=[]
		self.flank=maxFlank
		self.canvasLeftMin=None
		self.canvasRightMin=None
		self.samePix=[] # [ReadName]=[ readPix ]
		self.diffPix=[]
		self.twoClip=[]
		self.oneClip=[]
		self.mappedAln=[]
		self.insertionReads=[]
		self.pix={}
		self.readPix=[]
		self.order=[]
		### RED CHANNEL ###
		self.unmapped=0
		### BLUE CHANNEL ###
		self.mapqFunc=255.0/maxMapq
		### GREEN CHANNEL ###
		self.strandMatch=255
		self.strandDif=127.5
		#self.forward=127.5
		#self.reverse=255
		### IMAGE SCALING ###
		self.iwidth=800
		self.iheight=300
	def drawCanvas(self,leftClip,rightClip):
		for x in range(leftClip-self.flank,leftClip+self.flank): self.canvas.append(x)
		for x in range(rightClip-self.flank,rightClip+self.flank): self.canvas.append(x)
		self.canvasLeftMin,self.canvasRightMin=leftClip-self.flank, rightClip-self.flank
	def drawCanvasCI(self,leftCI,rightCI):
		for x in range(leftCI[0]-self.flank,leftCI[1]+self.flank): self.canvas.append(x)
		for x in range(rightCI[0]-self.flank,rightCI[1]+self.flank): self.canvas.append(x)
		self.canvasLeftMin,self.canvasRightMin=leftCI[0]-self.flank, rightCI[1]-self.flank
	def drawInsertionCanvas(self,leftClip,rightClip,iSize):
		for x in range(leftClip-self.flank,leftClip+1): self.canvas.append(x)
		for x in range((-1*iSize),0): self.canvas.append(x)
		for x in range(rightClip,rightClip+self.flank+1): self.canvas.append(x)
		self.canvasLeftMin, self.canvasRightMin =leftClip-self.flank, rightClip
	def transformMapq(self,mapq):
		mapq = self.mapqFunc*mapq
		if mapq > 255: mapq=255
		return mapq
	def zeroPix(self): return [[0,0,0] for x in self.canvas]
	def svPainter(self,reads):
		self.mapped=50
		self.clip=255
		for name in reads:
			Read=reads[name]
			Read.pixelPrep()
			if Read.leftClip!=None and Read.rightClip!=None: self.twoClip.append((Read.mapq+Read.score,name))
			elif (Read.leftClip!=None and Read.rightClip==None) or (Read.leftClip==None and Read.rightClip!=None): self.oneClip.append((Read.mapq+Read.score,name)) 
			else: self.mappedAln.append((Read.mapq,name))
			if Read.sameStrand==True: self.samePix.append(name) 
			else: self.diffPix.append(name)
			for Aln in Read.alignments:
				tmp=[]
				mapq=self.transformMapq(Aln.mapq)
				strandPix=self.strandMatch
				if Read.strandPix != Aln.strand: strandPix=self.strandDif
				for x in self.canvas:
					if x==Aln.leftClip or x==Aln.rightClip: tmp.append([self.clip,mapq,strandPix])
					elif x!=Aln.leftClip and x!=Aln.rightClip and x in Aln.pos: tmp.append([self.mapped,mapq,strandPix])
					elif x not in Aln.pos: tmp.append([self.unmapped,self.unmapped,self.unmapped])
				if self.pix.get(name)==None: self.pix[name]=tmp
				else: self.pix[name]=pixelUnion(tmp,self.pix[name])
	def insertionPainter(self,reads,Ins):
		self.mapped=63.75
		self.clip1=127.5
		self.clip2=191.25
		self.ins=255
		for name in reads:
			Read=reads[name]
			Read.readPosition()
			if Read.left >= self.canvasLeftMin and Read.right >= self.canvasRightMin: self.bridgingReads.append(name)
			masterStrand=Read.alignments[0].strand
			for Aln in Read.alignments:
				tmp=[]
				mapq=self.transformMapq(Aln.mapq)
				strandPix=self.forward
				if Aln.strand=='-':strandPix=self.reverse
				for x in self.canvas:
					if x<0 and Read.insertion==None: tmp.append([self.unmapped,self.unmapped,self.unmapped])
					elif x<0 and Read.insertion!=None:
						if x<Read.insertion+(-1*Ins.size):
							tmp.append([self.ins,self.unmapped,self.unmapped])
							self.insertionReads.append(name)
						else: tmp.append([self.unmapped,self.unmapped,self.unmapped]) 
					if x==Aln.rightClip or x==Aln.leftClip:
						tmp.append([self.clip1,mapq,strandPix])
						self.clippedReads.append(name)
					elif x!=Aln.leftClip and x!=Aln.rightClip and x in Aln.pos: tmp.append([self.mapped,mapq,strandPix])
					elif x not in Aln.pos and x >=0: tmp.append([self.unmapped,self.unmapped,self.unmapped])
				if masterStrand=='+':
					if self.samePix.get(name)==None: self.samePix[name]=tmp
					else: self.samePix[name]=pixelUnionInsertion(tmp,self.samePix[name],self.clip1,self.clip2)
				else:
					if self.diffPix.get(name)==None: self.diffPix[name]=tmp
					else: self.diffPix[name]=pixelUnionInsertion(tmp,self.diffPix[name],self.clip1,self.clip2)
	def orderPixelsDelDup(self,MAX):
		"""print the reads in order, for deletions and duplications"""
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.samePix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.samePix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.samePix)
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.diffPix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.diffPix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.diffPix)
		for x in self.order[0:MAX-1]: self.readPix.append(self.pix[x])
		for x in range(MAX-len(self.readPix)): self.readPix.append(self.zeroPix())
	def orderPixelsInversion(self,MAX):
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.diffPix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.diffPix)
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.samePix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.samePix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.samePix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.diffPix)
		for x in self.order[0:MAX-1]: self.readPix.append(self.pix[x])
		for x in range(MAX-len(self.readPix)): self.readPix.append(self.zeroPix())
	def orderPixelsInsertion(self):
		for x in self.samePix:
			if x in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.samePix:
			if x in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.samePix:
			if x in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.samePix:
			if x in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.samePix:
			if x not in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.samePix:
			if x not in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.samePix:
			if x not in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.samePix:
			if x not in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.samePix[x])
		for x in self.diffPix:
			if x in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.diffPix[x])
		for x in self.diffPix:
			if x in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.diffPix[x])
		for x in self.diffPix:
			if x in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.diffPix[x])
		for x in self.diffPix:
			if x in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.diffPix[x])
		for x in self.diffPix:
			if x not in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.diffPix[x])
		for x in self.diffPix:
			if x not in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.diffPix[x])
		for x in self.diffPix:
			if x not in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.diffPix[x])
		for x in self.diffPix:
			if x not in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.diffPix[x])
	def printPixels(self,SV,o):
		unscaled='{}_{}_{}_{}_{}_unscaled.png'.format(o,SV.chrom,SV.start,SV.end,SV.svtype)
		scaled='{}_{}_{}_{}_{}_scaled.png'.format(o,SV.chrom,SV.start,SV.end,SV.svtype)
		dat='{}_{}_{}_{}_{}_pixels.txt'.format(o,SV.chrom,SV.start,SV.end,SV.svtype)
		img_unscaled=smp.toimage(self.readPix)
		img_scaled=smp.toimage(smp.imresize(self.readPix,(self.iheight,self.iwidth)))
		img_unscaled.save(unscaled)
		img_scaled.save(scaled)
		ofh = open(dat,'w')
		for x in self.readPix:
			tmp=[]
			for y in x: tmp.append(','.join(map(str,y)))
			ofh.write('\t'.join(tmp)+'\n')
		ofh.close()

