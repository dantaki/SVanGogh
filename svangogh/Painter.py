#!/usr/env python
import scipy.misc as smp
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
class Painter():
	def __init__(self,maxFlank):
		self.canvas=[]
		self.flank=maxFlank
		self.canvasLeftMin=None
		self.canvasRightMin=None
		self.forwardPix={} # [ReadName]=[ readPix ]
		self.reversePix={}
		self.bridgingReads=[]
		self.clippedReads=[]
		self.insertionReads=[]
		self.readPix=[]
		### RED CHANNEL ###
		self.unmapped=0
		### BLUE CHANNEL ###
		self.mapqFunc=4.25
		### GREEN CHANNEL ###
		self.forward=127.5
		self.reverse=255
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
	def svPainter(self,reads):
		self.mapped=50
		self.clip=255
		for name in reads:
			Read=reads[name]
			Read.readPosition()
			if Read.left >= self.canvasLeftMin and Read.right >= self.canvasRightMin:self.bridgingReads.append(name)
			masterStrand=Read.alignments[0].strand
			for Aln in Read.alignments:
				tmp=[]
				mapq=self.transformMapq(Aln.mapq)
				strandPix=self.forward
				if Aln.strand=='-':strandPix=self.reverse  
				for x in self.canvas:
					if x==Aln.leftClip or x==Aln.rightClip:
						tmp.append([self.clip,mapq,strandPix])
						self.clippedReads.append(name)
					elif x!=Aln.leftClip and x!=Aln.rightClip and x in Aln.pos:
						tmp.append([self.mapped,mapq,strandPix])
					elif x not in Aln.pos:
						tmp.append([self.unmapped,self.unmapped,self.unmapped])
				if masterStrand=='+': 
					if self.forwardPix.get(name)==None: self.forwardPix[name]=tmp			
					else: self.forwardPix[name]=pixelUnion(tmp,self.forwardPix[name])
				else:
					if self.reversePix.get(name)==None: self.reversePix[name]=tmp
					else: self.reversePix[name]=pixelUnion(tmp,self.reversePix[name])
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
					if self.forwardPix.get(name)==None: self.forwardPix[name]=tmp
					else: self.forwardPix[name]=pixelUnionInsertion(tmp,self.forwardPix[name],self.clip1,self.clip2)
				else:
					if self.reversePix.get(name)==None: self.reversePix[name]=tmp
					else: self.reversePix[name]=pixelUnionInsertion(tmp,self.reversePix[name],self.clip1,self.clip2)
	def orderPixels(self):
		for x in self.forwardPix:
			if x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.reversePix:
			if x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.reversePix[x])
		for x in self.reversePix:
			if x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.reversePix[x])
		for x in self.reversePix:
			if x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.reversePix[x])
		for x in self.reversePix:
			if x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.reversePix[x])
	def orderPixelsInsertion(self):
		for x in self.forwardPix:
			if x in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x not in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x not in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x not in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.forwardPix:
			if x not in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.forwardPix[x])
		for x in self.reversePix:
                        if x in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.reversePix[x])
                for x in self.reversePix:
                        if x in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.reversePix[x])
                for x in self.reversePix:
                        if x in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.reversePix[x])
                for x in self.reversePix:
                        if x in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.reversePix[x])
                for x in self.reversePix:
                        if x not in self.insertionReads and x in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.reversePix[x])
                for x in self.reversePix:
                        if x not in self.insertionReads and x in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.reversePix[x])
                for x in self.reversePix:
                        if x not in self.insertionReads and x not in self.clippedReads and x in self.bridgingReads: self.readPix.append(self.reversePix[x])
                for x in self.reversePix:
                        if x not in self.insertionReads and x not in self.clippedReads and x not in self.bridgingReads: self.readPix.append(self.reversePix[x])
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

