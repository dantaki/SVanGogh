from Canvas import Canvas
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
def appendOrder(i,o,t):
	for x in sorted(i, key=itemgetter(0)):
		if x[1] in list(set(t)-set(o)): o.append(x[1])
class Painter():
	def __init__(self,Args=None):
		self.canvas=[]
		self.samePix=[] # [ReadName]=[ readPix ]
		self.diffPix=[]
		self.twoClip=[]
		self.oneClip=[]
		self.mappedAln=[]
		self.insAln=[]
		self.pix={}
		self.readPix=[]
		self.order=[]
		### RED CHANNEL ###
		self.unmapped=0
		### BLUE CHANNEL ###
		self.mapqFunc=255.0/Args.maxMapq
		### GREEN CHANNEL ###
		self.strandMatch=255
		self.strandDif=127.5
		### IMAGE SCALING ###
		wscale, hscale  = Args.scaling, Args.scaling
		if Args.wscaling != None: wscale=Args.wscaling
		if Args.hscaling != None: hscale=Args.hscaling
		self.iwidth, self.iheight =Args.maxFlank*wscale*4, Args.maxReads*hscale
	def drawCanvas(self,flank,startClip,endClip,hasIns): self.canvas=Canvas(flank,startClip,endClip,hasIns).coord
	def transformMapq(self,mapq):
		mapq = self.mapqFunc*mapq
		if mapq > 255: mapq=255
		return mapq
	def zeroPix(self): return [[0,0,0] for x in self.canvas]
	def svPainter(self,reads=None):
		self.mapped=127.5
		self.clip=255
		for name in reads:
			Read=reads[name]
			if Read.insertion!=None: self.insAln.append(name)			
			if Read.startClip!=None and Read.endClip!=None: self.twoClip.append((Read.mapq+Read.score,name))
			elif (Read.startClip!=None and Read.endClip==None) or (Read.startClip==None and Read.endClip!=None): self.oneClip.append((Read.mapq+Read.score,name)) 
			else: self.mappedAln.append((Read.mapq,name))
			if Read.sameStrand==True: self.samePix.append(name) 
			else: self.diffPix.append(name)
			for Aln in Read.alignments:
				tmp=[]
				mapq=self.transformMapq(Aln.mapq)
				strandPix=self.strandMatch
				if Read.strandPix != Aln.strand: strandPix=self.strandDif
				for x in self.canvas:
					if x<0:
						if Read.insertion==None: tmp.append([self.unmapped,self.unmapped,self.unmapped])
						elif Read.insertion!=None:
							if -1*x<=Read.insertion: tmp.append([self.clip,self.unmapped,self.unmapped])
							else: tmp.append([self.unmapped,self.unmapped,self.unmapped])
					else:
						if x==Aln.startClip or x==Aln.endClip: tmp.append([self.clip,mapq,strandPix])
						elif x!=Aln.startClip and x!=Aln.endClip and x in Aln.pos: tmp.append([self.mapped,mapq,strandPix])
						elif x not in Aln.pos: tmp.append([self.unmapped,self.unmapped,self.unmapped])
				if self.pix.get(name)==None: self.pix[name]=tmp
				else: self.pix[name]=pixelUnion(tmp,self.pix[name])
	def orderPixels(self,MAX):
		if len(self.insAln)>0: 
			appendOrder(self.twoClip,self.order,list(set(self.insAln)&set(self.samePix)))
			appendOrder(self.oneClip,self.order,list(set(self.insAln)&set(self.samePix)))
			appendOrder(self.mappedAln,self.order,list(set(self.insAln)&set(self.samePix)))
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.samePix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.samePix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.samePix)
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.diffPix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.diffPix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.diffPix)
		for x in self.order[0:MAX]: self.readPix.append(self.pix[x])
		for x in range(MAX-len(self.readPix)): self.readPix.append(self.zeroPix())
	def orderPixelsInversion(self,MAX,):
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.diffPix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.diffPix)
		if len(self.twoClip)>0: appendOrder(self.twoClip,self.order,self.samePix)
		if len(self.oneClip)>0: appendOrder(self.oneClip,self.order,self.samePix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.samePix)
		if len(self.mappedAln)>0: appendOrder(self.mappedAln,self.order,self.diffPix)
		for x in self.order[0:MAX]: self.readPix.append(self.pix[x])
		for x in range(MAX-len(self.readPix)): self.readPix.append(self.zeroPix())
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

