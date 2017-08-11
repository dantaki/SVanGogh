#!/usr/env python
from Bam import Bam
from Cigar import Cigar
from Alignment import Alignment
from Read import Read
from Insertion import Insertion
from Painter import Painter
from front_end import Arguments
from tqdm import tqdm
def main():
	Args=Arguments()
	if Args.verbose==True:
		for SV in tqdm(Args.regions): iterator(Args,SV)
	else: 
		for SV in Args.regions: iterator(Args,SV)
def iterator(Args,SV):
	##### ITERATE READS #####
	bam = Bam(Args.ifh)
	bam.leftFlank(SV,Args)
	if SV.svtype=='DEL' or SV.svtype =='DUP' or SV.svtype == 'INV': bam.rightFlank(SV,Args)
	#########################
	READS = bam.reads
	bam.pixelPrep(SV,Args)
	Bosch=Painter(Args.maxFlank,Args.maxMapq)
	if SV.svtype!='INS':
		#Bosch.drawCanvas(SV.leftCI,SV.rightCI)
		if SV.svtype=='DEL' or SV.svtype=='INV':
			if bam.medRightClip > bam.medLeftClip: Bosch.drawCanvas(bam.medLeftClip,bam.medRightClip)
			else: Bosch.drawCanvas(bam.medRightClip,bam.medLeftClip)
		elif SV.svtype=='DUP': 
			Bosch.drawCanvas(bam.medLeftClip,bam.medRightClip)
		Bosch.svPainter(READS)
		if SV.svtype=='DEL' or SV.svtype=='DUP': 
			Bosch.orderPixelsDelDup(Args.maxReads,Args.verbose)
		elif SV.svtype=='INV':
			Bosch.orderPixelsInversion(Args.maxReads,Argx.verbose)
		#Bosch.orderPixels()
	else: 
		Ins = Insertion()
		Ins.findInsertions(READS)
		if Ins.size==None: Ins.size=0
		if Ins.leftClip==None: Ins.leftClip=SV.start
		if Ins.rightClip==None: Ins.rightClip=SV.end
		Bosch.drawInsertionCanvas(Ins.leftClip,Ins.rightClip,Ins.size)
		Bosch.insertionPainter(READS,Ins)
		Bosch.orderPixelsInsertion()
	if len(Bosch.readPix)>0: Bosch.printPixels(SV,Args.ofh)
