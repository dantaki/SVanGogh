#!/usr/env python
from Arguments import Arguments
from Bam import Bam
from Painter import Painter
from tqdm import tqdm
def main():
	Args=Arguments()
	if Args.progress==True:
		for SV in tqdm(Args.regions): iterator(Args,SV)
	else: 
		for SV in Args.regions: iterator(Args,SV)
def iterator(Args,SV):
	bam = Bam(SV,Args)
	READS, Bosch=bam.reads,Painter(Args)
	Bosch.drawCanvas(Args.maxFlank,bam.medStart,bam.medEnd,bam.hasIns)
	Bosch.svPainter(READS,SV.svtype)
	if SV.svtype!='INV': Bosch.orderPixels(Args.maxReads)
	elif SV.svtype=='INV': Bosch.orderPixelsInversion(Args.maxReads)
	if len(Bosch.readPix)>0: Bosch.printPixels(SV,Args)
	if Args.verbose==True: print "=^..^= finished =^..^="