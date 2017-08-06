#!/usr/env python
from Bam import Bam
from Cigar import Cigar
from Alignment import Alignment
from Read import Read
from Insertion import Insertion
from Painter import Painter
from front_end import Arguments
import numpy as np
def main():
	Args=Arguments()
	for SV in Args.regions:
		##### ITERATE READS #####
		bam = Bam(Args.ifh)
		bam.leftFlank(SV,Args)
		if SV.svtype=='DEL' or SV.svtype =='DUP' or SV.svtype == 'INV': bam.rightFlank(SV,Args)
		#########################
		READS = bam.reads
		print np.median([x[0] for x in bam.clips if x[0] !=None]) 
		Bosch=Painter(Args.maxFlank)
		if SV.svtype!='INS':
			Bosch.drawCanvas(SV.leftCI,SV.rightCI)
			Bosch.svPainter(READS)
			Bosch.orderPixels()
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
