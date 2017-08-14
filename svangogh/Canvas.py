from math import ceil,floor
class Canvas():
	def __init__(self,flank=None,startClip=None,endClip=None):
		self.coord=[]
		leftMin, leftMax=startClip-flank,startClip+flank
		rightMin, rightMax=endClip-flank, endClip+flank
		if leftMax > rightMin:
			dif = endClip-startClip
			leftFlank, rightFlank= ceil(dif/2.),floor(dif/2.)
			leftMin, leftMax = int(startClip-(flank+(flank-leftFlank))),int(startClip+leftFlank)
			rightMin,rightMax = int(endClip-rightFlank), int(endClip+(flank+(flank-rightFlank)))
		for x in range(leftMin,leftMax): self.coord.append(x)
		for x in range(rightMin,rightMax): self.coord.append(x)
	def drawInsertionCanvas(self,startClip,endClip,iSize):
		for x in range(startClip-self.flank,startClip+1): self.canvas.append(x)
		for x in range((-1*iSize),0): self.canvas.append(x)
		for x in range(endClip,endClip+self.flank+1): self.canvas.append(x)
		self.leftMin, self.rightMin =startClip-self.flank, endClip
