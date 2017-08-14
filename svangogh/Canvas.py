from math import ceil,floor
class Canvas():
	def __init__(self,flank=None,startClip=None,endClip=None,svtype=None):
		self.coord=[]
		self.size=flank*4
		if svtype=='INS': flank=int((self.size/3.)/2.)
		leftMin, leftMax=int(startClip-flank),int(startClip+flank)
		rightMin, rightMax=int(endClip-flank),int(endClip+flank)
		if leftMax > rightMin:
			dif = endClip-startClip
			leftFlank, rightFlank= ceil(dif/2.),floor(dif/2.)
			leftMin, leftMax = int(startClip-(flank+(flank-leftFlank))),int(startClip+leftFlank)
			rightMin,rightMax = int(endClip-rightFlank), int(endClip+(flank+(flank-rightFlank)))
			if rightMin==leftMax: rightMin+=1
		for x in range(leftMin,leftMax): self.coord.append(x)
		if svtype=='INS':
			for x in reversed(range((-1*(self.size-((rightMax-rightMin)+(leftMax-leftMin)))),0)): self.coord.append(x)
		for x in range(rightMin,rightMax): self.coord.append(x)