from math import ceil,floor
class Canvas():
	def __init__(self,flank=None,startClip=None,endClip=None,medStart=None,medEnd=None,hasIns=None):
		self.coord=[]
		self.size=flank*4
		if hasIns==True: flank=int((self.size/3.)/2.)
		leftMin,leftMax = self.medianCoord(medStart,flank)
		rightMin,rightMax = self.medianCoord(medEnd,flank)
		if startClip[1]-startClip[0] <= flank*2: 
			leftMin = startClip[0]-25
			leftMax = (flank*2)+leftMin
		if endClip[1]-endClip[0] <= flank*2: 
			rightMax = endClip[1]+25
			rightMin = rightMax-(flank*2)
		if leftMax > rightMin:
			dif = endClip[0]-startClip[1]
			leftFlank, rightFlank= ceil(dif/2.),floor(dif/2.)
			leftMax, rightMin = int(startClip[1]+leftFlank), int(endClip[0]-rightFlank)
			leftMin, rightMax = int(leftMax-(flank*2)),int((flank*2)+rightMin)
			#leftMin, leftMax = int(startClip-(flank+(flank-leftFlank))),int(startClip+leftFlank)
			#rightMin,rightMax = int(endClip-rightFlank), int(endClip+(flank+(flank-rightFlank)))
			if rightMin==leftMax: rightMin+=1
		for x in range(leftMin,leftMax): self.coord.append(x)
		if hasIns==True:
			for x in reversed(range((-1*(self.size-((rightMax-rightMin)+(leftMax-leftMin)))),0)): self.coord.append(x)
		for x in range(rightMin,rightMax): self.coord.append(x)	
	def medianCoord(self,clip=None,flank=None): return int(clip-flank),int(clip+flank)
