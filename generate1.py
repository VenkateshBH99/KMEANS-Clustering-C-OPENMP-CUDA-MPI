import numpy as np
from PIL import Image
import random
method=input()
print(method)
k=5
for img in range(8):
	file=open(method+"/"+str(img)+"out.txt","r")
	colors=[[] for i in range(k)]
	for i in range(k):
		for j in range(3):
			colors[i].append(random.randint(0,255))



	lines=file.readlines()
	brr=[]
	y=int(lines[1].rstrip())
	for i in range(k+2,len(lines)):
		lines[i]=list(map(int,lines[i].rstrip().split(',')))
		brr.append(colors[lines[i][5]-1]+[lines[i][4]-1,lines[i][3]-1,255])


	#brr.sort(key=lambda k:k[3])
	array=[[] for i in range(y)]
	for i in range(len(brr)):
		array[brr[i][3]].append(brr[i])
	#for i in array:
	#	i.sort(key=lambda k:k[4])

	img_array=[]
	for i in range(len(array)):
		for j in range(len(array[0])):
			
			array[i][j].pop(3)
			array[i][j].pop(3)
	c_array = np.asarray(array)
	print(c_array)
	invimg = Image.fromarray(c_array.astype('uint8'))
	invimg.save(method+'/'+str(img)+'.png')
	file.close()

