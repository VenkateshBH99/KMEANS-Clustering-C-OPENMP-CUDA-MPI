import numpy as np
from PIL import Image

file=open("img_file.txt","r")
img_f=file.readlines()
for i in range(len(img_f)):
	img_f[i]=img_f[i].strip()
#img = Image.open('1.jpeg')
for i in range(len(img_f)):
	img=Image.open("Images/"+img_f[i])
	array = np.array(img)
	print(array.shape)


	file=open("Image_Detail/"+str(i)+".txt","w")
	print(i)
	file.write(str(len(array[0]))+"\n")
	file.write(str(len(array))+"\n")


	for i in range(len(array)):
		for j in range(len(array[0])):
			file.write(str(array[i][j][0])+","+str(array[i][j][1])+","+str(array[i][j][2])+","+str(j+1)+","+str(i+1)+"\n")

	file.close()
