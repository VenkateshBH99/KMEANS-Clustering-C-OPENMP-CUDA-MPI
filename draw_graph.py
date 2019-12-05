import matplotlib.pyplot as plt 


file1=open("TIME/C_Time","r")
C=file1.readlines()
for i in range(len(C)):
	C[i]=float(C[i].rstrip())
file1.close()

file1=open("TIME/OMP_Time","r")
OMP=file1.readlines()
for i in range(len(OMP)):
	OMP[i]=float(OMP[i].rstrip())
file1.close()

file1=open("TIME/MPI_Time","r")
MPI=file1.readlines()
for i in range(len(MPI)):
	MPI[i]=float(MPI[i].rstrip())
file1.close()
# line 1 points 
x1 = [i+1 for i in range(len(C)) ]  
# plotting the line 1 points 
plt.plot(x1, C, label = "Simple C")  
plt.plot(x1, OMP, label = "Openmp") 
plt.plot(x1, MPI, label = "MPI") 
plt.ylim(0,1) 
plt.xlim(1,8) 
# naming the x axis 
plt.xlabel('Image_no') 
# naming the y axis 
plt.ylabel('Time in sec') 
# giving a title to my graph 
plt.title('K-means') 

# show a legend on the plot 
plt.legend() 

# function to show the plot 
plt.show() 
