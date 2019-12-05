#!/bin/sh
> TIME/C_Time
> TIME/OMP_Time
> TIME/MPI_Time
echo "converting images to matrices......"
python3 img.py > Image_details
echo "conversion completes........."
echo "compiling codes......"
gcc -o kmeanc kmean.c -lm
clang++ -fopenmp -o kmeanomp kmean_omp.c -lm
mpicc -o kmeanmpi kmean_mpi.c  -lm
echo "compilation completes......"
echo "executing codes on different images........"
for((i=0;i<8;i++))
do
	img=Image_Detail/"$i".txt
	mpiout=MPI/"$i"out.txt
	cout=C/"$i"out.txt
	ompout=OMP/"$i"out.txt
	./kmeanc "$img" "$cout" 5 > Detail/C/"$i"infoc
	./kmeanomp "$img" "$ompout" 512 5 > Detail/OMP/"$i"infoomp  
	mpirun -np 2 ./kmeanmpi "$img" "$mpiout" 5 > Detail/MPI/"$i"infompi
done
echo "execution completes....."
echo "generating images from the output......"
for k in MPI C OMP
do
	echo "$k" | python3 generate1.py > G
done
echo "generation completes...."
echo "Plotting graphs...."
python3 draw_graph.py


