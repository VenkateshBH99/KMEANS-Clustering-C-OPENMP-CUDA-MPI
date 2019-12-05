#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

typedef struct
{
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;
	
} Point;

//Read dimensions of the image
void readImageSize(FILE *ifp,int* a,int* b)
{
	fscanf(ifp,"%d\n",a);
	printf("%d\n",*a);

	fscanf(ifp,"%d\n",b);
	printf("%d\n",*b);
}

//reads the ifp file and stores in structure
void readPoints(FILE* ifp,Point *points, int num_points)
{
	int i;
	for(i=0;i<num_points;i++)
	{
		fscanf(ifp,"%lf,%lf,%lf,%lf,%lf", &points[i]._r, &points[i]._g, &points[i]._b, &points[i]._m, &points[i]._n);
	}
}

//Initialize random points as assumed means
void initialize(Point* mean,int K, int num_points, Point* points)
{
	int i, a, p=2;
	srand(time(NULL));
	for(i=0;i<K;i++)
	{
		a = num_points/p;
		mean[i]._r = points[a]._r;
		mean[i]._g = points[a]._g;
		mean[i]._b = points[a]._b;
		mean[i]._m = points[a]._m;	
		mean[i]._n = points[a]._n;
		p++;
	}
}

//Initially all points doesn't belong any cluster
int IntClusterMem(int *cluster, int num_points)
{
	int i;
	for(i=0;i<num_points;i++)
	{
		cluster[i]=-1;
	}	
}

//Euclidean Distance
double calculateDistance(Point point1,Point point2)
{
	return sqrt((pow((point1._r-point2._r),2)+pow((point1._g-point2._g),2)+pow((point1._b-point2._b),2)));	
}

//Assign the points to the clusters
int pointsCluster(Point point,Point* mean,int K)
{
	int parent=0;
	double dist = 0;
	double minDist=calculateDistance(point,mean[0]);
	int i;
	for(i=1;i<K;i++)
	{	
		dist=calculateDistance(point,mean[i]);
		if(minDist>=dist)
		{
			parent=i;
			minDist=dist;
		}
	}
return parent;
}

//calculate new mean
void calcNewMean(Point* points,int* cluster,Point* mean,int K,int num_points)
{
	Point* newMean=malloc(sizeof(Point)*K);
	int* members=malloc(sizeof(int)*K);
	int i;
	for(i=0;i<K;i++)
	{
		members[i]=0;
		newMean[i]._r=0;
		newMean[i]._g=0;
		newMean[i]._b=0;
		newMean[i]._m=0;
		newMean[i]._n=0;
	}	
	for(i=0;i<num_points;i++)
	{
		members[cluster[i]]++;
		newMean[cluster[i]]._r+=points[i]._r;
		newMean[cluster[i]]._g+=points[i]._g;
		newMean[cluster[i]]._b+=points[i]._b;
	}
	for(i=0;i<K;i++)	
	{
		if(members[i]!=0.0)
		{
			newMean[i]._r/=members[i];
			newMean[i]._g/=members[i];
			newMean[i]._b/=members[i];
		}
		else
		{
			newMean[i]._r=0;
			newMean[i]._g=0;
			newMean[i]._b=0;
			newMean[i]._m=0;
			newMean[i]._n=0;
		}
	}
	for(i=0;i<K;i++)
	{
		mean[i]._r=newMean[i]._r;
		mean[i]._g=newMean[i]._g;
		mean[i]._b=newMean[i]._b;
		mean[i]._m=newMean[i]._m;
		mean[i]._n=newMean[i]._n;
	}	
}

//Termination Condition
int chkConvrg(int *before_clusters,int *after_cluster,int num_points, float tol)
{
	int i;
	tol = num_points*tol;
	for(i=0;i<num_points;i++)
		if((before_clusters[i]-after_cluster[i])>tol)
			return -1;
	return 0;
}

int main(int argc, char* argv[])
{
	int rank;
	int size;
	struct timespec start_t, stop_t;
	clock_gettime(CLOCK_MONOTONIC,&start_t);
	double span_t;
	int K=0;
	int num_points;
	int i,j,l=1;
	int job_size;
	int job_done=0;
	int x,y;
	float tol;
	double TIC,TOC,tic,ticIn, tocIn, ticC1, tocC1, ticC2, tocC2, ticC3, tocC3, ticNewMean, tocNewMean, ticChCon, tocChCon,ticOut,tocOut,tocPC4,ticPC4,ticPCal,tocPCal,ticPC5,tocPC5, ticPC6, tocPC6;
	double tspanNewMean = 0.0, tspanC2=0.0, tspanChCon=0.0, tspanC3=0.0, tspanOut=0.0,tspanPCal=0.0,tspanPC4=0.0,tspanPC5=0.0, tspanPC6=0.0;
	float check;

	Point* mean;
	Point* points;
	Point* get_points;
	int * formed_clusters;
	int * before_clusters;
	int * after_cluster;

	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	clock_gettime(CLOCK_MONOTONIC,&stop_t);
	span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
	TIC = span_t;

	//creating of derived MPI structure
	MPI_Datatype MPI_POINT;
	MPI_Datatype type=MPI_DOUBLE;
	int blocklen=2;
	MPI_Aint disp=0; 
	//C type that holds any valid address
	MPI_Type_create_struct(1,&blocklen,&disp,&type,&MPI_POINT);
	MPI_Type_commit(&MPI_POINT);

	if(rank!=0)
	{
		//Receiving the cluster
		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticPC4 = span_t;
		MPI_Recv(&job_size ,1 ,MPI_INT ,0,0,MPI_COMM_WORLD,&status);
		MPI_Recv(&K,1 ,MPI_INT ,0,0,MPI_COMM_WORLD,&status);
		mean =malloc(sizeof(Point)*K);
		MPI_Recv(mean ,K,MPI_POINT,0,0,MPI_COMM_WORLD,&status);
		points =(Point*)malloc(sizeof(Point)*job_size);
		after_cluster =(int*)malloc(sizeof(int)*job_size);
		
		for(i=0;i<job_size;i++)
		{
			
			MPI_Recv(&points[i]._r,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._g,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._b,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._m,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._n,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));

		tocPC4 = span_t;
		tspanPC4 = (double)(tocPC4-ticPC4);
			
		while(1)
		{
			ticPCal = MPI_Wtime();
			for(i=0;i<job_size;i++)
			{
				after_cluster[i]=pointsCluster(points[i],mean,K);
			}
			tocPCal = MPI_Wtime();
			if(rank == 1)
				tspanPCal += (double)(tocPCal-ticPCal);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticPC5 = span_t;
			
	        MPI_Send(after_cluster,job_size, MPI_INT,0, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocPC5 = span_t;
			tspanPC5 += (double)(tocPC5-ticPC5);

			MPI_Bcast(&job_done,1, MPI_INT,0,MPI_COMM_WORLD);

			if(job_done==1) //No more work to be done
			break;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticPC6 = span_t;

			//Receiving recently created mean from master
			//MPI_Recv(mean,K,MPI_POINT,0,0, MPI_COMM_WORLD,&status);
			MPI_Bcast(mean,K, MPI_POINT,0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocPC6 = span_t;
			tspanPC6 += (double)(tocPC6-ticPC6);
		}

	}
	else
	{
		
		
		FILE *ifp;
		ifp=fopen(argv[1],"r");
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
	    span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticIn = span_t;
		readImageSize(ifp,&x,&y);
		K = atoi(argv[3]);
		num_points = x*y;
		points =(Point*)malloc(sizeof(Point)*num_points);
		readPoints(ifp,points,num_points);
		fclose(ifp);
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
	    span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		tocIn = span_t;
		//Allocates memory
		before_clusters =(int*)malloc(sizeof(int)*num_points);
		after_cluster=(int*)malloc(sizeof(int)*num_points);
		mean =malloc(sizeof(Point)*K);
		//cuts job depending on number of threads

		check = num_points%(size-1);
		if(check==0.00)
		{
			job_size=num_points/(size-1);
		}
		else
		{
			printf("\n Enter no. of Processes as n+1 (where n divides  %d  in equal parts)\n\n",num_points);
			exit(1);
		}
		//initializing to default values
		initialize(mean,K,num_points,points);
		IntClusterMem(before_clusters,num_points);
		IntClusterMem(after_cluster,num_points);
		tol = 0;
		printf("Tolerance = %f\n",tol);

		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticC1=span_t;
	
		//Sending the essential cluster to other processors
		for(i=1;i<size;i++)
		{
			MPI_Send(&job_size ,1 , MPI_INT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&K ,1 , MPI_INT ,i,0,MPI_COMM_WORLD);
			MPI_Send(mean ,K, MPI_POINT ,i,0,MPI_COMM_WORLD);
			for(j=0;j<job_size;j++)
			{
			MPI_Send(&points[j]._r+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._g+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._b+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._m+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._n+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
	        span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));

		tocC1 = span_t;
		//master processor job
		while(1)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticC2 = span_t;
			ticC2 = span_t;
			for(i=1;i<size;i++)
				MPI_Recv(after_cluster+(job_size*(i-1)),job_size,MPI_INT,i,0,MPI_COMM_WORLD,&status);
	        MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocC2 = span_t;

			tspanC2 += (double)(tocC2-ticC2);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
	                span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticNewMean = span_t;		
			calcNewMean(points,after_cluster,mean,K,num_points);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
	                span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocNewMean = span_t;
			tspanNewMean += (double)(tocNewMean-ticNewMean);		
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
	                span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticChCon = span_t;

			if(chkConvrg(after_cluster,before_clusters,num_points,tol)==0)
			{
				printf("K-mean algorithm Converged at iteration %d\n",l);
				job_done=1;

			}
			else
			{
				l++;
				for(i=0;i<num_points;i++)
					before_clusters[i]=after_cluster[i];
			}
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
	                span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocChCon = span_t;
			tspanChCon += (double)(tocChCon-ticChCon);

			//Informing slaves that no more job to be done
			MPI_Bcast(&job_done,1, MPI_INT,0,MPI_COMM_WORLD);

			if(job_done==1)
				break;
			 MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tic = span_t;
			
			//Sending the recently created mean
			MPI_Bcast(mean,K, MPI_POINT,0, MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocC3 = span_t;
			tspanC3 += (double)(tocC3-tic);

		}
		
		//Outputting to the ofp file
		FILE* ofp=fopen(argv[2],"w");
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
	                span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticOut = span_t;
		
		fprintf(ofp,"%d\n",x);
		fprintf(ofp,"%d\n",y);
		for(i=0;i<K;i++)
		fprintf(ofp,"%d,%d,%d,%d,%d\n",(int)mean[i]._r,(int)mean[i]._g,(int)mean[i]._b,(int)mean[i]._m,(int)mean[i]._n);
		for(i=0;i<num_points;i++)
		fprintf(ofp,"%d,%d,%d,%d,%d,%d\n",(int)points[i]._r,(int)points[i]._g,(int)points[i]._b,(int)points[i]._m,(int)points[i]._n,(int)after_cluster[i]+1);
		fclose(ofp);
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
	    span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		tocOut = span_t;
		tspanOut = (double)(tocOut-ticOut);
	}
	clock_gettime(CLOCK_MONOTONIC,&stop_t);
	span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
	 TOC = span_t;


	if(rank==0)
	{
		printf("Input time : %f sec\n", tocIn - ticIn);
		printf("send initial data time : %f sec\n", tocC1 - ticC1);

		printf("Recieve results from slave time : %f sec\n", tspanC2);
		printf("NewMean calculation time : %f sec\n", tspanNewMean);
		printf("Convergence check time : %f sec\n", tspanChCon);
		printf("send newmeans to slaves time : %f sec\n", tspanC3);
		printf("Total time for iterative communication data (mean, after cluster) : %f\n",tspanC2+tspanC3);
		printf("Output time : %f sec\n", tocOut - ticOut);
		printf("\tMaster time (Communication) : %f sec\n", ((tocC1 - ticC1) +  tspanC2 + tspanC3));
		printf("\tMaster time (Calculation) : %f sec\n", (tspanNewMean+  tspanChCon));
		printf("\tMaster time (File operation) : %f sec\n", ((tocIn - ticIn)+  (tocOut - ticOut)));
		printf("\tTotal Master time : %f sec\n", ((tocIn - ticIn)+  (tocOut - ticOut)) + (tspanNewMean+  tspanChCon) + ((tocC1 - ticC1) +  tspanC2 + tspanC3));
	    printf("\nTotal Time : %f sec\n",TOC - TIC);
	    FILE* timef=fopen("TIME/MPI_Time","a");
		fprintf(timef,"%f\n",TOC-TIC);
		fclose(timef);
	}
	else if(rank==1)
	{
	
		printf("Calculate After_Cluster time : %f sec\n", tspanPCal);
	}
	//End of all

	MPI_Finalize();

	     return 0;
}




