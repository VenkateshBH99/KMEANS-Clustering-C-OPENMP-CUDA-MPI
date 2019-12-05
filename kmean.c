#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
typedef struct
{
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;
	
} Point;

//Reads the image dimensions a x b pixels
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

//All points having no clusters
void IntClusterMem(int *cluster, int num_points)
{
	int i;
	for(i=0;i<num_points;i++)
		{
		cluster[i]=-1;
		}	
}

// Euclidean Distance
double calculateDistance(Point point1,Point point2)
{
	return sqrt((pow((point1._r-point2._r),2)+pow((point1._g-point2._g),2)));	
}

//to calculate which cluster is the point belonging to.
int pointsCluster(Point point,Point* mean,int K)
{
	int parent=0;
	double dist=0;
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

//check for convergence
// it checks that is each points cluster remaining the same
int chkConvrg(int *before_clusters,int *after_cluster,int num_points, float tol)
{
	int i;
	tol = num_points*tol;
	for(i=0;i<num_points;i++)
		if(abs(before_clusters[i]-after_cluster[i])>tol)	
			return -1;
	return 0;
}

int main(int argc, char* argv[])
{
	int K;
	int num_points;
	int i;
	int job_done=0;
	int x,y;
	clock_t tic, toc;
	double tspan = 0.0, tspantemp=0.0;
	int iter=0;
	Point* mean;
	Point* points;
	Point* get_points;
	int * formed_clusters;
	int * before_clusters;
	int * after_cluster;
	float tol=0.0;

	K = atoi(argv[3]);

	printf("Tolerance = %.10f\n",tol);       
	//Readinf file
	FILE *ifp;
	ifp=fopen(argv[1],"r");
	readImageSize(ifp,&x,&y);
	num_points = x*y;
	points=(Point*)malloc(sizeof(Point)*num_points);
	readPoints(ifp,points,num_points);
	fclose(ifp);

	before_clusters=(int*)malloc(sizeof(int)*num_points);
	after_cluster=(int*)malloc(sizeof(int)*num_points);
	mean=malloc(sizeof(Point)*K);

	//initializing to default values
	initialize(mean,K,num_points,points);
	IntClusterMem(before_clusters,num_points);
	IntClusterMem(after_cluster,num_points);

	while(1)
	{	

		tic = clock();
		iter++;
		for(i=0;i<num_points;i++)
		{
			after_cluster[i]=pointsCluster(points[i],mean,K);
		}
		calcNewMean(points,after_cluster,mean,K,num_points);
		toc = clock();	
		tspantemp = (double)(toc-tic) / (double)CLOCKS_PER_SEC;
		tspan += tspantemp;
		if(chkConvrg(after_cluster,before_clusters,num_points,tol)==0)
		{
			printf("K-mean algorithm Converged!\n");
			job_done=1;
		}
		else
		{
			for(i=0;i<num_points;i++)
				before_clusters[i]=after_cluster[i];
		}
		
		if(job_done==1)
			break;

	}
	printf("Total Iterations = %d\n",iter);
	printf("Total time elapsed in forming clusters : %f sec\n", tspan);
	//Outputting to the ofp file

	FILE* ofp=fopen(argv[2],"w");
	fprintf(ofp,"%d\n",x);
	fprintf(ofp,"%d\n",y);
	for(i=0;i<K;i++)
	fprintf(ofp,"%d,%d,%d,%d,%d\n",(int)mean[i]._r,(int)mean[i]._g,(int)mean[i]._b,(int)mean[i]._m,(int)mean[i]._n);
	for(i=0;i<num_points;i++)
	fprintf(ofp,"%d,%d,%d,%d,%d,%d\n",(int)points[i]._r,(int)points[i]._g,(int)points[i]._b,(int)points[i]._m,(int)points[i]._n,after_cluster[i]+1);
	fclose(ofp);

	FILE* timef=fopen("TIME/C_Time","a");
	fprintf(timef,"%f\n",tspan);
	fclose(timef);

	//End of all
    return 0;
}
