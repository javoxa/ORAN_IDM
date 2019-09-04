/*generador gamma*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "C:\libreria\aleatorios.h"
#include "/Users/javo/libreria/aleatorios.h"
/*
 * Return a deviate distributed as a gamma distribution of integer ia,
 * i.e., a waiting time to the iath event poisson process of unit a 
 * mean, using ran1(idum) as the source of uniform deviates*/
float gamdev(int ia, int *idum)
{
	float ran1(long *idum);
	void nrerror(char error_tex []);
	int j;
	float am,e,s,v1,v2,x,y;
	
	//if (ia < 1)nerror("Error fatal en gamdev \n");
	if (ia < 6) 
	{
		x=1.0;
		for(j=1;j<=ia;j++) x*=ran1(idum);
		x=-log(x);
	}
	else
	{
		do
		{
			do
			{
				do
				{
					v1=ran1(idum);
					v2=2.0*ran1(idum)-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (ran1(idum) > e);
	}
	return x;
}

int main()
{
	float gam;
	long int seed;
	fprintf(stderr,"seed=? \n");
	scanf("%lu",&seed);
	seed=-fabs(seed);
	int i;
	while (i<50)
	{
		gam=gamdev(&seed);
		printf("%f \n",gam);
		i++;
	}
}

