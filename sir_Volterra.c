#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/Users/javo/libreria/integradores_volterra.h"
//#include "parameters.h"
/*=====revange======*/
int main(void)
{
	/*declaro un puntero y abro un archivo para hacer el plot*/
	FILE *plot1,*plot2;
	char file1[20],file2[20];
	
	/*declaro las tasas*/
	//double beta,mu,gama;
	double theta_gama,k_gama,theta_mu;
	
	/*declaro las varianzas*/
	//double var_g,var_mu;
	/*declaracion de variables*/
	double t;//,tf=50.,dt=0.01;
	long int count,size =(tf/dt);
	/*valores de las tasas medias*/
	/*
	gama=1.;//periodo medio infeccioso [1/T]
	mu=1./10.;//periodo medio de vida 75 years => 1/75 = 0.013
	beta=2.2;//suc
	*/
	/*en base a las tasas, calculo los k y theta para las gamma*/
	/*para el periodo infeccioso gama*/
	double k_mu=kmu;
	k_gama = 1. /*exponencial k_g = 1*/; 
	theta_gama = 1./(k_gama*gama);
	/*================================================================*/
	/*Comenzamos variando k_mu desde 1 a 100 (en forma entera)*/
	/*para k_mu=1 caso exponencial*/
	//k_mu=1;
	theta_mu = 1./(k_mu*mu);

	//escalado(100,100,Imp,Imag);
	//wpgmt(100,100,250,1,Imp);
	//system("echo '	load \"grafico.plt\"	'|gnuplot -persist");

	return 0;
}
