#include "SpecialFunctions.h"
/*pdf gamma*/
double pdfgamma(double tt,double s, double k, double theta)
{
	double salida,tau=tt-s;
	if (tt<=0){salida =0.;}
	else//{salida = exp(-tau/(k*theta))/(k*theta);}
	{salida = pow(tau,k-1)*exp(-tau/theta)/(exp(gammln(k))*pow(theta,k));}
	
	return salida;
}
/*Fbar - supervivencia*/
double Fbar(double tt,double s, double k, double theta)
{
	double salida,tau=tt-s;
	if (tt<0.)
	{
		salida=0.;
	}
	else
	{
		//salida = exp(-tau/(k*theta));
		salida = gammq(k,tau/theta);
		/*if (tau<k*theta){salida=1.;}
		else{salida=0.;}*/
	}
	return salida;
}
/*1-Fbar*/
double F(double tt,double s, double k, double theta)
{
	double salida,tau=tt-s;
	if (tt<0.)
	{
		salida = 0.;
	}
	else
	{
		salida = 1.0 - gammq(k,tau/theta);
		//salida = 1.0 - exp(-tau/(k*theta));
		/*if (tau<k*theta){salida=0.;}
		else{salida=1.;}*/
	}
	return salida;
}
/*fuerza de infeccion*/
double force(double S, double I, double N,double bbeta)
{
	return bbeta*S*I/N ;
}
/*Resuelve volterra para un kernel Fbar*/
double solve_Fbar_rk2(double y[],int n,double ddt, double k, double theta)
{
	double sum=0.,Kn,salida;
	int i=0;
	Kn=Fbar(n*ddt,(n-1)*ddt,k,theta);
	for (i=0;i<=n-2;i++)
	{
		sum+=Fbar(n*ddt,i*ddt,k,theta)*y[i]/Kn;
	}
	salida = ddt*Kn*(sum + y[n-1]);
	return salida;
	
}
/*resuelve volterra para un kernel 1-Fbar*/
double solve_F_rk2(double y[],int n,double ddt, double k, double theta)
{
	double sum=0.,Kn,salida;
	int i=0;
	Kn=F(n*ddt,(n-1)*ddt,k,theta);
	for (i=0;i<=n-2;i++)
	{
		sum+=F(n*ddt,i*ddt,k,theta)*y[i]/Kn;
	}
	salida = ddt*Kn*(sum + y[n-1]);
	return salida;
	
}
/*resuelve volterra para un kernel = 1.*/
double solve_rk2(double y[],int n,double ddt)
{
	double sum=0.,Kn,salida;
	int i=0;
	Kn=1.;
	for (i=0;i<=n-2;i++)
	{
		sum+=1.*y[i]/Kn;
	}
	salida = ddt*Kn*(sum + y[n-1]);
	return salida;
	
}
/*resuelve Volterra para dos kernel Fbar Fbar*/
double solve_Fbar_Fbar_rk2(double y[],int n, double ddt, double k1, double theta1, double k2, double theta2)
{
	double sum=0.,Kn,salida;int i=0;
	Kn=Fbar(n*ddt,(n-1)*ddt,k1,theta1)*Fbar(n*ddt,(n-1)*ddt,k2,theta2);
	for (i=0;i<=n-2;i++)
	{
		sum+=Fbar(n*ddt,i*ddt,k1,theta1)*Fbar(n*ddt,i*ddt,k2,theta2)*y[i]/Kn;
	}
	salida = ddt*Kn*(sum + y[n-1]);
	return salida;
	
}
/*resuelve volterra para los kernel F y Fbar, en ese orden entran los parametros*/
double solve_F_Fbar_rk2(double y[],int n,double ddt, double k1, double theta1, double k2, double theta2)
{
	double sum=0.,Kn,salida;int i=0;
	Kn=F(n*ddt,(n-1)*ddt,k1,theta1)*Fbar(n*ddt,(n-1)*ddt,k2,theta2);
	for (i=0;i<=n-2;i++)
	{
		sum+=F(n*ddt,i*ddt,k1,theta1)*Fbar(n*ddt,i*ddt,k2,theta2)*y[i]/Kn;
	}
	salida = ddt*Kn*(sum + y[n-1]);
	return salida;
	
}
/*resuelve volterra para un kernel pdf, a tasa constante*/
double solve_pdf_c(double y,int n,double ddt, double k, double theta)
{
	double sum=0.,Kn,salida;int i=0;
	Kn=pdfgamma(n*ddt,(n-1)*ddt,k,theta);
	for (i=0;i<=n-2;i++)
	{
		sum+=pdfgamma(n*ddt,i*ddt,k,theta)*y/Kn;
	}
	salida = ddt*Kn*(sum + y);
	return salida;
	
}
