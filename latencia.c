#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "/Users/javo/libreria/aleatorios_mio.h"
//#include "/Users/javo/libreria/SpecialFunctions.h"
#include "/Users/javo/libreria/integradores_volterra.h"

/* para la mierda de latencia*/
int main(int argc,char *argv[])
{
	int days;
    float t,dt;
    float theta, k, media, var;
    float failure_rate,comulative;
    
    //parameter gamma
    /*
    var = 4.;     // var    =  k * theta * theta;
    media = 6.;  // media  =  k * theta;
    
    theta = var/media;
    k = media*media/var;*/
    
    k = 3.;
    media = 3.;
    
    var = media*media/k;
    theta = media/k;
    
    /* entrada */
    t = 0.;
    dt = 0.1;
    
    while( t < 30. ){
    
    failure_rate = pdfgamma( t, 0., k, theta) / Fbar( t, 0., k, theta);
    
    comulative = F( t, 0., k, theta);
    
    printf("%.2f\t%.2f\t%.2f\n",t, comulative, failure_rate);
    
    t = t + dt;
    
    }

    
    
	
	
	exit(0);
}
