#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "/Users/javo/libreria/aleatorios_mio.h"
#include "/Users/javo/libreria/SpecialFunctions.h"
#include "/Users/javo/libreria/ranlib.h"
#include "/Users/javo/libreria/rnglib.h"
#include "/Users/javo/libreria/libmul.h"

#define RAD_CENSAL 56 //56 to Tartagal, 60 to Oran
#define MAX_VECI 10 // 10 to Tartagal, 11 to Oran
#define DAYS 3287 // 9 years two 29 feb
#define WEEKS 470
int main(int argc,char *argv[])
{
	FILE *archivo;
	int i,j,k;
    int imported[DAYS];
    int week,day;
    long int semla;
    int semana[7],rep_mes[];
    /* indice semilla */
    semla = -45;
	
	  /*imported people*/
    archivo=fopen("aguasblancas.txt","w");
    
    /*0 = domingo -> 6 = sabado*/
    /*prob en la semana
    lunes a jueves = 40%
    viernes a domingo = 60% */
    semana[0] = 0.2;
    semana[1] = 0.1;
    semana[2] = 0.1;
    semana[3] = 0.1;
    semana[4] = 0.1;
    semana[5] = 0.2;
    semana[6] = 0.2;
    
    
    rep_mes = 0.25;
    
    prob = 0.25/100.;
    inc = 0;
    prevalencia = 0;
    /*days*/
    i = 0;
    /*week*/
    j = 0;
    week = 6;
    while ( i<DAYS ){
    
        if (week < i){
            printf("%d\t%d\n", j,prevalencia);
            j++;
            prevalencia = 0;
            week = week + 7;
        }
        inf_imp = prob*( bolivia[j]*imported[i] )/POPULATION;
        inc = poidev(inf_imp,&semla);
        prevalencia = prevalencia + inc;
        
        i++;
    }
	
	
	exit(0);
}
