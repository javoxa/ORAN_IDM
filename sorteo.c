#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "/Users/javo/libreria/aleatorios_mio.h"
#include "/Users/javo/libreria/SpecialFunctions.h"

#define RAD_CENSAL 56 //56 to Tartagal, 60 to Oran
#define MAX_VECI 10 // 10 to Tartagal, 11 to Oran
#define DAYS 3287 // 9 years two 29 feb
#define WEEKS 470
#define POPULATION 74000.
int main(int argc,char *argv[])
{
	FILE *archivo;
	int i,j,k;
    int bolivia[WEEKS],imported[DAYS];
    int week,day;
    float prob,inf_imp;
    int inc,prevalencia;
    long int semla;
    /* indice semilla */
    semla = -45;
	
	  /*imported people*/
	  archivo=fopen("yacuiba.txt","r");
      for(i=0;i<DAYS-1;i++) {fscanf(archivo,"%d\n",&imported[i]);}
	  fclose(archivo);
    
      /* casos importados */
	  archivo=fopen("serie_bolivia.txt","r");
      for( j=0; j<WEEKS; j++ ){fscanf(archivo,"%d\n",&bolivia[j]);}
	  fclose(archivo);
    /*prob tart*/
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
