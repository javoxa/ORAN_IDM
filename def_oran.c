#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//cantidad de iteraciones
#define     ITERAR              50   //minimo ->1

//limite de control biologico
#define     LIMPIA              0
#define     LARVICIDA           0

//indices para control biologico
#define     INDICE_Larv         1.
#define     INDICE_Pup          1.

//efectividad del control biologico
#define     EFECTIV             0.4 // 0.4 ~ 0.6 Dorso paper mosquitos
#define     EFIC_LP             0.4 // eficiencia de larvas y pupas
#define     MATAR_VECTORES      10.//12.5//temp
#define     NO_INFECCION        19.75 //temp

/*paramites aedes en days */
#define     EGG_LIFE            90. // vida de los huevos [dias]

#define     beta_day            2.69 // or 0.25 or 1.19 or 0.47 or 0.5
#define     bite_rate           1.79 // 2.59

#define     MU_MOSQUITO_JOVEN   1./4. //0.091//1./10.//0.091 // 1./2.;//0.091;
          
#define     MADURACION_MOSQUITO 1./4. //1./4.;//maduracion
          
#define     MU_MOSQUITA_ADULTA  1./11.;// vida^{-1} del vector en optimas condiciones (22 grados) 0.091; //

#define     Ve_Lat              8 // cantidad de latentes mosquitos (variar Rate_VE()

#define     Kmax                75.//212. // or 170. or 400.  per house

#define     DIFUSION_VECTORES   0.10

#define     SUV                 0.6 // lo que sobrevive q esta en agua cuando T<10 C
#define     Cond_I              10 // condiciones a todos estados

#define     LAG                 7 // lag para actuar con el larvicida

#define     AV_EF_TEMP          0.5 // valor para el promedio de la temperatura
#define     AV_EF_HUM           0.5 //valor promedio de la temperatura
#define     Rthres              10.5 // 7.5 or 10.0 or 12.5
#define     Hmax                24.
#define     kmmC                3.3E-6//6.6E-5//3.9E-5//6.6E-5 //3.3E-6 // or 3.3E-6 or 6.6E-5 or 3.9E-5

//parametros tartagal
#define     POPULATION          75697 //6999 tartagal, 75697 to Oran
#define     RAD_CENSAL          60 //56 to Tartagal, 60 to Oran
#define     MAX_VECI            11 // 10 to Tartagal, 11 to Oran
#define     DAYS                3287 // 9 years two 29 feb
#define     WEEKS               470
#define     BEGIN_YEAR          2009
#define     YEARS               9 // 2009 -> 2017 = 9 years

//#define HORAH 78888//55536
#define     Dt 1.E-0

/* prevalencia cero y population*/
#define     ALPHA               0.85
#define     MIObh               0.75//0.0521*2. //0.0521 //1. //
#define     MIObv               0.75//0.4867*2. //0.4867 //1. //
#define     Remove_infect       7. //dias
#define     Remove_expose       7. //dias

/* tasa de gente importada infectada*/
#define     RATE_IMPORT         1./150.

float egg_wet( float rain ){

    return 0.8*pow(rain/Rthres,5.)/( 1. + pow(rain/Rthres,5.) );
    
    }

float C_Gillet(float L, float KL){

    float salida;
    
    if (L<KL){salida = 1. - L/KL;}
    else {salida = 0.;}
    
    return salida;
    
    }

float H_t_1(float rain, float Hum, float H_t, float Temp){

    float aux,salida;
    //DELTA_t = rain - k*(25. + Temp)*(25. + Temp)*(100. - Hum);
    aux = H_t + rain - kmmC*(25. + Temp)*(25. + Temp)*(100. - Hum); // H_t + DELTA_t
    
    if (aux <=0.){salida=0.;}
    else{
        if (Hmax<= aux){
            salida = Hmax;
        }
        else{
            salida = aux;
        }
    }
    
    return salida;
    }


float theta_T( float Temp) {
    float salida;
    
    if ( (11.7 < Temp) && (Temp < 32.7)) {
        salida = 0.1137*(-5.4 + 1.8*Temp - 0.2124*Temp*Temp + 0.01015*Temp*Temp*Temp - 0.0001515*Temp*Temp*Temp*Temp);
    }
    else {
        salida = 0.;
    }
    
    return salida;
}

float rate_mx( float Tm, float DHA, float DHH, float T12) {

    float salida,aux0,aux1;
    /*R gas es 8.314472 (J) \'o 1.987207 (cal)*/
    
    aux0 = (DHA/1.987207) * (1./298. - 1./(Tm + 273.15) );
    
    aux1 = (DHH/1.987207) * (1./T12 - 1./(Tm + 273.15) );
    
    salida = ( (Tm + 273.15)/298. ) * ( exp( aux0 )  / (1.+exp(aux1)) );
    
    return salida;
}

float rate_Ve(float Tm, float k){

    float salida, theta;
    if ( k <= 3.)   {
        /* to k == 3 ; R^2 = 1 */
        //if ( Tm > 34.) {  Tm = 34.;  }
        //theta  =  -0.0054*Tm*Tm*Tm + 0.3437*Tm*Tm - 6.6102*Tm + 42.272;
        theta  =  exp( -0.1659*Tm + 6.7031 );
    }
    else{
        /* to k = 10~7; R^2 = 1*/
        //if ( Tm > 34.) {  Tm = 35.;  }
        //theta  =  -0.0104*Tm*Tm*Tm + 0.4304*Tm*Tm - 6.1852*Tm + 35.765;
        theta  =  exp( -0.155*Tm + 6.5031 );
    }
    
    
    
    salida = k/theta;
 
    return salida;
}


float muerte_V(float Tm){

    float salida;
    float factor = 0.0360827; //factor a 22 grados
    salida = 8.692E-1 -(1.59E-1)*Tm + (1.116E-2)*Tm*Tm -(3.408E-4)*Tm*Tm*Tm + (3.809E-6)*Tm*Tm*Tm*Tm;
    return salida/factor;

}





