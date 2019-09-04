/* by me 2018 */


#include <math.h>
#ifndef _ambientales_
#define _ambientales_

#define Rthres 8.
#define Hmax 24.
#define kmmC 3.9E-5 // or 3.3E-6 or 6.6E-5


float egg_wet(float rain){
    return 0.8*pow(rain/Rthres,5.)/( 1. + pow(rain/Rthres,5.) );
    }

float C_Gillet(float L, float KL){
    float salida;
    if (L>KL){salida = 1. - L/KL;}
    else {salida = 0.;}
    return salida;
    }

float H_t_1(float rain, float Hum, float H_t, float Temp){
    float aux,salida;
    //DELTA_t = rain - k*(25. + Temp)*(25. + Temp)*(100. - Hum);
    aux = H_t + rain - kmmC*(25. + Temp)*(25. + Temp)*(100. - Hum); // H_t + DELTA_t
    if (aux <=0.){salida=0.;}
    else{
        if (Hmax<= aux){salida = Hmax;}
        else{salida = aux;}
        }
    return salida;
    }


float theta_T(float Temp){
    float salida;
    if ((Temp < 32.7) && (11.7 < Temp)) {
        salida = 0.1137*(-5.4 + 1.8*Temp - 0.2124*Temp*Temp + 0.01015*Temp*Temp*Temp - 0.0001515*Temp*Temp*Temp*Temp);
        }
    else {salida = 0.;}
    return salida;
}

float rate_mx(float Tm, float DHA, float DHH,float T12){
    float salida,aux0,aux1;
    /*R gas es 8.314472 (J) \'o 1.987207 (cal)*/
    aux0 = (DHA/1.987207) * (1./298. - 1./(Tm + 273.15) );
    aux1 = (DHH/1.987207) * (1./T12 - 1./(Tm + 273.15) );
    salida = kmmC*( (Tm + 273.15)/298. ) * exp( aux0 )  / (1.+exp(aux1));
    return salida;
    }

float pupa(float v[],float ED,float alpha1){
    float mL,mP,muP,mE,CG,muL,CL,muE;
    mL=v[0];mP=v[1];muP=v[2];mE=v[3];CG=v[4];muL=v[5];CL=v[6];muE=v[7];
    return  (mL / (mP + muP))* ( mE*CG / (muL+CL + mL) ) * ( alpha1/(mE*CG+muE) )*ED ;
    }

#endif 


