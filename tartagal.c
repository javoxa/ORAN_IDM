#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "/Users/javo/libreria/aleatorios_mio.h"
#include "/Users/javo/libreria/SpecialFunctions.h"
//#include "/Users/javo/libreria/ambientales.h"
/*
Ojo! En este modelo solo tengo los datos diarios de humedad, que a diferencia 
* Atettion! In this model only 
de oran si los tenia, por eso aqui la humedad tambien va por dias
*/

	/*macross seven trash*/
 
	/*seir - si : Tartagal, rain, temp, patch censal */
             
/* from  Jan 1, 2009 to Dec 31, 2017 */

//cantidad de iteraciones
#define ITERAR 75
//limite de control biologico
#define LIMPIA 30
#define LARVICIDA 100
//indices para control biologico
#define INDICE_Larv 5.
#define INDICE_Pup 5.
//efectividad del control biologico
#define EFECTIV 0.5 // 0.4 ~ 0.6 Dorso paper mosquitos
#define EFIC_LP 0.5 // eficiencia de larvas y pupas

/*paramites aedes*/
#define Kmax 212.//212. // or 170. or 400.  per house
#define beta_day 0.35 // or 0.25 or 1.19 or 0.47 or 0.5
#define SUV 0.5 // lo que sobrevive q esta en agua cuando T<10 C
#define Cond_I 10 // condiciones a todos estados

#define Rthres 12.5 // 7.5 or 10.0 or 12.5
#define Hmax 24.
#define kmmC 3.9E-5 //3.3E-6 // or 3.3E-6 or 6.6E-5 or 3.9E-5

//parametros tartagal
#define POPULATION 76000 //tartagal, 89000 to Oran
#define RAD_CENSAL 56 //56 to Tartagal, 60 to Oran
#define MAX_VECI 10 // 10 to Tartagal, 11 to Oran
#define DAYS 3287 // 9 years two 29 feb
#define WEEKS 470
#define BEGIN_YEAR 2009
#define YEARS 9 // 2009 -> 2017 = 9 years
//#define HORAH 78888//55536
#define Dt 1.E-0


/* prevalencia cero y population*/

#define ALPHA 0.75
#define MIObh 0.75//0.0521//0.0521*2. //0.0521 //1. //
#define MIObv 0.75//0.4867//0.4867*2. //0.4867 //1. //
/* tasa de gente importada infectada*/
#define RATE_IMPORT (1./100.)*0.4

float egg_wet(float rain);
float C_Gillet(float L, float KL);
float H_t_1(float rain, float Hum, float H_t, float Temp);
float theta_T(float Temp);
float rate_mx(float Tm, float DHA, float DHH,float T12);

int main(int argc,char *argv[])
{
	struct timespec timeStamp, oldTimeStamp;
	clock_gettime(CLOCK_REALTIME, &oldTimeStamp);
	
	static long int seed,seed0;	
	 seed = -atol(argv[1]);
	 seed0 = -atol(argv[2]);


      FILE *archivo;
      
      /*year, moth, week, day*/
      int day = 0,
      year = BEGIN_YEAR,
      week = 0,
      month = 365.;
      
      float inc_y=0.,
      inc_yp=0.;
      
      /*count */
      int i=0,j,kk,gg,hi,hf;
      float dt_W=Dt;
      //ambientales en dia
      static float Tmax[DAYS],Tmin[DAYS],Rain[DAYS];
      
      /*esto es para tartagal*/
      static float Hum_max[DAYS],Hum_min[DAYS];
      float rate_imp[YEARS];
      float inf_import;
      float inf_import_bol;
      int semana,se;
      
      /************************************************************/
      /*imported*/
      float Hog[RAD_CENSAL],n_cont[RAD_CENSAL];
      float importados[DAYS],imported[DAYS];
      float rate_week_imp[WEEKS];
      float aux0,aux1,aux2,aedes,perturbado=0.,MV,MH;
      float av_Hum,av_Temp,Tdeath;
      
      /* neighborhood contacts */
      int cc[RAD_CENSAL][MAX_VECI];
      int sorteo[RAD_CENSAL][RAD_CENSAL];
      int limpieza = LIMPIA;
      float inter_larv_pup=0.,inter_mosco=0.,cacharro;
	
	  /*variables for each census radius*/
	  float H_t[RAD_CENSAL];
	  float ED[2][RAD_CENSAL],EW[2][RAD_CENSAL],Larv[2][RAD_CENSAL],Pupa[2][RAD_CENSAL],MOSCO[2][RAD_CENSAL];
	  float H[2][RAD_CENSAL],Vs[2][RAD_CENSAL],Vi[2][RAD_CENSAL],Hs[2][RAD_CENSAL],Hi[2][RAD_CENSAL],He[2][RAD_CENSAL],Hr[2][RAD_CENSAL];//,MH[2][64],MV[2][64];
      float ED0[RAD_CENSAL],EW0[RAD_CENSAL],Larv0[RAD_CENSAL],Pupa0[RAD_CENSAL],MOSCO0[RAD_CENSAL];
      float H_t0[RAD_CENSAL];
      
      /*paramites de chimpokemon*/
	  float mL[RAD_CENSAL],mP[RAD_CENSAL],mE[RAD_CENSAL],muL[RAD_CENSAL],muP[RAD_CENSAL],muE[RAD_CENSAL],CG[RAD_CENSAL],CL[RAD_CENSAL],KL[RAD_CENSAL],muM[RAD_CENSAL];
	  
	  /*paramites sir - si*/
	  float m,c,bh,bv,gh,inc,sh;
      int ant=0;    
      int corr,year_epid=0;     
      float inc_ya[YEARS],inc_ypc[YEARS];
      float save_inc_ya[YEARS][ITERAR],save_inc_ypc[YEARS][ITERAR],inter[2][YEARS][ITERAR];
      
      /*armo el vector perturbacion de agentes */
      float inf_ext[DAYS][RAD_CENSAL];
      
      /*Archive read*/
         
      /*array census radius used ; as array*/
	  archivo=fopen("censo_tar.txt","r");
      for( j=0; j<RAD_CENSAL; j++ ){fscanf(archivo,"%f\n",&H[ant][j]);H[!ant][j]=H[ant][j];}
      fclose(archivo);
    
      /*vecinos de cada radio*/
      archivo=fopen("contactos_vecinos.txt","r");
      for( j=0; j<RAD_CENSAL; j++ ){for(hi=0;hi<MAX_VECI;hi++){fscanf(archivo,"%d\t",&cc[j][hi]);}fscanf(archivo,"\n");}
	  fclose(archivo);
    
      /*Hogares por radio censal*/
	  archivo=fopen("hogares_tar.txt","r");
      for( j=0; j<RAD_CENSAL; j++ ){fscanf(archivo,"%f\n",&Hog[j]);}
	  fclose(archivo);

      /*imported people*/
	  archivo=fopen("yacuiba.txt","r");
      for(i=0;i<DAYS-1;i++) {fscanf(archivo,"%f\n",&imported[i]);}
	  fclose(archivo);
	  
	  /*Contactos de vecinos*/
	  archivo=fopen("num_contactos_tar.txt","r");
      for( j=0; j<RAD_CENSAL; j++ ){fscanf(archivo,"%f\n",&n_cont[j]);}
	  fclose(archivo);

	  /*infection rate import per years */
      archivo=fopen("rate_import.txt","r");
      for(i=0;i<YEARS;i++){fscanf(archivo,"%f\n",&rate_imp[i]);}
      fclose(archivo); 

      /*infection rate import per week 0~1.*/
      archivo=fopen("serie_bolivia.txt","r");
      for(i=0;i<WEEKS;i++){fscanf(archivo,"%f\n",&rate_week_imp[i]);}
      fclose(archivo);
    
      /*eviromental values*/
      archivo=fopen("Tartago_2009_2017.txt","r");
      for(i=0;i<DAYS;i++) {fscanf(archivo,"%f\t%f\t%f\t%f\t%f\n",&Tmax[i],&Tmin[i],&Rain[i],&Hum_max[i],&Hum_min[i]);}
	  fclose(archivo);

	  /* initial conditions censal radius*/
	  archivo=fopen("bichos.txt","r");
      for(j=0;j<RAD_CENSAL;j++){fscanf(archivo,"%f\t%f\t%f\t%f\t%f\t%f\n",&ED0[j],&EW0[j],&Larv0[j],&Pupa0[j],&MOSCO0[j],&H_t0[j]);} //radios censales subidos
      fclose(archivo);
      
      /*archivo de sorteo para radios al azar*/
      archivo=fopen("sorteo.txt","r");
	  for( i=0; i<RAD_CENSAL; i++ ){
		  for( j=0; j<RAD_CENSAL; j++ ){
			  fscanf(archivo,"%d\t",&sorteo[i][j]);
		  }
		  fscanf(archivo,"\n");
	  }
	  fclose(archivo);
	  

      
      /* ***************************** *
       * ***************************** *
       * ***************************** *
       * Comienza el loop de los years *
       * ***************************** *
       * ***************************** *
       * */
      
      for ( corr=0; corr<ITERAR; corr++ ){
      
	  day = 0;
	  year = BEGIN_YEAR;
	  week = 0;
	  month = 365.;
	  inc_y = 0.;
	  inc_yp = 0.;
	  ant = 0;
	  inter_larv_pup = 0.;
	  inter_mosco = 0.;
	  year_epid = 0.;
	  cacharro=LARVICIDA;
          
      /*Siempre pongo los importados*/
		  
	  /*imported people*/
      for(i=0;i<DAYS-1;i++){
		  importados[i] = imported[i];
		  }
	  
	  /* initial conditions*/
	  for(i=0;i<YEARS;i++){
		  inc_ya[i] = 0.;
		  inc_ypc[i] = 0.;
		  save_inc_ya[i][corr] = 0.;
		  save_inc_ypc[i][corr] = 0.;
		  }
		  
	  /* initial conditions censal radius*/
      
    
	  for(j=0;j<RAD_CENSAL;j++){
		  ED[ant][j] = ED0[j];
		  EW[ant][j] = EW0[j];
		  Larv[ant][j] = Larv0[j];
		  Pupa[ant][j] = Pupa0[j];
		  MOSCO[ant][j] = MOSCO0[j];
          ED[!ant][j] = ED[ant][j];
		  EW[!ant][j] = EW[ant][j];
		  Larv[!ant][j] = Larv[ant][j];
		  Pupa[!ant][j] = Pupa[ant][j];
		  MOSCO[!ant][j] = MOSCO[ant][j];
          H_t[j] = H_t0[j];
          }
	  
	  
    
      /*initial conditions for vectors and populations*/
      for(j=0;j<RAD_CENSAL;j++){
		  Vs[ant][j] = MOSCO[ant][j];
		  Vi[ant][j] = 0.;
		  Vs[!ant][j] = Vs[ant][j];
		  Vi[!ant][j] = 0.;
		  Hs[ant][j] = (int)(ALPHA*H[ant][j]);
		  He[ant][j] = 0.;
		  Hi[ant][j] = 0.;
		  Hr[ant][j] = H[ant][j] - Hs[ant][j];
		  Hs[!ant][j] = Hs[ant][j];
		  He[!ant][j] = He[ant][j];
		  Hi[!ant][j] = Hi[ant][j];
		  Hr[!ant][j] = Hr[ant][j];
		  }
	  
	  inter[0][year_epid][corr] = 0;
	  inter[1][year_epid][corr] = 0.;
      
	  /*paramites anuj
	   * por el momento constante en todo*/
	  inc = 0.; //incidencia para la semana.
	  m = 1.; //mosquito per human -> ya no lo uso
	  c = 1.;// bite rate -> absorbido en los bh y bv
	  bh = MIObh;//anuj paper
	  bv = MIObv;//anuj paper
	  gh = 1./7.;// una semana
	  sh = 1./7.;// exposed 1 week
		
		
      /*Armo la perturbacion de infectados*/     
      /*indice 0 para todos*/
      srand(time(NULL));
      for(j=0;j<RAD_CENSAL;j++){
		  for(i=0;i<DAYS-1;i++){
			  inf_ext[i][j]=0;
			  }
		}
	 se = 0;//semana epi
	 semana = 7;
	  for(i=0;i<DAYS-1;i++){
		  if (i>semana){
			  se++;
			  semana+=7;
			  }
		  inf_import = RATE_IMPORT*(rate_week_imp[se]*importados[i])/POPULATION;
		  inf_import_bol = poidev(dt_W*inf_import,&seed0);
		  
		  /*Sorteo en un radio al azar*/
		  
		  while(0<inf_import_bol){
			  j = (int)( (RAD_CENSAL-1)*( rand()/(RAND_MAX+1.0) ) );
			  inf_ext[i][j]++;
			  inf_import_bol--;
			  }
	    }
	    
      
      
      
     
	  srand(time(NULL));//printf("Here!!! %d\n",year_epid);
	  // start jaunery , 2009
	  for ( i=0; i<DAYS-1; i++ ){
		  day++;
          /*********************************************/
          /*Esto es para tartago*/
          /*******************************************/
          
          Tdeath = Tmin[i];
          av_Temp = 0.5*(Tmin[i] + Tmax[i]);
          av_Hum = 0.5*(Hum_max[i] + Hum_min[i]);
          
          /*******************************************/
		  
		  /*Recorro todo los radio censales*/
		   
		  for(j=0; j<RAD_CENSAL; j++){
		  /*cosas del mosco en funcion de las ambientales KL lo divido por la cantidad de hogares*/
           
                        /*Aedes populations*/
          aux0 = Rain[i];
          aux1 = H_t[j];
		  H_t[j] = H_t_1(aux0,av_Hum,aux1,av_Temp);
		  KL[j] = Hog[j]*Kmax*H_t[j]/Hmax + 1.;
		  aux0 = Larv[ant][j];
		  aux1 = KL[j];
		  CG[j] = C_Gillet(aux0,aux1);
		  CL[j] = 1.5*(aux0/aux1);
          
		  /*maturations rates*/
		  
		  mE[j] = 0.24*rate_mx(av_Temp, 10798.,100000.,14184.);
		  
		  mL[j] = 0.2088*rate_mx(av_Temp, 26018.,55990.,304.6);
		  
		  if(Tdeath<13.4){mL[j]=0.;} //larvario
		  
		  mP[j] = 0.384*rate_mx(av_Temp, 14931.,-472379.,148.);
		  
		  /*mortality rates*/
		  
		  muE[j] = 0.011;//0.011;
		  
		  muL[j] = 0.01 + 0.9725*exp(- (av_Temp - 4.85)/2.7035);
		  
		  muP[j] = 0.01 + 0.9725*exp(- (av_Temp - 4.85)/2.7035);
		  
		  muM[j] = 0.091;
		  
		  /*Eggs dry*/
		  
		  aux0 = ED[ant][j]*egg_wet(Rain[i]);
		  
          /*deposito de egg*/
          
		  aux1 = poidev(dt_W*beta_day*theta_T(av_Temp)*MOSCO[ant][j],&seed0);
		  aux2 = poidev(dt_W*muE[j]*ED[ant][j],&seed0);
		  ED[!ant][j] = ED[ant][j]+ (aux1 - aux2-aux0);
		  
		  if( ED[!ant][j]<0. ){ ED[!ant][j] = 0.;}
		  
		  /*Eggs wet*/
		  
		  aux1 = poidev(dt_W*muE[j]*EW[ant][j],&seed0);
		  aux2 = poidev(dt_W*mE[j]*CG[j]*EW[ant][j],&seed0);
		  EW[!ant][j] = EW[ant][j] + aux0 - aux1 - aux2;
		  
		  if ( EW[!ant][j] < 0. ){
			  
			  EW[!ant][j] = 0.;
			  aux2 = 0.5*( EW[ant][j] + aux0 );
			  
			  }
		  
		  if ( Tdeath < 10. ){
			  
			  EW[!ant][j] = SUV*EW[!ant][j];
			  
			  }
		  
		  aux0=aux2;
          
		  /*Larvitar*/
		  
		  aux1 = poidev(dt_W*(muL[j]+CL[j])*Larv[ant][j],&seed0);
		  
		  aux2 = poidev(dt_W*mL[j]*Larv[ant][j],&seed0);
		  
		  Larv[!ant][j] = Larv[ant][j] + aux0 - aux1 - aux2;
		  
		  if( Larv[!ant][j] < 0. ){
			  
			  Larv[!ant][j] = 0.;
			  aux2 = 0.5*(Larv[ant][j] + aux0);
			  
			  }
			  
		  if ( Tdeath < 10. ){
			  
			  Larv[!ant][j] = SUV*Larv[!ant][j];
			  
			  }
		 
		  aux0=aux2;
          
		  /*Pupitar*/
		  
		  aux1 = poidev(dt_W*muP[j]*Pupa[ant][j],&seed0);
		  aux2 = poidev(dt_W*mP[j]*Pupa[ant][j],&seed0);
		  Pupa[!ant][j] = Pupa[ant][j] + aux0 - aux1 - aux2;
		  
		  if( Pupa[!ant][j] < 0. ){
			  Pupa[!ant][j] = 0.;
			  aux2 = 0.5*( Pupa[ant][j] + aux0 );
			  }
			  
		  if ( Tdeath < 10. ){
			  Pupa[!ant][j] = SUV*Pupa[!ant][j];
			  }
		  
		  /*control de larvitar y pupitar */
		  
		  if (cacharro>0){
			  if ( (  Larv[!ant][j]/(float)(Hog[j]) > INDICE_Larv  ) || (  Pupa[!ant][j]/(float)(Hog[j]) > INDICE_Pup)  ) {Larv[!ant][j] = EFIC_LP*Larv[!ant][j];
				  Pupa[!ant][j] = EFIC_LP*Pupa[!ant][j];
				  inter_larv_pup++;
			     }
			  cacharro--;
		     }
		 /*fin control*/ 
		  
		  aux0 = aux2;
          
		  /*Tyranitar*/
          
		   /*campo medio con los vecinos censales
		    * el indice j es el correspondiente al numero
		    * de radio censal con el que se cuenta */
		    /* i es el dia, j es el radio censal*/
		    /*estoy agregando los infectados externos
		     * al termino de trasmision
		     * */
		     
		    MH = H[ant][j] + inf_ext[i][j] ;
		    MV = Hi[ant][j] + inf_ext[i][j];
		    
		    for( hi=0; hi<MAX_VECI; hi++ ){
				
				hf = cc[j][hi];
				
				if(hf ==! 0){
					
					MH+=H[ant][hf-1];
					MV+=Hi[ant][hf-1];
					
					}
					
				else{break;}
			}
			
			
			/*Agrego uno radio al azar que no puede ser de cc[j][hi] y 
			 * tiene que ser distinto de "j"
			 * parece que no funciona*/
			hi = n_cont[j];
			kk = RAD_CENSAL - hi - 1;
			gg = (int)( kk*ran1(&seed));
			hf = sorteo[j][gg];
			MH+=H[ant][hf-1];
			MV+=Hi[ant][hf-1];
			
			/* *******************************************************
			 * *******************************************************
			 * *******************************************************
			 * *******************************************************
			 */
			  
            /* Vectores S-I with demography */
           /*pica mosco*/
		   aux1 = poidev(theta_T(av_Temp)*dt_W*beta_day*bv*Vs[ant][j]*MV/MH,&seed); //con temperatura theta_T(av_Temp)*
           aux2 = poidev(muM[j]*Vs[ant][j],&seed);
           
		   Vs[!ant][j] = Vs[ant][j] + aux0 - aux1 - aux2;
		   
		   if( Vs[!ant][j]<0. ){
			   
            Vs[!ant][j] = 0.;
            aux1 = aux0;
            aux2 = Vs[ant][j];
            
            }
            
           aux2 = poidev(muM[j]*Vi[ant][j],&seed);
		   Vi[!ant][j] = Vi[ant][j] + ( aux1 - aux2 );
		   
		   
		   if(Vi[!ant][j]<0.){
			   
			   Vi[!ant][j] = 0.;
			   aux1 = Vi[ant][j];
			   
			   }
			   
			   
           MOSCO[ant][j] = Vs[ant][j] + Vi[ant][j];
           
		   
		   /*Hospedadores S-E-I-R*/
		   /* *******************************************************
			* *******************************************************
		    * *******************************************************
		    * *******************************************************
			*/
			  
		   MH = H[ant][j];
           MV = Vi[ant][j];
		    for(hi=0;hi<MAX_VECI;hi++){
				hf=cc[j][hi];
				if(hf==!0){
					MH+= H[ant][hf-1];
					MV+= Vi[ant][hf-1];
					}
				else{break;}
			}

			/*Agrego uno radio al azar no es ni j, ni un vecino 
			 * sorteo[i][gg] es un arreglo de sorteo
			 * */
			hi = n_cont[j];
			kk = RAD_CENSAL - hi - 1;
			gg = (int)( kk*ran1(&seed));
			hf = sorteo[j][gg];
			MH+=H[ant][hf-1];
			MV+=Vi[ant][hf-1];
            
            /*pica mosco*/
		   aux0=poidev(theta_T(av_Temp)*beta_day*bh*MV*Hs[ant][j]/MH,&seed); // con temp theta_T(av_Temp)*
		   
		   Hs[!ant][j] = Hs[ant][j] - aux0;
		   
		   if( Hs[!ant][j] < 0. ) {
			   
			   Hs[!ant][j] = 0.;
			   aux0 = Hs[ant][j];
			   }
			   
		   aux1 = poidev(sh*He[ant][j],&seed);
		   
		   He[!ant][j] = He[ant][j] + aux0 - aux1;
		   
		   if(He[!ant][j]<0.){
			   
               He[!ant][j] = 0.;
			   aux1 = He[ant][j] +aux0;
			   
               }
		   
		   aux1 = poidev(dt_W*gh*Hi[ant][j],&seed);
		   
           /* inf_import*inf_ext[i][j] es el infectado externo*/
           aux2 = inf_ext[i][j];
           
		   Hi[!ant][j] = Hi[ant][j] + aux0 - aux1 + inf_ext[i][j];
		   if(Hi[!ant][j]<0.){Hi[!ant][j]=aux2;aux1=aux0+Hi[ant][j];}		   
		   /*control biologico por reduccion de moscos*/
		   
		   if((aux1>0) || (aux2>0)){
			   
			   if(0<limpieza){
				   
				   MOSCO[!ant][j] = EFECTIV*MOSCO[!ant][j];
				   Vs[!ant][j] = EFECTIV*Vs[!ant][j];
				   Vi[!ant][j] = EFECTIV*Vi[!ant][j];
				   Larv[!ant][j] = EFECTIV*Larv[!ant][j];
				   Pupa[!ant][j] = EFECTIV*Pupa[!ant][j];
				   limpieza--;
				   inter_mosco++;
				   }
			   }
		  
		   /*end of biological control fin de control biologico*/
		   
		   
		   aux0 = aux1;
		   H[!ant][j] = H[ant][j] +inf_ext[i][j];
		   perturbado+=inf_ext[i][j];
		   inc_yp+=inf_ext[i][j];
		   Hr[!ant][j] = H[!ant][j] - Hi[!ant][j] - He[!ant][j] - Hs[!ant][j];
		   
		  inc+=aux0;//
		  
		  inc_y+=aux0;
		  	  
		//if (day == 7){aedes=aedes+MOSCO[!ant][j];}
	  }
	  
	  if (day == 7){
			  week++;aedes=0.;day=0;inc=0.;perturbado=0.;
			}
			
	  if ( i>=month ){ //month no es month es un contador para los dias
		  
		  inc_ya[year_epid] = inc_ya[year_epid]+inc_y;
		  //if (year_epid==1){printf("%.2f\n",inc_ya[year_epid]);}
		  
		  inc_ypc[year_epid] = inc_ypc[year_epid]+inc_yp;
		  
		  save_inc_ya[year_epid][corr] = inc_y;
		  
		  save_inc_ypc[year_epid][corr] = inc_yp;
		  
		  inc_y = 0.;
		  
		  inc_yp = 0.;
		  
		  year++;
		  
		  year_epid++;
		  
		  month+=365;
		  
		  limpieza = LIMPIA;
		  cacharro = LARVICIDA;
		  
		  inter[0][year_epid][corr] = inter_larv_pup;
		  inter[1][year_epid][corr] = inter_mosco;
		  inter_larv_pup = 0.;
		  inter_mosco = 0.;
		  }
		  
		  /*if ((year_epid<2010) && (corr> (int)(ITERAR - 2 ) )){
			 // printf("larv=%.2f\tmosco=%f\tcorr=%d\n",inc_ypc[year_epid],inc_ya[year_epid],corr);
			  
		  }*/
			 
	  ant=!ant;	  
	  if (YEARS<year_epid){year_epid=0;}
	  
	  }
	  
  } printf("tartagal_model\n");
  

  printf("#year\t<inc>\tS_inc\t<imp>\tS_imp\tr_imp\t<larv>\t<mosq>\n");
  aux2=(float)ITERAR;

  for(i=0;i<YEARS;i++){
	  aux0=0.;aux1=0.;inc_y=0.;inc_yp=0.;MV=0.;MH=0.;
	  for(j=0;j<ITERAR;j++){
		  inc_y = inc_y + save_inc_ya[i][j];
		  inc_yp = inc_yp + save_inc_ypc[i][j];
		  MH = MH + inter[0][i][j];
		  MV = MV + inter[1][i][j];
	  }
	  
	  /*mean value*/
	  
	  inc_y = inc_y/aux2;
	  inc_yp = inc_yp/aux2;
	  MV = MV/aux2;
	  MH = MH/aux2;
	  
	  /* Calculo de la desviacion estandar */
	  
	  for(j=0;j<ITERAR;j++){
		  aux0 = aux0 + (save_inc_ya[i][j] - inc_y)*(save_inc_ya[i][j] - inc_y);
		  aux1 = aux1 + (save_inc_ypc[i][j] - inc_yp)*(save_inc_ypc[i][j] - inc_yp);
	  }
	  
	  /* ******** error absoluto por orde de -> JPA ******* */
	  /* Para mi es mas representativo la varianza pero bue... */
	  
	  aux0 = sqrt(aux0/( (aux2-1.)*(aux2) ) );
	  aux1 = sqrt(aux1/( (aux2-1.)*(aux2) ) );
	  
	  printf("%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",BEGIN_YEAR+i,inc_y,aux0,inc_yp,aux1,rate_imp[i],MV,MH);
	  
  }
	  
	  printf("Iterations per year %d\n",ITERAR);
	  printf("bite rate per day %.2f\n",beta_day);
	  printf("larvicida max %d\n",LARVICIDA);
	  printf("descacharrado max %d\n",LIMPIA);
	  printf("Intervencion para larva %.1f\n",INDICE_Larv);
	  printf("Intervencion para pupa %.1f\n",INDICE_Pup);
	  clock_gettime(CLOCK_REALTIME, &timeStamp);
      printf("Tiempo de computo = %ld seg\n", 1000 * (timeStamp.tv_sec - oldTimeStamp.tv_sec) + (timeStamp.tv_nsec - oldTimeStamp.tv_nsec)/ 1000000000);
      
      
exit(0);
}


float egg_wet(float rain){
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
        if (Hmax<= aux){salida = Hmax;}
        else{salida = aux;}
        }
    return salida;
    }


float theta_T(float Temp){
    float salida;
    if ( (11.7 < Temp) && (Temp < 32.7)) {
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
    salida = ( (Tm + 273.15)/298. ) * ( exp( aux0 )  / (1.+exp(aux1)) );
    return salida;
    }

