#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "/Users/javo/libreria/aleatorios.h"
#include "/Users/javo/libreria/SpecialFunctions.h"
#include "/Users/javo/libreria/ranlib.h"
#include "/Users/javo/libreria/rnglib.h"
#include "/Users/javo/libreria/libmul.h"



#include "def_oran.c"

/*
Ojo! En este modelo solo tengo los datos diarios de humedad, que a diferencia 
* Atettion! In this model only 
de oran si los tenia, por eso aqui la humedad tambien va por dias
*/

	/*macross seven trash*/
 
	/*seir - si : Tartagal, rain, temp, patch censal */
             
/* from  Jan 1, 2009 to Dec 31, 2017 */


int main(int argc,char *argv[])
{
	struct timespec timeStamp, oldTimeStamp;
	clock_gettime(CLOCK_REALTIME, &oldTimeStamp);
	
	static long int seed,seed0;	
	 seed = -atol(argv[1]);
	 seed0 = -atol(argv[2]);


      FILE *archivo;
      FILE *guardar;
      FILE *semanal;
    
      guardar = fopen("larvas_oran.txt","w");
      semanal = fopen("semanal_oran.txt","w");
    
      
      /*year, moth, week, day*/
      int day = 0,
      year = BEGIN_YEAR,
      week = 0,
      month = 365;
      
      float inc_y=0.,
      inc_yp=0.;
      
      /*count */
      int i,j,kk,gg,hi,hf,pp,oo=0;
      float dt_W=Dt;
      //ambientales en dia
      float Tmin[DAYS],
           Tmax[DAYS],
           Tmean[DAYS],
           Tmedian[DAYS],
           Rain[DAYS],
           Hum_mean[DAYS],
           Hum_median[DAYS];
    
      float rate_imp[YEARS];
      float inf_import;
      float inf_import_bol;
      int semana,se;
      
      /************************************************************/
      /*imported*/
      float Hog[RAD_CENSAL],n_cont[RAD_CENSAL];
      float importados[DAYS],imported[DAYS];
      float rate_week_imp[WEEKS];
      float aux0,aux1,aux2,aux3,aux4,MV,MH;
      float av_Hum,av_Temp,Tdeath;
      
      /* neighborhood contacts */
      int cc[RAD_CENSAL][MAX_VECI];
      int sorteo[RAD_CENSAL][RAD_CENSAL];
      int inter_larv_pup,inter_mosco;
	
	  /*variables for each census radius*/
	  float H_t[RAD_CENSAL];
	  float ED[2][RAD_CENSAL],EW[2][RAD_CENSAL],Larv[2][RAD_CENSAL],Pupa[2][RAD_CENSAL],MOSCO[2][RAD_CENSAL];
	  float H[2][RAD_CENSAL],Vector[2][RAD_CENSAL],Vs[2][RAD_CENSAL],Vi[2][RAD_CENSAL],Hs[2][RAD_CENSAL],Hi[2][RAD_CENSAL],He[2][RAD_CENSAL],Hr[2][RAD_CENSAL];//,MH[2][64],MV[2][64];
      float ED0[RAD_CENSAL],EW0[RAD_CENSAL],Larv0[RAD_CENSAL],Pupa0[RAD_CENSAL],MOSCO0[RAD_CENSAL];
      float H_t0[RAD_CENSAL];
      int delay_LP[RAD_CENSAL][DAYS];
      
      /*paramites de chimpokemon*/
	  float mL[RAD_CENSAL],mP[RAD_CENSAL],mE[RAD_CENSAL],muL[RAD_CENSAL],muP[RAD_CENSAL],muE[RAD_CENSAL],CG[RAD_CENSAL],CL[RAD_CENSAL],KL[RAD_CENSAL],muM[RAD_CENSAL];
	  
	  /*paramites sir - si*/
	  float m,c,bh,bv,gh,inc,incimp,sh;
      int ant=0;    
      int corr,year_epid;
      //float inc_ya[YEARS],inc_ypc[YEARS];
    
      float save_inc_ya[YEARS],save_inc_ypc[YEARS],
      inter_LP[YEARS],
      inter_M[YEARS];
      
      /*armo el vector perturbacion de agentes */
      float inf_ext[DAYS][RAD_CENSAL];
    
      /*para la latencia y dispersion*/
      int st_mn,im,num_ev;
      float pr,pr_mn[MAX_VECI];
      int *VectorState;
      float Ve[2][RAD_CENSAL][Ve_Lat];
      float sv;
      float muMad[RAD_CENSAL];
      float muV[RAD_CENSAL];
      float VD[RAD_CENSAL];
      float D_Vector = DIFUSION_VECTORES; //difusion
      float RVs[RAD_CENSAL];
      int st;
    
      
      /*Archive read*/
    
      /*array census radius used ; as array*/
	  archivo = fopen("censo_oran.txt","r");
      for(  j = 0; j < RAD_CENSAL; j++  ){
        fscanf( archivo, "%f\n", &H[ant][j] );
        H[!ant][j] = H[ant][j];
      }
      fclose(archivo);
    
      /*vecinos de cada radio*/
      archivo = fopen("contactos_vecinos_oran.txt","r");
      for(  j = 0; j < RAD_CENSAL; j++  ){
        for( hi = 0; hi < MAX_VECI; hi++){
            fscanf( archivo, "%d\t", &cc[j][hi]);
            }
            fscanf(archivo,"\n");
        }
	  fclose(archivo);
    
      /*Hogares por radio censal*/
	  archivo = fopen("hogares_oran.txt","r");
      for(  j = 0; j < RAD_CENSAL; j++  ){
        fscanf( archivo, "%f\n", &Hog[j]);
        }
	  fclose(archivo);

      /*imported people*/
	  archivo = fopen("yacuiba.txt","r");
      for(i=0;i<DAYS;i++) {fscanf(archivo,"%f\n",&imported[i]);}
	  fclose(archivo);
	  
	  /*Contactos de vecinos*/
	  archivo = fopen("num_contactos_oran.txt","r");
      for( j=0; j<RAD_CENSAL; j++ ){fscanf(archivo,"%f\n",&n_cont[j]);}
	  fclose(archivo);

	  /*infection rate import per years */
      archivo = fopen("rate_import.txt","r");
      for(i=0;i<YEARS;i++){fscanf(archivo,"%f\n",&rate_imp[i]);}
      fclose(archivo); 

      /*infection rate import per week 0~1.*/
      archivo = fopen("serie_bolivia.txt","r");
      for(i=0;i<WEEKS;i++){fscanf(archivo,"%f\n",&rate_week_imp[i]);}
      fclose(archivo);
    
      /*eviromental values*/
      archivo = fopen("Oran_2009_2017.txt","r");
      for(i=0;i<DAYS;i++) {fscanf(archivo,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",&Tmin[i],&Tmax[i],&Tmean[i],&Tmedian[i],&Rain[i],&Hum_mean[i],&Hum_median[i]);}
	  fclose(archivo);

	  /* initial conditions censal radius*/
	  archivo = fopen("bichos.txt","r");
      for(j=0;j<RAD_CENSAL;j++){fscanf(archivo,"%f\t%f\t%f\t%f\t%f\t%f\n",&ED0[j],&EW0[j],&Larv0[j],&Pupa0[j],&MOSCO0[j],&H_t0[j]);} //radios censales subidos
      fclose(archivo);
    
      /*archivo de sorteo para radios al azar*/
      archivo = fopen("sorteo_oran.txt","r");
	  for( i=0; i < RAD_CENSAL; i++ ){
		  for( j=0; j < RAD_CENSAL; j++ ){
			  fscanf(archivo,"%d\t",&sorteo[i][j]);
		  }
		  fscanf(archivo,"\n");
	  }
	  fclose(archivo);
	  
      	  /* initial conditions TO SAVE*/
	  
        for( i = 0; i < YEARS -1 ; i++){
			
		  save_inc_ya[i]  =  0.;
		  save_inc_ypc[i] =  0.;
          inter_LP[i]     =  0.;
          inter_M[i]      =  0.;
          
		  }
      
    
      	  /*paramites anuj
	   * por el momento constante en todo*/
	  //inc = 0.; //incidencia para la semana.
      
      
	  bh = MIObh;//anuj paper
	  bv = MIObv;//anuj paper
	  gh = 1./Remove_infect;
	  sh = 1./Remove_expose;
      //sv = 1./5.;// five days to mosquito

    
      /* ***************************** *
       * ***************************** *
       * ***************************** *
       * Comienza el loop de los years *
       * ***************************** *
       * ***************************** *
       * */
      corr = 0;
      //while( corr < ITERAR ){
      //printf("%d\t",corr);
	  day = 0;
	  year = BEGIN_YEAR;
	  week = 0;
	  month = 365.;
	  inc_y = 0.;
	  inc_yp = 0.;
      inc = 0.;
      incimp = 0.;
	  ant = 0;
	  inter_larv_pup = 0;
	  inter_mosco = 0;
	  year_epid = 0;
          
      /*Siempre pongo los importados*/
		  
	  /*imported people*/
      for (  i = 0; i < DAYS ; i++  ) {
		  importados[i] = imported[i];
		  }

          /* indice 0 larvicida */
          for ( i = 0; i < DAYS; i++ ){
                for ( j = 0; j < RAD_CENSAL; j++ ){
                delay_LP[j][i] = 0;
                }
          }
		  
	  /* initial conditions censal radius*/
      
    
	  for(  j = 0; j < RAD_CENSAL; j++ ){
      
		  ED[ant][j]        =   ED0[j];
		  EW[ant][j]        =   EW0[j];
		  Larv[ant][j]      =   Larv0[j];
		  Pupa[ant][j]      =   Pupa0[j];
		  MOSCO[ant][j]     =   MOSCO0[j];
          ED[!ant][j]       =   ED[ant][j];
		  EW[!ant][j]       =   EW[ant][j];
		  Larv[!ant][j]     =   Larv[ant][j];
		  Pupa[!ant][j]     =   Pupa[ant][j];
		  MOSCO[!ant][j]    =   MOSCO[ant][j];
          H_t[j]            =   H_t0[j];
		  Vs[ant][j]        =   (int)(0.5*MOSCO[ant][j]);
		  Vi[ant][j]        =   0.;
		  Vs[!ant][j]       =   Vs[ant][j];
		  Vi[!ant][j]       =   0.;
		  Hs[ant][j]        =   (int)(ALPHA*H[ant][j]); //usar random
		  He[ant][j]        =   0.;
		  Hi[ant][j]        =   0.;
		  Hr[ant][j]        =   H[ant][j] - Hs[ant][j];
		  Hs[!ant][j]       =   Hs[ant][j];
		  He[!ant][j]       =   He[ant][j];
		  Hi[!ant][j]       =   Hi[ant][j];
		  Hr[!ant][j]       =   Hr[ant][j];
          Vector[ant][j]    =   Vs[ant][j];
          Vector[!ant][j]   =   Vector[ant][j];
          
          // exposed vector
          for ( st = 0; st < Ve_Lat; st++){
          
                Ve[ant][j][st]  = 0.;
                Ve[!ant][j][st] = 0.;
          }
      }
		
      /*Armo la perturbacion de infectados*/     
      /*indice 0 para todos*/
      srand(time(NULL));
      for( j = 0; j < RAD_CENSAL; j++ ){
		  for( i = 0; i < DAYS-1 ; i++){
			  inf_ext[i][j] = 0;
			  }
		}
        
	 se = 0;//semana epi
	 semana = 7;
     
	  for( i = 0; i < DAYS - 1; i++ ){
		  if ( i > semana ){
			  se++;
			  semana+=7;
              //printf("%d\t%d\t%.0f\n",corr,se,importados[i]);
			  }
          
		  inf_import = RATE_IMPORT*(rate_week_imp[se]*importados[i])/POPULATION;
  //       if (i == DAYS - 1) { printf("%d\t%d\t%.0f\n",corr,i,importados[i]);
		  inf_import_bol = poidev( dt_W*inf_import, &seed0 );
          
     
		  
		  /*Sorteo en un radio al azar*/
		//printf("antes %d\t%d\t%d\n",corr,i,inf_import_bol);
		  while( inf_import_bol > 0 ){
          
			  j = (int)( (RAD_CENSAL-1)*( rand()/(RAND_MAX+1.0) ) );
			  inf_ext[i][j]++;
			  inf_import_bol--;
              //printf("%d\t%d\t%d\n",corr,i,inf_import_bol);
			  }
          
	    }
      
        
	  srand(time(NULL));
	  // start jaunery , 2009
	  for ( i = 0; i < DAYS; i++ ){
      //printf(" sale i = %d, j = %d\n",i,week);
		  day = day + 1;
          /*********************************************/
          /*Esto es para tartago*/
          /*******************************************/
          
          Tdeath  = Tmin[i];
          /*Esto cambia con oran que tengo por hs y puedo
                calcular las medias, medianas, etc...*/
          //av_Temp = AV_EF_TEMP*(Tmin[i] + Tmax[i]);
          //av_Temp = Tmedian[i]; // con la mediana
          av_Temp = Tmean[i]; // con la media
          
          //av_Hum  = Hum_median[i]; // con la mediana
          av_Hum  = Hum_mean[i]; // con la media
          
          //av_Hum  = AV_EF_HUM*(Hum_max[i] + Hum_min[i]);
          
          /*******************************************/
		  //printf ("%d\t%d\t%d\n",i,day,week);
		  /*Recorro todo los radio censales*/
		   
		  for( j = 0; j < RAD_CENSAL; j++ ){
		  /*cosas del mosco en funcion de las ambientales KL lo divido por la cantidad de hogares*/
           
                        /*Aedes populations*/
          aux0 = Rain[i];
          aux1 = H_t[j];
          
		  H_t[j] = H_t_1( aux0, av_Hum, aux1, av_Temp );
		  KL[j]  = Hog[j]*Kmax*H_t[j]/Hmax + 1.;
          
		  aux0 = Larv[ant][j];
		  aux1 = KL[j];
          
		  CG[j] = C_Gillet( aux0, aux1 );
		  CL[j] = 1.5*(aux0/aux1);
          
		  /*maturations rates*/
		  
		  mE[j] = 0.24*rate_mx(av_Temp, 10798.,100000.,14184.);
		  
		  mL[j] = 0.2088*rate_mx(av_Temp, 26018.,55990.,304.6);
		  
		  if( Tdeath < 13.4 ){  mL[j] = 0.; } //larvario
		  
		  mP[j] = 0.384*rate_mx(av_Temp, 14931.,-472379.,148.);
		  
		  /*mortality rates*/
		  
		  muE[j] = 1./EGG_LIFE;//0.008355;//0.011;
		  
		  muL[j] = 0.01 + 0.9725*exp(- (av_Temp - 4.85)/2.7035);
		  
		  muP[j] = 0.01 + 0.9725*exp(- (av_Temp - 4.85)/2.7035);
		  
		  muM[j] = MU_MOSQUITO_JOVEN ;//+ 0.9725*exp(- (av_Temp - 4.85)/2.7035) ;// 1./2.;//0.091;
          
          muMad[j] = MADURACION_MOSQUITO ;//+ 0.384*rate_mx(av_Temp, 14931.,-472379.,148.);//1./4.;//maduracion
          
          muV[j] = MU_MOSQUITA_ADULTA ;//+ 0.9725*exp(- (av_Temp - 4.85)/2.7035) ;//1./21.;// vida del vector
		  
		  /*Eggs dry*/
		  
		  aux0 = ED[ant][j]*egg_wet(Rain[i]);
		  
          /*deposito de egg*/
          
		  aux1 = poidev( dt_W*bite_rate*beta_day*theta_T(av_Temp)*Vector[ant][j], &seed0);
          
		  aux2 = poidev( dt_W*muE[j]*ED[ant][j], &seed0);
          
		  ED[!ant][j] = ED[ant][j] + (aux1 - aux2 - aux0);
		  
		  if( ED[!ant][j]<0. ){ ED[!ant][j] = 0.;}
		  
		  /*Eggs wet*/
		  
		  aux1 = poidev( dt_W*muE[j]*EW[ant][j] ,&seed0 );
          
		  aux2 = poidev( dt_W*mE[j]*CG[j]*EW[ant][j], &seed0 );
          
		  EW[!ant][j] = EW[ant][j] + aux0 - aux1 - aux2;
		  
		  if ( EW[!ant][j] < 0. ){
			  
			  EW[!ant][j] = 0.;
			  aux2 = (int)(0.5*( EW[ant][j] + aux0 ));
			  
			  }
		  
		  if ( Tdeath < 10. ){
			  
			  EW[!ant][j] = (int)( SUV*EW[!ant][j] );
			  
			  }
		  
		  aux0 = aux2;
          
		  /*Larvitar*/
		  
		  aux1 = poidev(dt_W*(muL[j]+CL[j])*Larv[ant][j],&seed0);
		  
		  aux2 = poidev(dt_W*mL[j]*Larv[ant][j],&seed0);
		  
		  Larv[!ant][j] = Larv[ant][j] + aux0 - aux1 - aux2;
		  
		  if( Larv[!ant][j] < 0. ){
			  
			  Larv[!ant][j] = 0.;
			  aux2 = (int)( 0.5*(Larv[ant][j] + aux0) );
			  
			  }
			  
		  if ( Tdeath < 10. ){
			  
			  Larv[!ant][j] = (int)( SUV*Larv[!ant][j] );
			  
			  }
		 
		  aux0 = aux2;
          
		  /*Pupitar*/
		  
		  aux1 = poidev( dt_W*muP[j]*Pupa[ant][j], &seed0 );
          
		  aux2 = poidev( dt_W*mP[j]*Pupa[ant][j], &seed0 );
          
		  Pupa[!ant][j] = Pupa[ant][j] + aux0 - aux1 - aux2;
		  
		  if( Pupa[!ant][j] < 0. ){
          
			  Pupa[!ant][j] = 0.;
			  aux2 = (int)(0.5*( Pupa[ant][j] + aux0 ));
              
			  }
			  
		  if ( Tdeath < 10. ){
			  Pupa[!ant][j] = (int)( SUV*Pupa[!ant][j] );
			  }
		  
		  /*control de larvitar y pupitar */
		  
		  if ( inter_larv_pup < LARVICIDA ){
          
                m = Larv[ant][j]/Hog[j];
                c = Pupa[ant][j]/Hog[j];
          
			  if (  m > INDICE_Larv  ) {
                  /*delay en el radio*/
                  delay_LP[j][i] = 1;
			     }
              else { if ( c > INDICE_Pup ) {
                  /*delay en el radio*/
                  delay_LP[j][i] = 1;}
                  }
              
                /*control con lag*/
                 if ( i > LAG ){
                 
                    if ( delay_LP[j][i-LAG] == 1){
                    
                        Larv[!ant][j] = (int)(EFIC_LP*Larv[!ant][j]);
                        Pupa[!ant][j] = (int)(EFIC_LP*Pupa[!ant][j]);
                        
                        inter_larv_pup++;
                        m = (i-LAG);
                        for ( pp = m ; pp <= i ; pp++){
                                delay_LP[j][pp] = 0;
                            }
                        
                     }
                 }
              
            }
            
             /*FIN CONTROL BIOLOGICO*/
              
             /*Tyranitar*/
              
             aux3 = poidev( dt_W*muM[j]*MOSCO[ant][j], &seed0);
             aux1 = poidev( dt_W*muMad[j]*MOSCO[ant][j], &seed0);
             
             
             MOSCO[!ant][j] = MOSCO[ant][j] + aux2 - aux3 - aux1;
             
             if (  MOSCO[!ant][j] < 0  ){
                MOSCO[!ant][j] = aux2;
                aux1 = (int)(0.5*MOSCO[ant][j]);
             }
             
             /* arriba el feminismo */
             aux1 = (int)(0.5*aux1);
             if ( aux1<0 ){ aux1 = 0.;}
             aux0 = aux1;
             
             /* muerte al macho*/
             if ( Tdeath < MATAR_VECTORES ){
			  
			  //EW[!ant][j] = (int)( SUV*EW[!ant][j] );
              aux3 =2.*dt_W*muV[j]*Vector[ant][j];
			  
			  }
              else {aux3 = dt_W*muV[j]*Vector[ant][j];}
            
             aux3 = poidev(aux3,&seed0);
             //aux0 = aux2; // para evitar los mosquitos jovenes
             
             /*Tyranitar Gx*/
             
             /* Vector[][] :  mosquito hembra */
             
             Vector[!ant][j] = Vector[ant][j] + aux0 - aux3;
             
             if ( Vector[!ant][j] < 0. ){
             
                   Vector[!ant][j] = aux1;
                   aux3 = Vector[ant][j];
                 
             }
             /* cantidad de Vectores muertos (aux4) */
             
             aux4 = aux3;
             
             /* reparto los muertos en los tres estados*/
             // num_ev -> numero de eventos
             // st_mn -> numero de estados a repartir multinomial
             // pr_mn -> prob de repartir multinomial
             
             num_ev = (int)(aux4);
             if (  num_ev < 0  ){ num_ev = 0;}
             
             st_mn = 2 + Ve_Lat; // : Vs,Ve1,Ve2,...,Ve3,Vi
             
             pr = 1./(float)(st_mn);
             
             //prob como vector para dist la multinomial
             for (  im = 0; im < st_mn -1; im++){ pr_mn[im] = pr; }
             
             // las mosquitas (si muerte al machop)
             VectorState = genmul(num_ev,pr_mn, st_mn);
             
			/* *******************************************************
			 * *******************************************************
			 * *******************************************************
			 * *******************************************************
			 */
              
		   /*campo medio con los vecinos censales
		    * el indice j es el correspondiente al numero
		    * de radio censal con el que se cuenta */
            
		    /* i es el dia, j es el radio censal*/
		    /*estoy agregando los infectados externos
		     * al termino de trasmision
		     * */
		     
		    MH          =   H[ant][j]  + inf_ext[i][j];
		    MV          =   Hi[ant][j] + inf_ext[i][j];
		    //if (inf_ext[i][j]==!0.){printf("inf: %d \t num: %.0f\t patch: %d\n",oo,inf_ext[i][j],j+1);oo++; }
		    for( hi = 0; hi < MAX_VECI; hi++ ){
				
				hf = cc[j][hi];
				
				if(hf ==! 0){
					
					MH+=H[ant][hf-1]  + inf_ext[i][hf-1];
                    
					MV+=Hi[ant][hf-1] + inf_ext[i][hf-1];
					
					}
					
				else{ break; }
			}
			//if (MV==!0.){printf("inf: %d \t num: %.0f\t patch: %d\n",oo,MV,j);oo++; }
			
			/*Agrego uno radio al azar que no puede ser de cc[j][hi] y 
			 * tiene que ser distinto de "j"
			 * parece que no funciona (si junio2019)*/
			hi = n_cont[j];
			kk = RAD_CENSAL - hi - 1;
			gg = (int)( kk*ran1(&seed));
			hf = sorteo[j][gg];
			MH+=H[ant][hf-1]  + inf_ext[i][hf-1];
			MV+=Hi[ant][hf-1] + inf_ext[i][hf-1];
			//if (MV==!0.){printf("inf: %d \t num: %.0f\t patch: %d\n",oo,MV,j);oo++; }
			/* *******************************************************
			 * *******************************************************
			 * *******************************************************
			 * *******************************************************
			 */
			  
            /* epidemiologia moscos */
             
           //vector expuesto
             
		   aux0 = dt_W*bite_rate*bv*Vs[ant][j]*MV/MH; //con temperatura theta_T(av_Temp)*
           aux0 = poidev(aux0,&seed);
           //if (aux0 ==!0.){printf("inf: %d \t num: %.0f\t patch: %d\n",oo,aux0,j);oo++; }
           //rate exposed vector (temperature) shape parameter k=10 or 7 or 3 change def_tartagal.c
           sv = rate_Ve(av_Temp,(float)Ve_Lat);
           //if (aux0 ==!0.){printf("inf: %d \t rate: %.2f\t patch: %d\n",oo,sv,j);oo++; }
           /*
             aux3 = dt_W*sv*Ve[ant][j][0];
             
             aux3 = poidev(aux3,&seed0);
             
             Ve[!ant][j][0] = Ve[ant][j][0] + aux0 - aux3 -(float)(VectorState[0]);
             if ( Ve[!ant][j][0] < 0.){
             
                Ve[!ant][j] = 0.;
                aux3 = Ve[ant][j][0] + aux0 - (float)(VectorState[0]);
                
                if ( aux3<0. ){ aux3 = 0.;}
             }
             */
           //linear triks Exposed Vector
             for ( st = 0; st < Ve_Lat; st++){
             
                    aux3 = dt_W*sv*Ve[ant][j][st];
                    aux3 = poidev(aux3,&seed0);
                    if (aux3 ==!0.){printf("inf: %d \t rate: %.2f\t patch: %d\n",week,sv,j);oo++; }
             
                    Ve[!ant][j][st] = Ve[ant][j][st] + aux0 - aux3 - (float)(VectorState[st]);
             
                    if ( Ve[!ant][j][st] < 0.){
             
                    Ve[!ant][j][st] = 0.;
                    aux3 = Ve[ant][j][st] + aux0 - (float)(VectorState[st]);
                
                    if ( aux3<0. ){ aux3 = 0.;}
                    
                    }
                 
                    aux0 = aux3;
             
             
             }
             
             // vector infectado
             
             Vi[!ant][j] = Vi[ant][j] + aux3 - (float)(VectorState[10]);
             
             if( Vi[!ant][j] < 0. ){
             
                Vi[!ant][j] = 0.;
                
             }
             
           //Vector suceptible
             
           Vs[!ant][j] = Vector[!ant][j] - Vi[!ant][j];
           //if ( Vs[!ant][j] < 0. ) { Vs[!ant][j] = 0.;}
           
           //resto todo los expuestos
            for ( st = 0; st < Ve_Lat; st++){
                
                Vs[!ant][j] = Vs[!ant][j] - Ve[!ant][j][st];
             
            }
             
            if ( Vs[!ant][j] < 0. ){
                Vs[!ant][j] = 0.;
            }
             
            //re calculo los vectores del prox periodo por si todos fueron 0
            Vector[!ant][j] = Vs[!ant][j] +Vi[!ant][j];
             
            //sumo a todo los expuestos
            for ( st = 0; st < Ve_Lat; st++){
                
            Vector[!ant][j] = Vector[!ant][j] + Ve[!ant][j][st];
             
            }
            
		   /*Hospedadores S-E-I-R*/
		   /* *******************************************************
			* *******************************************************
		    * *******************************************************
		    * *******************************************************
			*/
			  
		   MH = H[ant][j];
           MV = Vi[ant][j];
           
		    for(  hi = 0; hi < MAX_VECI; hi++  ){
            
				hf = cc[j][hi];
				if( hf==!0 ){
                
					MH+= H[ant][hf-1];
					MV+= Vi[ant][hf-1];
                    
					}
				else{ break; }
			}

			/*Agrego uno radio al azar no es ni j, ni un vecino 
			 * sorteo[i][gg] es un arreglo de sorteo
			 * */
			hi = n_cont[j];
			kk = RAD_CENSAL - hi - 1;
			gg = (int)( kk*ran1(&seed));
			hf = sorteo[j][gg];
			MH+= H[ant][hf-1];
			MV+= Vi[ant][hf-1];
            
            /*pica mosco */
		   //aux0=poidev(theta_T(av_Temp)*bite_rate*bh*MV*Hs[ant][j]/MH,&seed); // con temp theta_T(av_Temp)*
           
		   if ( Tdeath < NO_INFECCION ){
			  
			  aux0 = 0.;//if(Tdeath <10.){printf("infect %d\t Tmin = %.2f\n",oo,Tdeath);oo++;}
			  
			  }
              else { aux0 = bite_rate*bh*MV*Hs[ant][j]/MH;} //theta_T(av_Temp)*
              
           aux0 = poidev(aux0, &seed);
          // if (aux0 == 0.){printf("infect %d\t Tmin%.2f\n",oo,Tdeath);oo++;}
		   Hs[!ant][j] = Hs[ant][j] - aux0;
		   
		   if( Hs[!ant][j] < 0. ) {
			   
			   Hs[!ant][j] = 0.;
			   aux0 = Hs[ant][j];
               
			   }
			   
		   aux1 = poidev(sh*He[ant][j],&seed);
		   
		   He[!ant][j] = He[ant][j] + aux0 - aux1;
		   
		   if( He[!ant][j] < 0. ){
			   
               He[!ant][j] = aux0;
			   aux1 = He[ant][j];
			   
               }
		   
		   aux3 = poidev( dt_W*gh*Hi[ant][j], &seed );
		   
           /* inf_ext[i][j] es el infectado externo */
           
           aux2 = inf_ext[i][j];
           
		   Hi[!ant][j] = Hi[ant][j] + aux1 - aux3 + aux2;
           
		   if( Hi[!ant][j] < 0. ){
           
                Hi[!ant][j] = aux1 + aux2;
               
                }
              
           inc_y   =  inc_y + aux1; // contador de caso nuevo autoctono
           inc_yp  =  inc_yp + aux2; // contador de caso nuevo importado
           inc     =  inc + aux1 ;
           incimp  =  incimp + aux2;
           /*
           if ( (aux1 > 0.) || (aux2 > 1000.) ){
           oo++;
           printf ("inc=%.0f\timp=%.0f\tdia=%d\trad=%d\t%d\t%d\t%d\n",inc_y,inc_yp,i,j,corr,year_epid+2009,oo);
           }
           */
           
		   /*control biologico por reduccion de moscos*/
           
		   
		   if( inter_mosco < (LIMPIA) ){
              
                if ( (aux1 ==! 0) || (aux2 ==! 0) ){
   //printf ("MOSQUITO: %d\t%d\n",LIMPIA,inter_mosco);
				   MOSCO[!ant][j] = (int)(EFECTIV*MOSCO[!ant][j]);
				   Vs[!ant][j] = (int)(EFECTIV*Vs[!ant][j]);
				   Vi[!ant][j] = (int)(EFECTIV*Vi[!ant][j]);
                   Vector[!ant][j] = Vs[!ant][j] + Vi[!ant][j];
				   Larv[!ant][j] = (int)(EFECTIV*Larv[!ant][j]);
				   Pupa[!ant][j] = (int)(EFECTIV*Pupa[!ant][j]);
				   inter_mosco++;
                   for ( st = 0; st < Ve_Lat; st++){
                
                        Vector[!ant][j] = Vector[!ant][j] + (int)(EFECTIV*Ve[!ant][j][st]);
             
                        }
                   
                   }
				   
			   }
              
        //recuperados
           
      Hr[!ant][j] = H[!ant][j] - Hi[!ant][j] - He[!ant][j] - Hs[!ant][j];
		   
        //inc+=aux0;
	  }//j
		  
		   /*end of biological control fin de control biologico*/
		   
           /* INICIO DE LA DISPERSION*/
          
      /*    Dispersion de vectores suceptibles    */
      //reclutas a cero
      for ( j = 0; j < RAD_CENSAL; j++ ){ RVs[j] = 0.;}
      
      //Calculo los vectores suceptibles que se van VD[][] //menos al ultimo que tiene un solo contacto
      for ( j = 0; j < RAD_CENSAL; j++ ){
      
        VD[j] = bnldev(D_Vector,Vs[!ant][j],&seed);
        //calculo los Vs que se quedan
        aux0 = Vs[!ant][j] - VD[j];
        
        if (  aux0 < 0  ){
        
            VD[j]       =  Vs[!ant][j];
            Vs[!ant][j] =  0;
        }
        //Vs[][] que se quedan
        Vs[!ant][j] = aux0;
    
        num_ev = VD[j];
        
        st_mn = n_cont[j];
        
        pr = 1./(float)(st_mn);
        
        //prob como vector para dist la multinomial
        for (  im = 0; im < st_mn -1; im++){
          pr_mn[im] = pr;
        }
        //vectores repartidos a cada contacto
        VectorState = genmul(num_ev,pr_mn, st_mn);
        //Sumo los reclutas a cada radio
        for(  im = 0; im < st_mn ; im++  ){
        
            hf = cc[j][im];
            
            if( hf==!0 ){
            
					RVs[hf] = RVs[hf] + VectorState[hf-1];

                }
            else{ break;}
        }
      
      }
        
        
        //actualizo
        for ( j = 0; j < RAD_CENSAL; j++ ){
        
            Vs[!ant][j]      =  RVs[j]  +  Vs[!ant][j];
            Vector[!ant][j]  =  Vs[!ant][j]  +  Vi[!ant][j];
            
            for ( st = 0; st < Ve_Lat; st++){
            
                Vector[!ant][j] = Vector[!ant][j] + Ve[!ant][j][st];
            
            }
        
        }
      
      
      /* FIN DE LA DISPERSION*/
      /*Imprime a partir de 2009*/
      //if (i > 368 ){
      
      if (day == 7){
			  week++;
              day=0;
              //if ( year_epid == 4 ){
              fprintf(semanal,"%d\t%.0f\t%.0f\n",week,inc,incimp);
              fprintf(guardar,"%d\t",week);
          
              inc    =  0.;
              incimp =  0.;
          
              for ( j = 0; j < RAD_CENSAL; j++ ){
              
                        fprintf(guardar,"%.2f\t",Vs[ant][j]/(H[ant][j]));
                  
                }
          
              fprintf(guardar,"\n");
              //}
			}
		  //}
		  

          //{printf ("%d\t%.0f\t%.0f\t%.0f\t%.0f\t\n",2009 + year_epid,save_inc_ya[year_epid][corr],save_inc_ypc[year_epid][corr],inter_LP[year_epid][corr],inter_M[year_epid][corr]);}
       
	  ant=!ant;
      	  //printf("y = %d\n",year_epid);
	  //if (  YEARS < year_epid  ) {year_epid = 0 ;}
          
        
	  //printf(" sale i = %d, week = %d\n",i,week);
	  } //dia: i

  //corr++ ; printf("corr = %d\n",corr);
  //}
  
  //printf("\n");
  
  printf("LARVA: Oran_week_model\n");
  fclose(guardar);
  fclose(semanal);

	  clock_gettime(CLOCK_REALTIME, &timeStamp);
      printf("Tiempo de computo = %ld seg\n", 1000 * (timeStamp.tv_sec - oldTimeStamp.tv_sec) + (timeStamp.tv_nsec - oldTimeStamp.tv_nsec)/ 1000000000);
     
      
exit(0);
}
