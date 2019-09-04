#include <math.h>
#include "nrutil.h"
#include "nrutil.c"
#include "SpecialFunctions.h"
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
//parameter estimation
//lineal
void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a,
 float *b, float *siga, float *sigb, float *chi2, float *q)
{
	int i;
	float wt,t,sxoss,sx=0.,sy=0.,st2=0.,ss,sigdat;
	wt=0.;
	*b=0.0;
	if(mwt)
	{
		ss=0.0;
		for(i=1;i<=ndata;i++)
		{
			wt=1.0/sqrt(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	}
	else
	{
		for(i=1;i<=ndata;i++)
		{
			sx += x[i]*wt;
			sy += y[i]*wt; 
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt)
	{
		for(i=1;i<=ndata;i++)
		{
			t=(x[i]-sxoss)/sig[i];
			st2 +=t*t;
			*b += t*y[i];
		}
	}
	else
	{
		for(i=1;i<=ndata;i++)
		{
			t=x[i]-sxoss;
			st2 +=t*t;
			*b += t*y[i];
		}
	}
	*b /= st2; //solve for a,b,siga and sigb.
	*a = (sy - sx*(*b))/ss;
	*siga = sqrt( ( 1.0 + sx*sy/(ss*st2) )/ss );
	*sigb = sqrt( 1.0/st2 );
	*chi2 = 0.0;   // calculate chi square
	*q=1.0;
	*q=1.0;
	if (mwt==0)
	{
		for(i=1;i<=ndata;i++){*chi2 += sqrt(y[i]-(*a)-(*b)*x[i]);}
		sigdat = sqrt( (*chi2)/(ndata-2) );
		*siga *= sigdat;
		*sigb *= sigdat;
	}
	else
	{
		for (i=1;i<=ndata; i++)
		{*chi2 += sqrt(  (y[i]-(*a)-(*b)*x[i])/sig[i] );}
		if(ndata>2)
		{
			 *q= gammq(0.5*(ndata-2),0.5*(*chi2));
		}
	}
}
//non linear model darft!!
/* Now i've to add other void */
//#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
void covsrt(float **covar, int ma, int ia[], int mfit)
/* Expand in storage the covariance matrix covar, so as to take into 
account parameters that are being held fixed
 (For the latter, return zero covariances.)*/
{
    int i,j,k;
    float swap;
    for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--)
	{
		if (ia[j])
		{
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		} 
	}
}

void gaussj(float **a, int n, float **b, int m)
/*Linear equation solution by Gauss-Jordan elimination
a[1..n][1..n] is the input matrix. b[1..n][1..m] is input containing 
the m right-hand side vectors. On output, a is replaced by its matrix 
inverse, and b is replaced by the corresponding set of solution vectors.*/
{
	int *indxc,*indxr,*ipiv; int i,icol,irow,j,k,l,ll; float big,dum,pivinv,swap;
	//The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting.
	indxc=ivector(1,n); 
	indxr=ivector(1,n); 
	ipiv=ivector(1,n);
	//This is the main loop over the columns to be reduced.
	//This is the outer loop of the search for a pivot element.
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++)
	{
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++)
				{
					if (ipiv[k] == 0)
					{
						if (fabs(a[j][k]) >= big)
						{
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		//We now have the pivot element, so we interchange rows, if needed, 
		//to put the pivot element on the diagonal. The columns are not 
		//physically interchanged, only relabeled: indxc[i], the column 
		//of the ith pivot element, is the ith column that is reduced, while
		// indxr[i] is the row in which that pivot element was originally located.
		//If indxr[i] ̸= indxc[i] there is an implied column interchange. 
		//With this form of bookkeeping, the solution b’s will end up in the 
		//correct order, and the inverse matrix will be scrambled by columns.
		if (irow != icol) 
		{
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]);
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]);
		}
		indxr[i]=irow; 
		indxc[i]=icol;
		//We are now ready to divide the pivot row by the  
		//pivot element, located at irow and icol. 
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix"); 
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++) 	//Next, we reduce the rows...
		if (ll != icol)			//...except for the pivot one, of course.
		{
			dum=a[ll][icol];
			a[ll][icol]=0.0;
			for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
			for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
		}
	}
	//This is the end of the main loop over columns of the reduction. 
	//It only remains to unscram- ble the solution in view of the column interchanges.
	//We do this by interchanging pairs of columns in the reverse order that the permutation was built up.
	for (l=n;l>=1;l--) 
	{
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	} 
	//And we are done. 
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
// 

//Used by mrqmin to evaluate the linearized fitting matrix alpha, 
//and vector beta as in (15.5.8), and calculate χ2
void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int))
{
	int i,j,k,l,m,mfit=0;
	float ymod,wt,sig2i,dy,*dyda;
	dyda=vector(1,ma);
	for(j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for(j=1;j<=mfit;j++)	//initialize (symmetric) alpha, beta
	{
		for (k=1;k<=j;k++) alpha[j][k]=0.0; 
		beta[j]=0.0;
	}
	*chisq=0.0;
	for(i=1;i<=ndata;i++)	//summation loop over all data
	{
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]); 
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) 
		{
			if (ia[l])
			{
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;   // and find chi square	

	}
	
	for (j=2;j<=mfit;j++)			//fill in the symmetric side
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}

void mrqmin(float x[], float y[], float sig[], int ndata, float a[], 
int ia[], int ma, float **covar, float **alpha, float *chisq,
void (*funcs)(float, float [], float *, float [], int), float *alamda)
{
	//void covsrt();
	//void gaussj();
	//void mrqcof();
	
	int j,k,l;
	static int mfit;
	static float ochisq,*atry,*beta,*da,**oneda;
	if (*alamda < 0.0) //initialization 
	{
		atry=vector(1,ma); beta=vector(1,ma); da=vector(1,ma);
		for (mfit=0,j=1;j<=ma;j++) if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs); 
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=1;j<=mfit;j++)  //Alter linearized fitting matrix, by augmenting diagonal elements
	{
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda)); 
		oneda[j][1]=beta[j];
	}
	gaussj(covar,mfit,oneda,1); //Matrix solution. 
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0)  //once converged, evaluate covariance matrix 
	{
		covsrt(covar,ma,ia,mfit); 
		covsrt(alpha,ma,ia,mfit); //spread out alpha to its full size too
		free_matrix(oneda,1,mfit,1,1); 
		free_vector(da,1,ma); 
		free_vector(beta,1,ma); 
		free_vector(atry,1,ma); 
		return;
	}
	
	for (j=0,l=1;l<=ma;l++)  //Did the trial succeed?
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) // success, accept the new solution.
		{
			*alamda *= 0.1; ochisq=(*chisq);
			for (j=1;j<=mfit;j++) 
			{
				for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
				beta[j]=da[j];
			}
			for (l=1;l<=ma;l++) a[l]=atry[l];
		} 
	else //failure, increase alamda and return.
	{ 
		*alamda *= 10.0;
        *chisq=ochisq;
    }
}
