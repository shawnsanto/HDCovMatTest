#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define MATHLIB_STANDALONE 1
#define xmalloc(nbytes) malloc_or_exit(nbytes, __FILE__, __LINE__)

void *malloc_or_exit(size_t nbytes, const char *file, int line)
{
	void *p=malloc(nbytes);
	if (p==NULL){
		fprintf(stderr, "%s, line %d:" "unable to allocate %lu bytes, calling exit()\n", file, line, ((unsigned long)nbytes));
	}
	return p;
}

double *dvector(int n)
{return xmalloc(n * sizeof(double));}

double trChat(double *Yvec, int n, int p, int T, int s, int s1, int h, int h1)
{
	int i,j,k,l,d;
	double temp,temp1,temp2,temp3,temp4,Yjh1,Yks1,Yis1,Ykh1,Yls1,Ys=0.0,Yh=0.0,Ys1=0.0,Yh1=0.0,trCov=0.0,trCov0=0.0,trCov1=0.0,trCov2=0.0,trCov3=0.0;

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			if (i!=j)
			{
				temp=0.0;
				temp1=0.0;
				for (d=0; d<p;d++)
				{
					Ys=Yvec[(unsigned long)(s-1)*p*n+i*p+d];
					Yh=Yvec[(unsigned long)(h-1)*p*n+j*p+d];
					Ys1=Yvec[(unsigned long)(s1-1)*p*n+i*p+d];
					Yh1=Yvec[(unsigned long)(h1-1)*p*n+j*p+d];
					temp+=Ys*Yh;
					temp1+=Ys1*Yh1;
				}
				trCov+=temp*temp1/(((double) n)*(n-1.0));
				for (k=0;k<n;k++)
                		{
                 		if (k!=i && k!=j)
                 		{
		                   temp2=0.0;
                		   temp3=0.0;
                   		for (d=0;d<p;d++)
                   		{
                    		Yjh1=Yvec[(unsigned long)(h1-1)*p*n+j*p+d];
                    		Yks1=Yvec[(unsigned long)(s1-1)*p*n+k*p+d];
                    		temp2+=Yjh1*Yks1;
                    		Yis1=Yvec[(unsigned long)(s1-1)*p*n+i*p+d];
                    		Ykh1=Yvec[(unsigned long)(h1-1)*p*n+k*p+d];
                    		temp3+=Yis1*Ykh1;
                   		}
                    		trCov1+=temp*temp2/((double) n*(n-1.0)*(n-2.0));
                    		trCov2+=temp*temp3/((double) n*(n-1.0)*(n-2.0));
                   		for (l=0;l<n;l++)
                   		{
                    		if (l!=i && l!=j && l!=k)
                    		{
                     		temp4=0.0;
                     		for (d=0;d<p;d++)
                     		{
                      		Ykh1=Yvec[(unsigned long)(h1-1)*p*n+k*p+d];
                      		Yls1=Yvec[(unsigned long)(s1-1)*p*n+l*p+d];
                      		temp4+=Ykh1*Yls1;
                      		}
                      		trCov3+=temp*temp4/((double) n*(n-1.0)*(n-2.0)*(n-3.0));
                     }
                    }
                  }
                }
            }
        }
    }

    trCov0=trCov-trCov1-trCov2+trCov3;
    return(trCov0);
}

double varTnk12(double *trCovvec, int n, int T, int k1, int k2)
{

    int s, s1, h, h1;
    double tr0,tr1,tr2,tr3,tr4,tr5,tr6,tr7,tr8,tr9,Varkvec = 0.0;

	for (s=1;s<(k1+1);s++)
	{
		for (s1=1; s1<(k2+1); s1++)
		{
		  tr0=trCovvec[(s-1)+T*(s1-1)+(unsigned long)T*T*(s-1)+(unsigned long)T*T*T*(s1-1)];

			for (h=(k1+1); h<(T+1); h++)
			{
			  tr2=trCovvec[(h-1)+T*(s1-1)+(unsigned long)T*T*(h-1)+(unsigned long)T*T*T*(s1-1)];
			  tr6=trCovvec[(s-1)+T*(s1-1)+(unsigned long)T*T*(h-1)+(unsigned long)T*T*T*(s1-1)];

				for (h1=(k2+1); h1<(T+1); h1++)
				{

                           tr1=trCovvec[(s-1)+T*(h1-1)+(unsigned long)T*T*(s-1)+(unsigned long)T*T*T*(h1-1)];

                		   tr3=trCovvec[(h-1)+T*(h1-1)+(unsigned long)T*T*(h-1)+(unsigned long)T*T*T*(h1-1)];

		                   tr4=trCovvec[(s-1)+T*(s1-1)+(unsigned long)T*T*(s-1)+(unsigned long)T*T*T*(h1-1)];

                		   tr5=trCovvec[(h-1)+T*(s1-1)+(unsigned long)T*T*(h-1)+(unsigned long)T*T*T*(h1-1)];

                		   tr7=trCovvec[(s-1)+T*(h1-1)+(unsigned long)T*T*(h-1)+(unsigned long)T*T*T*(h1-1)];

		                   tr8=trCovvec[(s-1)+T*(s1-1)+(unsigned long)T*T*(h-1)+(unsigned long)T*T*T*(h1-1)];

                		   tr9=trCovvec[(s-1)+T*(h1-1)+(unsigned long)T*T*(h-1)+(unsigned long)T*T*T*(s1-1)];

		                   Varkvec+=(tr0*tr0+tr1*tr1+tr2*tr2+tr3*tr3-2*tr4*tr4-2*tr5*tr5-2*tr6*tr6-2*tr7*tr7+2*tr8*tr8+2*tr9*tr9);
				}
			}
		}
	}

    Varkvec=4*Varkvec/((double) n*(n-1.0)*k1*k2*(T-k1)*(T-k2));

    return(Varkvec);
}


int Dkhat(double *Yvec,int n, int p, int T, double *Dk, int *khat)
{
	int k,j,t,s,s0,t0,k0,h,h1;
	double *trSigk=dvector(T),tempCk,tempCk1,tempCk2,Ck,Ck1,Ck2,mk,mTk,largest;

	for (k=0;k<T;k++)
	{
		j=k+1;
		trSigk[k]=trChat(Yvec, n, p, T, j, j, j, j);
	}

	for (k=0;k<1;k++)
	{
		Ck=0;
		k0=k+1;
		for (t=0;t<k+1;t++)
		{
			for (s=0;s<T-k0;s++)
			{
				s0=s+k0+1;
				t0=t+1;
				tempCk=trChat(Yvec, n, p, T, t0, t0, s0, s0);
				Ck+=tempCk;
			}
		}
		mk=0.0;
		for (h=0;h<k+1;h++)
		{
			mk+=trSigk[h]/(k+1.0);
		}
		mTk=0.0;
		for (h1=0;h1<T-k0;h1++)
		{
			mTk+=trSigk[h1+k0]/(T-k0);
		}
		Dk[k]=mk+mTk-2*Ck/(k0*(T-k0)+0.0);

	}

	for (k=1;k<T-1;k++)
    {
        Ck1=0;
        Ck2=0;
		k0=k+1;

        mk=0.0;
		for (h=0;h<k;h++)
		{
			mk+=trSigk[h];
			tempCk1=trChat(Yvec, n, p, T, k+1, k+1, h+1, h+1);
            		Ck1+=tempCk1;
		}
		mTk=0.0;
		for (h1=0;h1<T-k0;h1++)
		{
			mTk+=trSigk[h1+k0];
			tempCk2=trChat(Yvec, n, p, T, h1+k0+1, h1+k0+1, k+1, k+1);
            		Ck2+=tempCk2;
		}

        Dk[k]=((k0-1)*(T-(k0-1))*Dk[k-1]+(T-2*k0+1)*trSigk[k]-mk+mTk+2*Ck1-2*Ck2)/(k0*(T-k0));
    }

	khat[0]=1;
	largest=Dk[0];
	if (T>2)
        {
         for (k=0;k<T-2;k++)
            {
                k0=k+1;
                if (Dk[k0]>largest)
                {khat[0]=k0+1;
                largest=Dk[k0];}
            }
        }

	free(trSigk);
	return 0;
}


int testandCP(double *Yvec, int n, int p, int T, int *khat, double *stdDk, double *maxstdDk, double *CorrMat)
{
    int T0=T-1, k1, k2, nk2, s, s1, h, h1;
    double *trCovvec=dvector((unsigned long)T*T*T*T), *Dk=dvector(T-1), *CovMat=dvector(T0*(T0+1)/2), Varkvec, trvecval;

    Dkhat(Yvec,n, p, T, Dk, khat);

    for (s=1;s<(T+1);s++)
    {
            for (s1=1; s1<(T+1); s1++)
            {
                for (h=1; h<(T+1); h++)
                {
                    for (h1=1; h1<(T+1); h1++)
                    {
                        trvecval=trChat(Yvec, n, p, T, s, s1, h, h1);
                        trCovvec[(h1-1) + T*(h-1) + (unsigned long)T*T*(s1-1) + (unsigned long)T*T*T*(s-1)]=trvecval;
                    }
                }
            }
    }

    for (k1=0;k1<T0;k1++)
    {
        for (k2=0;k2<T0-k1;k2++)
        {
               nk2=k2+k1;
               Varkvec=varTnk12(trCovvec, n, T, k1+1, nk2+1);
               CovMat[(2*T0-k1+1)*k1/2+k2]=Varkvec;
        }
     }


    maxstdDk[0]=Dk[0]/sqrt(CovMat[0]);

    for (k1=0;k1<T0;k1++){
        stdDk[k1]=Dk[k1]/sqrt(CovMat[(2*T0-k1+1)*k1/2]);
        if (stdDk[k1]>maxstdDk[0])
            {maxstdDk[0]=stdDk[k1];}
    }

    for (k1=0;k1<T0;k1++)
    {
        for (k2=0;k2<T0-k1;k2++)
        {
            CorrMat[(2*T0-k1+1)*k1/2+k2]=CovMat[(2*T0-k1+1)*k1/2+k2]/sqrt(CovMat[(2*T0-k1+1)*k1/2]*CovMat[(2*T0-k1-k2+1)*(k1+k2)/2]);
        }
    }

    free(Dk);
    free(trCovvec);
    free(CovMat);
    return 0;
}

void testandCP4r(double *Yvec, int *n, int *p, int *T, int *khat, double *stdDk, double *maxstdDk, double *CorrMat)
{
  testandCP(Yvec, n[0], p[0], T[0], khat, stdDk, maxstdDk, CorrMat);
}


