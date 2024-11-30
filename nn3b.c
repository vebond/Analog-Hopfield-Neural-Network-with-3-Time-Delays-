/* Program which calculate Lyapunov exponent          */ 
/* Neural network with external force fe=e*sin(fe*t)  */   
   
#define TESTTIME   
   
#ifdef TESTTIME   
#include <time.h>   
#endif   
   
#include<stdio.h>   
/* #include<conio.h> */   
#include<math.h>   
#include<stdlib.h>   
#define N 12
/* #define M 524288 */
#define M 100000
#define L1 1000
#define L2 1000
#define L3 1000
#define A 10
#define ID 32
#define K1 100000
   
void main()   
{   
   double f(long i,double t,double y[N],double ytt[N],double a[N][N],double c[N],double d,double e[K1],long j) ; 
   long i,lag ;
   long j,k,i1,i2 ;
   long m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12 ; 
   double sy[1],fe ;   
   double yt[L3+1][N],ytt[N] ;   
   double y[N],yy[N],k1[N],k2[N],k3[N],k4[N] ;   
   double sqy,lsqy,ly[1] ;   
   double a[N][N],c[N] ;   
   double t,dt ; 
   double e[K1],d ; 
   double y1[1],y2[1],y3[1],y4[1],y5[1],y6[1] ;   
   double y7[1],y8[1],y9[1],y10[1],y11[1],y12[1] ;   
   
   double yet[L3+1][N],yett[N] ;   
   double ye[N],yye[N],ke1[N],ke2[N],ke3[N],ke4[N] ;   
   
   double y00[1],y01[1] ;   
   double kk,ll ;

   FILE *FA ;   
   FILE *FY ;
   FILE *FO ;   
   FILE *FOI ;
   FILE *FS ;
   
   FILE *FE ; 
   FILE *FEE ;
   FILE *FR ;   
   
   FILE *FO1 ;   
   FILE *FO2 ;   
   FILE *FO3 ;   
   FILE *FO4 ;   
   FILE *FO5 ;   
   FILE *FO6 ;   
   FILE *FO7 ;   
   FILE *FO8 ;   
   FILE *FO9 ;   
   FILE *FO10 ;
   FILE *FO11 ;
   FILE *FO12 ;
   FILE *FYTT ;
   
#ifdef TESTTIME   
 time_t timer;   
 FILE *log1;   
 log1=fopen("time.log","a");   
 time(&timer);   
 fprintf(log1,"Begin %s\n",asctime(localtime(&timer)));   
#endif   

   FA = fopen("amatrix","r") ;
   FEE = fopen("gauss","r") ;
   FY = fopen("inity","r") ;
   FO = fopen("ampl","wb") ;   
   FOI = fopen("ampl_ex","w") ;   
   FS = fopen("energy","wb") ;
   FO1 = fopen("chn1","w") ;   
   FO2 = fopen("chn2","w") ;   
   FO3 = fopen("chn3","w") ;   
   FO4 = fopen("chn4","w") ;   
   FO5 = fopen("chn5","w") ;   
   FO6 = fopen("chn6","w") ;   
   FO7 = fopen("chn7","w") ;   
   FO8 = fopen("chn8","w") ;   
   FO9 = fopen("chn9","w") ;   
   FO10 = fopen("chn10","w") ;
   FO11 = fopen("chn11","w") ;
   FO12 = fopen("chn12","w") ;  
   FYTT = fopen("ytt","w") ;

   m1=0 ;   
   m2=1 ;   
   m3=2 ;   
   m4=3 ;   
   m5=4 ;   
   m6=5 ; 
   m7=6 ; 
   m8=7 ; 
   m9=8 ; 
   m10=9 ;
   m11=10 ;
   m12=11 ; 
   t=0.0 ; 
   dt=0.01 ; 
   lag=0 ; 
   d=0.0 ; 

   for (i=0; i<N; i++)   { 
/*       y[i]=((rand() % 1000)/250.-2.0)*1.0e-100 ;  */
       y[i]=0.000 ;
       ytt[i]=0.000 ; 
       ye[i]=0.000 ; 
       yett[i]=0.000 ; }

   for (i=0; i<N/3; i++) {
       c[i]=2.0;}
   for (i=N/3; i<2*N/3; i++) {
       c[i]=2.0;}
   for (i=2*N/3; i<N; i++) {
       c[i]=2.0;}

   for (i=0; i<K1; i++)  {
   fscanf(FEE," %lf ",&e[i]) ;
                           }
/*   for (i=0; i<K1; i++)  {
   fprintf(FET,"%10.6lf\n",e[i]) ;
                           }  /
 
/*   for (i1=0; i1<N; i1++) { 
       for(i2=0; i2<N; i2++)  {  
           a[i1][i2]=(rand() % 1000)/250.-2.0 ; 
            } } 
  
   a[0][3]=1.0 ;                          */ 
   
/*   for (i1=0; i1<N; i1++) {   
       for (i2=0; i2<N; i2++)  {   
       fwrite(a,sizeof(double),1,FA) ; 
            }}                               */

     for (i1=0; i1<N; i1++) {
       for (i2=0; i2<N; i2++)  {
       fscanf(FA," %lf ",&a[i1][i2]) ;
/*       a[i1][i2] = d*a[i1][i2]+1.0 ;  */
                               }
/*              fprintf(FA," \n") ;          */
                             }
     
     for (i1=0; i1<N; i1++) {
       for (i2=0; i2<N; i2++)  {
       printf(" %6.3lf ",a[i1][i2]) ;
                                }
       printf(" \n") ;
                            }
                                                                               
/*       for (i=0; i<N; i++)  {
       fwrite(y,sizeof(double),1,FY) ; 
                           }                 */
                             
       for (i=0; i<N; i++)  {                    
       fscanf(FY," %le \n",&y[i]) ;
       		ye[i] = y[i] ;
                           }

       for (i=0; i<N; i++)  {
       printf(" %le \n",y[i]) ;
                             }

       for (i=0; i<N; i++)  {
       printf(" %le \n",ye[i]) ;
                             }

       for (i1=0; i1<L3+1; i1++) {
       for (i2=0; i2<N; i2++) {
            yt[i1][i2]=0.0 ;
            yet[i1][i2]=0.0 ;
                              } }
       for (i=0; i<N; i++)  {
            yt[0][i]=y[i] ;
                            }        
                            
/*   y[0]=0.000001 ;    */   
   
   for (i=0; i<M; i++)    {   
        t=t+dt ;   
   
   if (i%ID==0) {
      fprintf(FOI,"%12ld  %20.10le ",i,y[0]) ;
      fprintf(FOI,"  %20.10le ",y[1]) ;   
      fprintf(FOI,"  %20.10le ",y[2]) ;   
      fprintf(FOI,"  %20.10le ",y[3]) ; 
      fprintf(FOI,"  %20.10le ",y[4]) ;   
      fprintf(FOI,"  %20.10le ",y[5]) ;   
      fprintf(FOI,"  %20.10le ",y[6]) ;   
      fprintf(FOI,"  %20.10le ",y[7]) ;   
      fprintf(FOI,"  %20.10le ",y[8]) ;   
      fprintf(FOI,"  %20.10le ",y[9]) ;
      fprintf(FOI,"  %20.10le ",y[10]) ;
      fprintf(FOI,"  %20.10le \n",y[11]) ;
      
      fprintf(FYTT,"%12ld  %20.10le \n",i,ytt[0]) ;
                             }

/* unperturbed trajectory */
   
       for (j=0; j<N/3; j++)        {   
           k1[j]=f(i,t,y,ytt,a,c,d,e,j)*dt ;}

       for (j=N/3; j<2*N/3; j++)        {
           k1[j]=f(i,t,y,ytt,a,c,d,e,j)*dt ;}
           
       for (j=2*N/3; j<N; j++)        {
            k1[j]=f(i,t,y,ytt,a,c,d,e,j)*dt ;}

       for (j=0; j<N; j++)        {
           yy[j]=y[j]+k1[j]+sqrt(dt)*d*e[i] ;}

        for (j=0; j<N/3; j++)       {   
            k2[j]=f(i,t,yy,ytt,a,c,d,e,j)*dt ;}

        for (j=N/3; j<2*N/3; j++)       {
            k2[j]=f(i,t,yy,ytt,a,c,d,e,j)*dt ;}
            
        for (j=2*N/3; j<N; j++)       {
            k2[j]=f(i,t,yy,ytt,a,c,d,e,j)*dt ;}

/* neuron numbers 0-3*/

  for (j=0; j<N/3; j++)          {   
      y[j] += (k1[j]+k2[j])/2.0+sqrt(dt)*d*e[i] ;   
      if (i<L1) {   
          yt[i+1][j]=y[j] ;   
          ytt[j]=0.0 ;   
               }   
      else     { 
              ytt[j]=yt[0][j] ;   
          for (k=0; k<L1; k++)   
              yt[k][j]=yt[k+1][j] ;   
              yt[L1][j]=y[j] ;
          }   
                                }
                 
/* neuron numbers 4-7*/

  for (j=N/3; j<2*N/3; j++)          {
      y[j] += (k1[j]+k2[j])/2.0+sqrt(dt)*d*e[i] ;
      if (i<L2) {
          yt[i+1][j]=y[j] ;
          ytt[j]=0.0 ;
               }
      else     {
              ytt[j]=yt[0][j] ;
           for (k=0; k<L2; k++)
               yt[k][j]=yt[k+1][j] ;
               yt[L2][j]=y[j] ;
           }
                                }
      
/* neuron numbers 8-11*/ 

  for (j=2*N/3; j<N; j++)          {
      y[j] += (k1[j]+k2[j])/2.0+sqrt(dt)*d*e[i] ;
      if (i<L3) {
          yt[i+1][j]=y[j] ;
          ytt[j]=0.0 ;
               }
      else     {
              ytt[j]=yt[0][j] ;
           for (k=0; k<L3; k++)
               yt[k][j]=yt[k+1][j] ;
               yt[L3][j]=y[j] ;
           }
                                }

/* perturbed trajectory */
/* neuron numbers 0-4 */

       for (j=0; j<N/3; j++)        {   
           ke1[j]=f(i,t,ye,yett,a,c,d,e,j)*dt ;}
           
       for (j=N/3; j<2*N/3; j++)        {
           ke1[j]=f(i,t,ye,yett,a,c,d,e,j)*dt ;}
           
       for (j=2*N/3; j<N; j++)        {
           ke1[j]=f(i,t,ye,yett,a,c,d,e,j)*dt ;}
           
       for (j=0; j<N; j++)        {
           yye[j]=ye[j]+ke1[j]+sqrt(dt)*d*e[i] ;}
   
        for (j=0; j<N/3; j++)       {   
            ke2[j]=f(i,t,yye,yett,a,c,d,e,j)*dt ;}   

        for (j=N/3; j<2*N/3; j++)       {
            ke2[j]=f(i,t,yye,yett,a,c,d,e,j)*dt ;}
            
        for (j=2*N/3; j<N; j++)       {
            ke2[j]=f(i,t,yye,yett,a,c,d,e,j)*dt ;}
            
/* neuron numbers 0-3 */

  for (j=0; j<N/3; j++)          { 
      ye[j] += (ke1[j]+ke2[j])/2.0+sqrt(dt)*d*e[i] ;   
      if (i==2*lag) ye[j]=ye[j]*(1.0+0.00000000000001) ;   
      if (i<L1) {   
          yet[i+1][j]=ye[j] ;   
          yett[j]=0.0 ;   
               }   
      else     {   
              yett[j]=yet[0][j] ;   
          for (k=0; k<L1; k++)   
              yet[k][j]=yet[k+1][j] ;   
              yet[L1][j]=ye[j] ;   
          }   
                                }   

/* neuron numbers 4-7 */

  for (j=N/3; j<2*N/3; j++)          {
      ye[j] += (ke1[j]+ke2[j])/2.0+sqrt(dt)*d*e[i] ;
      if (i==2*lag) ye[j]=ye[j]*(1.0+0.00000000000001) ;
      if (i<L2) {
          yet[i+1][j]=ye[j] ;
          yett[j]=0.0 ;
               }
      else     {
              yett[j]=yet[0][j] ;
          for (k=0; k<L2; k++)
              yet[k][j]=yet[k+1][j] ;
              yet[L2][j]=ye[j] ;
          }
                                }
/* neuron numbers 8-11 */

  for (j=2*N/3; j<N; j++)          {
      ye[j] += (ke1[j]+ke2[j])/2.0+sqrt(dt)*d*e[i] ;
      if (i==2*lag) ye[j]=ye[j]*(1.0+0.00000000000001) ;
      if (i<L3) {
          yet[i+1][j]=ye[j] ;
          yett[j]=0.0 ;
               }
      else     {
              yett[j]=yet[0][j] ;
          for (k=0; k<L3; k++)
              yet[k][j]=yet[k+1][j] ;
              yet[L3][j]=ye[j] ;
          }
                                }

 if (i>=2*lag) 
    {   
    if ((i-2*lag)%ID==0)   
       {   
      lsqy=0.0 ;   
      sqy=0.0 ;   
      for (j=0; j<N; j++)   
          { 
          lsqy+=(y[j]-ye[j])*(y[j]-ye[j]) ;   
          sqy+=y[j]*y[j] ;   
          }   
      ly[0]=A*sqrt(lsqy) ;   
      sy[0]=A*sqrt(sqy) ;   
   
         fwrite(ly,sizeof(double),1,FO) ;   
         fwrite(sy,sizeof(double),1,FS) ;
   
       }   
    }
   
 if (i>=2*lag)   
    {   
    if ((i-2*lag)%ID==0)   
       { 
       y1[0]=y[m1];   
       y2[0]=y[m2];   
       y3[0]=y[m3];   
       y4[0]=y[m4];   
       y5[0]=y[m5];   
       y6[0]=y[m6];   
       y7[0]=y[m7]; 
       y8[0]=y[m8];   
       y9[0]=y[m9];   
       y10[0]=y[m10];
       y11[0]=y[m11];
       y12[0]=y[m12];
         fprintf(FO1,"%20.10le\n",y1[0]);
         fprintf(FO2,"%20.10le\n",y2[0]);
         fprintf(FO3,"%20.10le\n",y3[0]);
         fprintf(FO4,"%20.10le\n",y4[0]);
         fprintf(FO5,"%20.10le\n",y5[0]);
         fprintf(FO6,"%20.10le\n",y6[0]);
         fprintf(FO7,"%20.10le\n",y7[0]);
         fprintf(FO8,"%20.10le\n",y8[0]);
         fprintf(FO9,"%20.10le\n",y9[0]);
         fprintf(FO10,"%20.10le\n",y10[0]);
         fprintf(FO11,"%20.10le\n",y11[0]);
         fprintf(FO12,"%20.10le\n",y12[0]);
       }   
    }
 
   
if (i%1024==0) printf("i= %ld \n",i);   
   
    } /* End of cycle for M=8192  */

      fclose(FA) ;
      fclose(FEE) ;
      fclose(FY) ;
      fclose(FO) ; 
      fclose(FOI) ;
      fclose(FS) ; 
      fclose(FO1) ;   
      fclose(FO2) ;   
      fclose(FO3) ;   
      fclose(FO4) ;   
      fclose(FO5) ;   
      fclose(FO6) ;   
      fclose(FO7) ;   
      fclose(FO8) ;   
      fclose(FO9) ;   
      fclose(FO10) ;
      fclose(FO11) ;
      fclose(FO12) ;
      fclose(FYTT) ;

   FE = fopen("ampl","rb") ;   
   FR = fopen("lapres.dat","w") ; 
   
         fread(y00,sizeof(double),1,FE) ;   
   
      kk=fabs(y00[0]) ;   
   
      fprintf(FR,"y0= %20.15lf  ",y00[0]) ;   
      fprintf(FR,"kk= %20.15lf \n",kk) ; 
   
   for (i=1; i<M; i++)       {   
         fread(y01,sizeof(double),1,FE) ;   
   
   ll=fabs(y01[0])/kk ;   
     if (ll<0.00000000000001) ll=0.00000000000001 ;   
     ll=10.0*log(ll)/i ;   
   
/*  printf("i=%5ld   %20.15f  \n",i,ll) ; */
    if(i%25==0)   fprintf(FR,"%5ld   %20.15lf  \n",i,ll) ;   
   
      }
   
      fclose(FE) ;   
      fclose(FR) ;
 
#ifdef TESTTIME 
 time(&timer); 
 fprintf(log1,"End %s\n",asctime(localtime(&timer))); 
 fclose(log1); 
#endif 

} 
 
 double f(long i,double t,double y[N],double ytt[N],double a[N][N],double c[N],double d,double e[K1],long j) 
      { 
      long k ; 
      double fe,b,p,s ; 
/*      e=0.0 ; 
      fe=1000 ;   */
      b=1.0 ; 
      d=0.0 ;
      s=0.0 ; 
      p=0.0 ; 
      for (k=0; k<N; k++)      { 
          s += c[k]*a[j][k]*tanh(b*ytt[k]-p) ;  } 
/*      printf("i = %ld\n",i) ;  */
/*      return (s-y[j]+d*e[i]) ;   */
      return (s-y[j]+d*e[i]) ;

      } 
 
