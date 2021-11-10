/* Implementing the non-vectorized ABM to Fokker Planck mapping for 1D implementation of CD Driving.
Control implemented through dynamic variation of selection coefficient protocols
s(t) and sCD(t)
Shamreen Iram, Hinczewski Biotheory Group
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

float Th(float r)
{
  float theta=(r>=0.) ? (1.) : (0.);
  return theta;
}

float sVarCD(double x)
{
  float s,ds,scd;
  float a,b,mu;
  a=0.02;b=2000.; mu=0.0025;
  s=(float)(a*tanh(x/b));
  ds=a/(b*0.05*pow(cosh(x/b),2.));
  scd=s+(ds/pow((pow((2.*mu-s),2.)+4.*mu*s),0.5));
  return scd;
}

float sVar(double x)
{
  float s;
  float a,b;
  a=0.02;b=2000.;
  s=(float)(a*tanh(x/b));
  return s;
}

float sVarL(double x)
{
  float s;
  float K,Q,B,t0,r;
  K=0.02; Q=100.;t0=700.;B=0.003;
  s=-0.00002466+(float)K/(1.0 + Q*exp(-B*(x-t0)));
  return s;
}


float sVarCDL(double x)
{
float s,ds,scd;
  float K,Q,B,t0,r,mu;
  K=0.02; Q=100.;t0=700.;B=0.003;mu=0.0025;
  s=-0.00002466+(float)K/(1.0 + Q*exp(-B*(x-t0)));
  ds=(float)B*K*Q*exp(-B*(x-t0))/pow((1.+Q*exp(-B*(x-t0))),2.0);
  scd=s+(ds/(0.05*pow((pow((2.*mu-s),2.)+4.*mu*s),0.5)));
  return scd;
}


float X0() //Using Box Muller transform for initial inverse sampling of p_eq(t=0)
{
  float U1,U2,Z,X;
  float mu=0.5,sigma=0.0720188;//parameter values obtained from Gaussian approx fit to p_eq (see notebook random.nb)


  U1= ((float)rand()/(RAND_MAX));
  U2= ((float)rand()/(RAND_MAX));

  Z= pow(-(2.*log(U1)),0.5)*cos(2.*3.14159*U2);
  X=Z*sigma + mu;
  return X;
}

int main()
{

	/* Defining simulation constants and variables
	   K ~ carrying capacity
       d ~ death rate (constant)
       b0~ maximum birth rate
       mu, nu ~mutation rates A->B, B->A
       s: 1+s is relative fitness of A over B
       Ai ~ pop.no of mutant A in i-th generation
       Bi ~ pop.no of mutant B in i-th gen
       Ni ~ total population no Ai+Bi in i-th gen (constrained by K)
bi ~ birth rate of mutant A in i-th gen (will determine Ai+1)

	*/




	//FILE *fp1;
	//FILE *fp2;
  FILE *fp3;
	//FILE *fp4;
	//fp1 = fopen("VxVar.csv","w+");
	//fp2 = fopen("DxVar.csv","w+");
  fp3 = fopen("popAsVarCD.csv","w+");
	//fp4 = fopen("popBVarCD.csv","w+");


	float tAi,tBi,Ai,Bi,K=10000.;
	float b0=2.,d=0.05,s,mu=0.0025,nu=0.0025,x;


  float bi,sma,sna,smb,snb,rm[3],rn[3];
	int i,j,m,n,pos,k,l;//i,m,n,l counters. pos variable stores dx in the correct bin
  //no.of bins to divide span of x 0->1





  int imax; //no.of generations to run for
  printf("Enter imax \n");
  scanf("%d", &imax);
  //printf("Enter s \n");
  //scanf("%f",&s);
        //printf("\nyour choice of imax= %d \n",imax);

  float **trajxA = (float **)malloc((imax+1) * sizeof(float *));
  for (i=0; i<(imax+1); i++)
      trajxA[i] = (float *)malloc(1000 * sizeof(float));

  //float *trajA=(float*)malloc((imax+1)*sizeof(float));// arrays storing population nos. of each mutant  in each generation */
  //float *trajB=(float*)malloc((imax+1)*sizeof(float));
  //double *trajA=(double*)malloc((imax+1)*sizeof(double));//xi=Ai/(Ai+Bi)


	srand(time(NULL));

  clock_t begin=clock();

	//trajA[0]=Ai; trajB[0]=Bi;
  //trajxA[0][0]=(double)(Ai/(Ai+Bi));
	k=0;

while(k<1000){
  i=0;
  x=X0();
  Ai=(float)floor(x*K);
  Bi=(float)floor((1.-x)*K);
	while(i<=imax){
        trajxA[i][k]=(double)(Ai/(Ai+Bi));
        s=sVarCDL((double)i);
        tAi=Ai;tBi=Bi;
        bi=( (Ai+Bi)<=K) ? (b0*(1.-(Ai+Bi)/K)) : (0.);//calculating birth rate of mut. A in i-th gen


        sma=0.;sna=0.;smb=0.;snb=0.;
        for ( m=1; m<=tAi; m++ ) {
                for(j=0;j<=2;j++){rm[j]= ((float)rand()/(RAND_MAX));};
                sma=(Th(rm[0]-d)*(1.+Th(bi-rm[1])*Th(rm[2]-mu)))+sma;
                smb=(Th(rm[0]-d)*Th(bi-rm[1])*Th(mu-rm[2]))+smb;
                };

        for ( n=1; n<=tBi; n++ ) {
                for(j=0;j<=2;j++){rn[j]= ((float)rand()/(RAND_MAX));};
                sna=(Th(rn[0]-d)*Th(bi/(1.+s)-rn[1])*Th(nu-rn[2]))+sna;
                snb=(Th(rn[0]-d)*(1.+Th(bi/(1.+s)-rn[1])*Th(rn[2]-nu)))+snb;
                };

        Ai=sna+sma;
        Bi=snb+smb;

        //printf("bi= %lf  Ai+1= %lf Bi+1= %lf \n",bi,Ai,Bi);

            //trajA[i][k] =Ai;
            //trajB[i] =Bi;
            //trajxA[i][k]=(double)(Ai/(Ai+Bi));
            //pos=(int)(floor(trajx[i-1]/intv));
            //dx[pos]+=trajx[i]-trajx[i-1];
            //dx2[pos]+=pow((trajx[i]-trajx[i-1]),2.);
            //count[pos]++;
            i++;
		};k++;
  };



		/* Validation of Fokker Planck mapping: Calculation of D(x) and v(x) on the fly*/
		 //for( l=0; l<=imax; l++){
		//fprintf(fp1,"%lf \n",trajA[l]);
                //fprintf(fp2,"%lf \n",trajB[l]);
                 //};

    //for( j=0; j<50; j++){
      //dx[j]=(count[j]>0) ?(dx[j]/(double)count[j]) : (0.);//calc v(x)= <dx>
      //dx2[j]=(count[j]>0) ? (((dx2[j]/(double)count[j])-pow(dx[j],2.0))/2.) : (0.);//calc D(x)=(<dx^2>-<dx>^2 )/2
      //printf("dx=%.5f , dx^2=%.8f, count=%.1f \n",dx[j],dx2[j],count[j]);
      //fprintf(fp1,"%.2f,%.5f\n",intv*((float)j+1.),dx[j]);
      //fprintf(fp2,"%.2f,%.8f\n",intv*((float)j+1),dx2[j]);

      //printf(fp2,"%lf \n",trajB[l]);
       //};


    for( j=0; j<=imax; j++){
      for(l=0;l<1000;l++){
      fprintf(fp3,"%.5f,",trajxA[j][l]);
      //fprintf(fp4,"%.5f\n",trajB[j]);
       };
       fprintf(fp3,"\n");
     };

     clock_t end=clock();
     double t=((double)(end-begin)/CLOCKS_PER_SEC)/60.;
     printf("run time=%lf mins \n",t);

     //fclose(fp1);
     //fclose(fp2);
     fclose(fp3);
     //fclose(fp4);
     free(trajxA);
     //free(trajA);
     //free(trajB);


    }
