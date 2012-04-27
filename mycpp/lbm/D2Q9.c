/* This is a standard two dimensional lattice Boltzmann program with
   periodic boundary conditions that should perform reasonably efficiently.
   It has been written for teaching purposes by Alexander Wagner
   at NDSU (March 2003). 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> /* for the sleep function. this may not be standard */
#include <math.h>
#include <mygraph.h>

#define xdim 31
#define ydim 21

#define DEBUG

double f0[xdim][ydim],
  f1[xdim][ydim],f2[xdim][ydim],f3[xdim][ydim],f4[xdim][ydim],
  f5[xdim][ydim],f6[xdim][ydim],f7[xdim][ydim],f8[xdim][ydim];

#ifdef DEBUG
double b1[xdim][ydim],b3[xdim][ydim],b5[xdim][ydim],b6[xdim][ydim],
  bb[xdim][ydim];
#endif

double omega=1;

/* Some constants that appear often and don't need to be calculated
   over and over again. So we calculate them just once here. */

double
  fourOnine=4.0/9.0,
  oneOnine=1.0/9.0,
  oneOthirtysix=1.0/36.0;

/* Some variables for the suspended particle */
typedef double *doublep;
typedef doublep pair[2];
pair *Bf1=NULL,*Bf3=NULL,*Bf5=NULL,*Bf6=NULL;
int Bf1c,Bf1m=0,Bf3c,Bf3m=0,Bf5c,Bf5m=0,Bf6c,Bf6m=0;
double Cx=15,Cy=10,Mx=0,My=0,Zxx=0,Zxy=0,Zyy=0,Ux=0,Uy=0,R=5,M=100;

int Repeat=1,iterations=0,FrequencyMeasure=100,Graphics=1;
int Pause=1,sstep=0,done=0;

double Amp=0.001,Amp2=0.1,Width=10;

/* Some special data types that are required to view the graphics */
int Xdim=xdim,Ydim=ydim,noreq;
double rho[xdim][ydim];
int rhoreq=0;
double u[xdim][ydim][2];
int ureq=0;

void BCmem(pair **p, int c, int *cm){
  if (c<*cm-2) return;
  *cm+=100;
  *p=(pair *) realloc(*p,*cm*sizeof(pair));
}

void circleBC(){
  /* Sets up the links that are cut by a circular object. */
  int x,y,Px,Py,Pxp,Pyp,R2;
  double dx,dy;
#ifdef DEBUG
  for (x=0;x<xdim;x++) for (y=0;y<ydim;y++)
    b1[x][y]=b3[x][y]=b5[x][y]=b6[x][y]=bb[x][y]=0;
#endif

  Bf1c=Bf3c=Bf5c=Bf6c=0;
  R2=R*R;
  for (x=Cx-R-2;x<Cx+R+2;x++){
    Px=(x+xdim)%xdim;
    Pxp=(x+1+xdim)%xdim;
    dx=x-Cx;
    for (y=Cy-R-2;y<Cy+R+2;y++){
      Py=(y+ydim)%ydim;
      Pyp=(y+1+ydim)%ydim;
      dy=y-Cy;
#ifdef DEBUG
      if ((dx*dx+dy*dy-R2)<=0) bb[Px][Py]=1;else bb[Px][Py]=-1;

#endif

      if ((dx*dx+dy*dy-R2)*((dx+1)*(dx+1)+dy*dy-R2)<=0){
	BCmem(&Bf1,Bf1c,&Bf1m);
	Bf1[Bf1c][0]=&(f2[Px][Py]);
	Bf1[Bf1c++][1]=&(f1[(Pxp)][Py]);
#ifdef DEBUG
	b1[Px][Py]=1;
      } else {
	b1[Px][Py]=-1;
#endif 
      }
      if ((dx*dx+dy*dy-R2)*((dx)*(dx)+(dy+1)*(dy+1)-R2)<=0){
	BCmem(&Bf3,Bf3c,&Bf3m);
	Bf3[Bf3c][0]=&(f4[Px][Py]);
	Bf3[Bf3c++][1]=&(f3[Px][Pyp]);
#ifdef DEBUG
	b3[Px][Py]=1;
      } else {
	b3[Px][Py]=-1;
#endif 
      }
      if ((dx*dx+dy*dy-R2)*((dx+1)*(dx+1)+(dy+1)*(dy+1)-R2)<=0){
	BCmem(&Bf5,Bf5c,&Bf5m);
	Bf5[Bf5c][0]=&(f7[Px][Py]);
	Bf5[Bf5c++][1]=&(f5[Pxp][Pyp]);
#ifdef DEBUG
	b5[Px][Py]=1;
      } else {
	b5[Px][Py]=-1;
#endif 
      }
      if (((dx)*(dx)+(dy+1)*(dy+1)-R2)*((dx+1)*(dx+1)+(dy)*(dy)-R2)<=0){
	BCmem(&Bf6,Bf6c,&Bf6m);
	Bf6[Bf6c][0]=&(f8[Pxp][Py]);
	Bf6[Bf6c++][1]=&(f6[Px][Pyp]);
#ifdef DEBUG
	b6[Px][Py]=1;
      } else {
	b6[Px][Py]=-1;
#endif 
      }
    }
  }
}

void initParticle(){
 Cx=10;Cy=10;Mx=0;My=0;Ux=0;Uy=0;R=5;M=100;
}

void init(){
  int x,y;
  double n,ux,uy,uxx,uyy,uxy,usq;

  printf("initialize\n");
  iterations=0;
  for (x=0;x<xdim;x++)
    for (y=0;y<ydim;y++){
      /* here we can define the macroscopic properties of the 
	 initial condition. */

      n=1+Amp2*exp(-(pow(x-xdim/2,2)+pow(y-ydim/2,2))/Width);
      ux=0/*.05*sin((y-ydim/2)*2*M_PI/xdim)*/;
      uy=0;
      /*n=1;ux=0;uy=0;*/

      /* The following code initializes the f to be the local equilibrium
	 values associated with the density and velocity defined above.*/

      uxx=ux*ux;
      uyy=uy*uy;
      uxy=2*ux*uy;
      usq=uxx+uyy;

      f0[x][y]=fourOnine*n*(1-1.5*usq);
      f1[x][y]=oneOnine*n*(1+3*ux+4.5*uxx-1.5*usq);
      f2[x][y]=oneOnine*n*(1-3*ux+4.5*uxx-1.5*usq);
      f3[x][y]=oneOnine*n*(1+3*uy+4.5*uyy-1.5*usq);
      f4[x][y]=oneOnine*n*(1-3*uy+4.5*uyy-1.5*usq);
      f5[x][y]=oneOthirtysix*n*(1+3*(ux+uy)+4.5*(uxx+uxy+uyy)-1.5*usq);
      f6[x][y]=oneOthirtysix*n*(1+3*(-ux+uy)+4.5*(uxx-uxy+uyy)-1.5*usq);
      f7[x][y]=oneOthirtysix*n*(1+3*(-ux-uy)+4.5*(uxx+uxy+uyy)-1.5*usq);
      f8[x][y]=oneOthirtysix*n*(1+3*(ux-uy)+4.5*(uxx-uxy+uyy)-1.5*usq);
      
    }
}

void TotMomentum(){
  int x,y;
  double Momx,Momy;

  Momx=Momy=0;
  for (x=0;x<xdim;x++)
    for (y=0;y<ydim;y++){
	Momx+=f1[x][y]-f2[x][y]+f5[x][y]-f6[x][y]-f7[x][y]+f8[x][y];
	Momy+=f3[x][y]-f4[x][y]+f5[x][y]+f6[x][y]-f7[x][y]-f8[x][y];
    }
  printf("MomF = (%e,%e), MomP = (%e,%e), MomT = (%e,%e)\n",
	 Momx,Momy,M*Ux,M*Uy,Momx+M*Ux,Momy+M*Uy);
}


void iterationColloid(){
#define alp 1  /* the degree of implicitness */
  double zxx,zxy,zyy; /* The inverse Z tensor */

  Cx+=Ux+xdim;
  Cx=fmod(Cx,xdim);
  Cy+=Uy+ydim;
  Cy=fmod(Cy,ydim);
  
  zxx=(M+alp*Zyy)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
  zxy=   alp*Zxy /((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
  zyy=(M+alp*Zxx)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);

  Ux+=zxx*(Mx-Zxx*Ux-Zxy*Uy)+zxy*(My-Zxy*Ux-Zyy*Uy);
  Uy+=zxy*(Mx-Zxx*Ux-Zxy*Uy)+zyy*(My-Zxy*Ux-Zyy*Uy);
#undef alp
}

void iteration(){
  int x,y,i;
  register double tmp1,tmp2;
  register  double n,ux,uy,uxx,uyy,uxy,usq,Fx,Fy,Fxx,Fyy,Fxy,Fsq;
  double f1y[ydim],f2y[ydim],f5y[ydim],f6y[ydim],f7y[ydim],f8y[ydim];
  double f3x[xdim],f4x[xdim],f5x[xdim],f6x[xdim],f7x[xdim],f8x[xdim];

  /* first we perform the collision step */
  
  for (x=0;x<xdim;x++)
    for (y=0;y<ydim;y++){
      n=f0[x][y]+f1[x][y]+f2[x][y]+f3[x][y]+f4[x][y]
	+f5[x][y]+f6[x][y]+f7[x][y]+f8[x][y];
      ux=f1[x][y]-f2[x][y]+f5[x][y]-f6[x][y]-f7[x][y]+f8[x][y];
      uy=f3[x][y]-f4[x][y]+f5[x][y]+f6[x][y]-f7[x][y]-f8[x][y];
      ux/=n;
      uy/=n;
      uxx=ux*ux;
      uyy=uy*uy;
      uxy=2*ux*uy;
      usq=uxx+uyy;
      /* We now need to implement any forcing terms. The term included
	 here is just an example. You should not calculate a constant
	 force term in the interation routine, but outside where it only
	 gets calculated once.*/
      Fx=Amp*sin(y*2*M_PI/ydim);
      Fy=0;
      Fxx=2*n*Fx*ux;
      Fyy=2*n*Fy*uy;
      Fxy=2*n*(Fx*uy+Fy*ux);
      Fsq=Fxx+Fyy;
      Fx*=n;
      Fy*=n;
      f0[x][y]+=omega*(fourOnine*n*(1-1.5*usq)-f0[x][y])
	-fourOnine*1.5*Fsq;
      f1[x][y]+=omega*(oneOnine*n*(1+3*ux+4.5*uxx -1.5*usq)-f1[x][y])
	+oneOnine*(3*Fx+4.5*Fxx-1.5*Fsq);
      f2[x][y]+=omega*(oneOnine*n*(1-3*ux+4.5*uxx -1.5*usq)-f2[x][y])
	+oneOnine*(-3*Fx+4.5*Fxx-1.5*Fsq);
      f3[x][y]+=omega*(oneOnine*n*(1+3*uy+4.5*uyy -1.5*usq)-f3[x][y])
	+oneOnine*(3*Fy+4.5*Fyy-1.5*Fsq);
      f4[x][y]+=omega*(oneOnine*n*(1-3*uy+4.5*uyy -1.5*usq)-f4[x][y])
	+oneOnine*(-3*Fy+4.5*Fyy-1.5*Fsq);
      f5[x][y]+=omega*(oneOthirtysix*n*(1+3*(ux+uy)+4.5*(uxx+uxy+uyy)
				      -1.5*usq)-f5[x][y])
	+oneOthirtysix*(3*(Fx+Fy)+4.5*(Fxx+Fxy+Fyy)-1.5*Fsq);
      f6[x][y]+=omega*(oneOthirtysix*n*(1+3*(-ux+uy)+4.5*(uxx-uxy+uyy)
				      -1.5*usq)-f6[x][y])
	+oneOthirtysix*(3*(-Fx+Fy)+4.5*(Fxx-Fxy+Fyy)-1.5*Fsq);
      f7[x][y]+=omega*(oneOthirtysix*n*(1+3*(-ux-uy)+4.5*(uxx+uxy+uyy)
				      -1.5*usq)-f7[x][y])
	+oneOthirtysix*(3*(-Fx-Fy)+4.5*(Fxx+Fxy+Fyy)-1.5*Fsq);
      f8[x][y]+=omega*(oneOthirtysix*n*(1+3*(ux-uy)+4.5*(uxx-uxy+uyy)
				      -1.5*usq)-f8[x][y])
	+oneOthirtysix*(3*(Fx-Fy)+4.5*(Fxx-Fxy+Fyy)-1.5*Fsq);
    }
  
  /* now we need to move the densities according to their velocities
     we are using periodic boundary conditions */

  /* since we are only using one lattice, we need to save some data on 
     the boundaries so that it does not get overwritten */

  for (y=0;y<ydim;y++){
    f1y[y]=f1[xdim-1][y];
    f2y[y]=f2[0][y];
    f5y[y]=f5[xdim-1][y];
    f6y[y]=f6[0][y];
    f7y[y]=f7[0][y];
    f8y[y]=f8[xdim-1][y];
  }
  for (x=0;x<xdim;x++){
    f3x[x]=f3[x][ydim-1];
    f4x[x]=f4[x][0];
    f5x[x]=f5[x][ydim-1];
    f6x[x]=f6[x][ydim-1];
    f7x[x]=f7[x][0];
    f8x[x]=f8[x][0];
  }
  
  /* Now we can move the densities along the lattice.
     You can also do this using loops, but that is actually
     more complicated */

  memmove(&f1[1][0],&f1[0][0],(xdim-1)*ydim*sizeof(double));
  memmove(&f2[0][0],&f2[1][0],(xdim-1)*ydim*sizeof(double));
  memmove(&f3[0][1],&f3[0][0],(xdim*ydim-1)*sizeof(double));
  memmove(&f4[0][0],&f4[0][1],(xdim*ydim-1)*sizeof(double));

  memmove(&f5[1][1],&f5[0][0],((xdim-1)*ydim-1)*sizeof(double));
  memmove(&f7[0][0],&f7[1][1],((xdim-1)*ydim-1)*sizeof(double));
  memmove(&f6[0][1],&f6[1][0],((xdim-1)*ydim)*sizeof(double));

  memmove(&f8[1][0],&f8[0][1],((xdim-1)*ydim)*sizeof(double));

  /* Now we need to fix the boundaries that have not yet been correctly
     updated */

  for (y=0;y<ydim;y++){
    f1[0][y]              =f1y[y];
    f2[xdim-1][y]         =f2y[y];
    f5[0][(y+1)%ydim]     =f5y[y];
    f6[xdim-1][(y+1)%ydim]=f6y[y];
    f7[xdim-1][y]         =f7y[(y+1)%ydim];
    f8[0][y]              =f8y[(y+1)%ydim];
  }
  for (x=0;x<xdim;x++){
    f3[x][0]              =f3x[x];
    f4[x][ydim-1]         =f4x[x];
    f5[(x+1)%xdim][0]     =f5x[x];
    f6[x][0]              =f6x[(x+1)%xdim];
    f7[x][ydim-1]         =f7x[(x+1)%xdim];
    f8[(x+1)%xdim][ydim-1]=f8x[x];
  }
  /* Objects in flow */
  Mx=My=Zxx=Zxy=Zyy=0;
  /* Really I am missing a factor of n here for the velocity corrections */
  for (i=0;i<Bf1c;i++){
    tmp1=*(Bf1[i][0]);
    tmp2=*(Bf1[i][1]);
    *(Bf1[i][0])=tmp2-2./3.*Ux;
    *(Bf1[i][1])=tmp1+2./3.*Ux;
    Mx+=2*(tmp2-tmp1);
    Zxx+= 2*2./3.;
  }
  for (i=0;i<Bf3c;i++){
    tmp1=*(Bf3[i][0]);
    tmp2=*(Bf3[i][1]);
    *(Bf3[i][0])=tmp2-2./3.*Uy;
    *(Bf3[i][1])=tmp1+2./3.*Uy;
    My+=2*(tmp2-tmp1);
    Zyy+= 2*2./3.;
  }
  for (i=0;i<Bf5c;i++){
    tmp1=*(Bf5[i][0]);
    tmp2=*(Bf5[i][1]);
    *(Bf5[i][0])=tmp2-1./6.*(Ux+Uy);
    *(Bf5[i][1])=tmp1+1./6.*(Ux+Uy);
    Mx+=2*(tmp2-tmp1);
    My+=2*(tmp2-tmp1);
    Zxx+= 2*1./6.;
    Zxy+= 2*1./6.;
    Zyy+= 2*1./6.;
  }
  for (i=0;i<Bf6c;i++){
    tmp1=*(Bf6[i][0]);
    tmp2=*(Bf6[i][1]);
    *(Bf6[i][0])=tmp2-1./6.*(-Ux+Uy);
    *(Bf6[i][1])=tmp1+1./6.*(-Ux+Uy);
    Mx+=-2*(tmp2-tmp1);
    My+= 2*(tmp2-tmp1);
    Zxx+= 2*1./6.;
    Zxy+=-2*1./6.;
    Zyy+= 2*1./6.;
  }
}

void analysis(int iterations){
  if (FrequencyMeasure%iterations==0){
    
  }
}

void GetGraphics(){
  double n;
  int x,y;

  if (rhoreq){
    rhoreq=0;
    for (x=0;x<xdim;x++)
      for (y=0;y<ydim;y++){
	rho[x][y]=f0[x][y]+f1[x][y]+f2[x][y]+f3[x][y]+f4[x][y]
	  +f5[x][y]+f6[x][y]+f7[x][y]+f8[x][y];
      }
  }
  if (ureq){
    ureq=0;
    for (x=0;x<xdim;x++)
      for (y=0;y<ydim;y++){
	n=f0[x][y]+f1[x][y]+f2[x][y]+f3[x][y]+f4[x][y]
	  +f5[x][y]+f6[x][y]+f7[x][y]+f8[x][y];
	u[x][y][0]=f1[x][y]-f2[x][y]+f5[x][y]-f6[x][y]-f7[x][y]+f8[x][y];
	u[x][y][1]=f3[x][y]-f4[x][y]+f5[x][y]+f6[x][y]-f7[x][y]-f8[x][y];
	u[x][y][0]/=n;
	u[x][y][1]/=n;
      }
  } 
}

void GUI(){
  DefineGraphNxN_R("rho",&rho[0][0],&Xdim,&Ydim,&rhoreq);
#ifdef DEBUG
  DefineGraphNxN_R("bb",&bb[0][0],&Xdim,&Ydim,&noreq);
  DefineGraphNxN_R("b1",&b1[0][0],&Xdim,&Ydim,&noreq);
  DefineGraphNxN_R("b3",&b3[0][0],&Xdim,&Ydim,&noreq);
  DefineGraphNxN_R("b5",&b5[0][0],&Xdim,&Ydim,&noreq);
  DefineGraphNxN_R("b6",&b6[0][0],&Xdim,&Ydim,&noreq);
#endif
  DefineGraphNxN_RxR("u",&u[0][0][0],&Xdim,&Ydim,&ureq);

  StartMenu("Lattice Boltzmann",1);
  DefineDouble("1/tau",&omega);
  DefineInt("Measurement freq.",&FrequencyMeasure);
  StartMenu("Reinitialize",0);
  DefineDouble("Amplitude",&Amp2);
  DefineDouble("Width",&Width);
  DefineFunction("reinitialize",&init);
  DefineFunction("Particle init",&initParticle);
  EndMenu();
  DefineDouble("Velocity amplitude",&Amp);
  StartMenu("Particle",0);
  DefineDouble("R",&R);
  DefineDouble("M",&M);
  DefineDouble("Ux",&Ux);
  DefineDouble("Uy",&Uy);
  DefineDouble("Cx",&Cx);
  DefineDouble("Cy",&Cy);
  DefineDouble("Mx",&Mx);
  DefineDouble("My",&My);
  EndMenu();
  DefineGraph(contour2d_,"Density&Vector plots");
  DefineInt("Repeat",&Repeat);
  DefineBool("Pause",&Pause);
  DefineBool("Single Step",&sstep);
  DefineBool("Done",&done);
  EndMenu();
}

int main(){
  int i,newdata;
  
  if (Graphics) GUI();

  init();
  while (!done){
    if (Graphics){
      Events(newdata);
      GetGraphics();
      DrawGraphs();
    } else {done=1;Pause=0;}
    if (!Pause||sstep){
      sstep=0;
      newdata=1;
      for (i=0;i<Repeat;i++){
	iterations++;
	circleBC();
	iteration();
	iterationColloid();
	TotMomentum();
	analysis(iterations);
      }
    } else sleep(1);
  }
  return 0;
}
