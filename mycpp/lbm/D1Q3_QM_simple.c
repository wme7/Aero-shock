/* An implementation of the D1Q3 model for non-ideal systems */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mygraph.h"
int xdim,xdim_init=100;

double *f0=NULL,*f1=NULL,*f2=NULL,*n=NULL,*sn=NULL,*ninit=NULL,*p=NULL;

double n0=1, u0=0, Amp=0.00001, omega=1;
double kappa=0.1;

double *dp=NULL,*j=NULL;
int dpreq=0,jreq=0;
int next=0,Pause=1,done=0,Repeat=1,iterations,GalileanMeasure=0;


#define Pi 3.14159265358979323846264338327950288419716939937510582097

/* This is the derivative of the potential that you put into your code.
   This will be the main aspect of your manipulation of the code.
   If you want to include a delta-function, that is clearly not quite so 
   easily achieved in c.
   You can either put a delta function directly into the iteration routine
   using an if statment (as in if (i=22) F=0.005 else F=0) or you can 
   still use the dV(x) macro using the 
   (i==22?0.005:0)
   c-construct here. Two delta function would be given by
   (i==22?0.005:i==44?-0.005:0)*/

/*#define dV(x) (i==xdim/4?Amp:i==3*xdim/4?-Amp:0)*/ /* Square well */

#define dV(x) (Amp*sin(2*Pi*x/xdim)) /* cos potential */

#define P(sn,dsn,ddsn) (-kappa*(sn*ddsn-dsn*dsn))

void initmem(){
  /* This is a routine that initializes the memory for the different fields*/
  xdim=xdim_init;
  f0=realloc(f0,xdim*sizeof(double));
  f1=realloc(f1,xdim*sizeof(double));
  f2=realloc(f2,xdim*sizeof(double));

  n=realloc(n,xdim*sizeof(double));
  sn=realloc(sn,xdim*sizeof(double));
  ninit=realloc(ninit,xdim*sizeof(double));
  p=realloc(p,xdim*sizeof(double));

  dp=realloc(dp,xdim*sizeof(double));
  j=realloc(j,xdim*sizeof(double));
}

void iteration(){
  double a,A,F,tmp,dn,dsn,ddsn,dp,ddp;
  int i,ip,im;

  /* This routine is the actual heart of the simulation. Here we calculate
     the macroscopic quantities of the probability density n[i] and the 
     current, which appears as "a" in this algorithm. */

  iterations++;
  for (i=0;i<xdim;i++){
    n[i]=f0[i]+f1[i]+f2[i];
    sn[i]=sqrt(n[i]);
  }
  for (i=0;i<xdim;i++){
    /* Here we calculate the local Schroedinger pressure */
    ip=(i+1)%xdim;
    im=(i+xdim-1)%xdim;
    dn=0.5*(n[ip]-n[im]);
    dsn=0.5*(sn[ip]-sn[im]);
    ddsn=(sn[ip]-2*sn[i]+sn[im]);
    p[i]=P(sn[i],dsn,ddsn);

    /* Now we perform the collision of the densities. */
    a=f1[i]-f2[i];
    F=n[i]*dV(i);
    A=a*a/n[i]+p[i]+(1/omega-0.5)*(dn*(a+0.5*F)/n[i]);
    f0[i]+=omega*(n[i]-A-f0[i])-2*F*a/n[i];
    f1[i]+=omega*(0.5*( a+A)-f1[i]) +0.5*( F+2*F*a/n[i]);
    f2[i]+=omega*(0.5*(-a+A)-f2[i]) +0.5*(-F+2*F*a/n[i]);
  }

  /* This is the streaming step, assuming periodic boundary conditions,
     written in a computationally efficient manner. */
  tmp=f1[xdim-1];
  memmove(&f1[1],&f1[0],(xdim-1)*sizeof(double));
  f1[0]=tmp;
  tmp=f2[0];
  memmove(&f2[0],&f2[1],(xdim-1)*sizeof(double));
  f2[xdim-1]=tmp;
}

void init(){
  int i;
  double a,A;
  /* This routine allows you to initialize the simulation by inputting 
     an arbitrary density and current at time zero of the simulation.
     Here you could also put in your analytical solutions to make sure 
     that they are consistent with the solution of the algorithm. */
  initmem();
  iterations=0;
  for (i=0;i<xdim;i++){
    n[i]=ninit[i]=n0+Amp*(1.0*rand()/RAND_MAX-0.5); /* initial density */
    a=n[i]*u0; /* initial current */
    A=4*a*a/n[i]+n[i]/3.;
    f0[i]=n[i]-A;
    f1[i]=0.5*( a+A);
    f2[i]=0.5*(-a+A);
  }
}



void GetData(){
  int i,di,ip,im,ii;
  double N,dsn,ddsn;
  /* This routine allows you to calculate quantities to be displayed 
     in the graphics routine, even if you do not need them during each
     iteration. */
  if (dpreq) {
    for (i=0;i<xdim;i++) {
      ip=(i+1)%xdim;
      im=(i+xdim-1)%xdim;
      dp[i]=-0.5*(p[ip]-p[im])+n[i]*dV(i);
    }
    dpreq=0;
  }
  if (jreq) {
    for (i=0;i<xdim;i++) {
      j[i]=f1[i]-f2[i]+0.5*n[i]*dV(i);
    }
    jreq=0;
  }
}

void GUI(){
  /* This routine initializes the Grapical User Interface (GUI). First we 
     define the fields we want to visualize, and in the second part we 
     structure the menu where we can manipulate variables during the 
     simulation. */
  SetDefaultColor(1);
  SetDefaultShape(1);
  DefineGraphN_Rp("n",&n,&xdim,NULL);
  DefineGraphN_Rp("ninit",&ninit,&xdim,NULL);
  SetDefaultColor(4);
  SetDefaultShape(1);
  DefineGraphN_Rp("-dp+nF",&dp,&xdim,&dpreq);
  DefineGraphN_Rp("p",&p,&xdim,NULL);
  SetDefaultColor(2);
  SetDefaultShape(2);
  DefineGraphN_Rp("j",&j,&xdim,&jreq);
  StartMenu("D1Q3",1);
  DefineInt("iterations",&iterations);
  DefineDouble("kappa",&kappa);
  DefineDouble("omega",&omega);
  DefineDouble("Amp",&Amp);
  DefineDouble("n0",&n0);
  DefineDouble("u0",&u0);
  StartMenu("Initializations",0);
  DefineInt("Xdim init",&xdim_init);
  DefineFunction("init",&init);
  EndMenu();
  DefineGraph(curve2d_,"Graphs");
  DefineInt("Repeat",&Repeat);
  DefineBool("next",&next);
  DefineBool("pause",&Pause);
  DefineBool("done",&done);
  EndMenu();
}



int main(int argc, char *argv[]){
  int newdata=1;
  int i;
  /* This is the main loop of the program. The execution starts here. 
     This routine is responsible of executing the program and calling 
     the Graphical User Interface (GUI). When the variable "done" is changed
     to a non-zero value in the GUI the program ends. Otherwise the
     program loops through the iterations and displays the new graphics
     through the "DrawGraphs" routine. The GUI is also able to change
     the program execution by calling functions, e.g. the init function
     in the Initializations sub-menu. 
     The program execution can be accelerated by increasing the value of
     the "Repeat" variable, because less time is spent re-drawing the 
     graphics. */
  
  init();
  GUI();

  while (done==0){
    Events(newdata);
    GetData();
    DrawGraphs();
    if (next|| !Pause){
      newdata=1;
      next=0;
      for (i=0;i<Repeat;i++){
	iteration();
      }
    }
    else sleep(1);
  }

  return 0;
}
