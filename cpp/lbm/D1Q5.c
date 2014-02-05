/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/* An implementation of the D1Q5 model for non-ideal systems */

#include <stdlib.h>
#include <math.h>

#include "mygraph.h"

#define xdim 100

double f0[xdim],f1[xdim],f2[xdim],f3[xdim],f4[xdim],n[xdim];

double n0=1, T0=0.33333333, Amp=0.01, omega=1;
double pc=1, theta=0.1, nc=1, kappa=0.1;
double ug[xdim],pg[xdim],tg[xdim];
int ugreq=0,pgreq=0,tgreq=0;
int next=0,pause=1,done=0,Repeat=1,iterations;

void init(){
  int i;
  iterations=0;
  for (i=0;i<xdim;i++){
    n[i]=n0+Amp*sin(2*M_PI*i/xdim);
    f0[i]=n[i]/4 *(4-5*T0+5*T0*T0);
    f1[i]=n[i]/6 *(  4*T0-5*T0*T0);
    f2[i]=n[i]/6 *(  4*T0-5*T0*T0);
    f3[i]=n[i]/24*(   -T0+5*T0*T0);
    f4[i]=n[i]/24*(   -T0+5*T0*T0);
  }
}

void init_analytical(){
  int i;
  iterations=0;
  for (i=0;i<xdim;i++){
    n[i]=n0+Amp*(1.0*random()/RAND_MAX-0.5);
    f0[i]=n[i]/4 *(4-5*T0+5*T0*T0);
    f1[i]=n[i]/6 *(  4*T0-5*T0*T0);
    f2[i]=n[i]/6 *(  4*T0-5*T0*T0);
    f3[i]=n[i]/24*(   -T0+5*T0*T0);
    f4[i]=n[i]/24*(   -T0+5*T0*T0);
  }
}

void (*iteration)();

void iteration_ideal(){
  double u,T,tmp,tmp2;
  int i;

  iterations++;
  for (i=0;i<xdim;i++){
    n[i]=f0[i]+f1[i]+f2[i]+f3[i]+f4[i];
    u=f1[i]-f2[i]+2*f3[i]-2*f4[i];
    T=(f1[i]+f2[i]+4*f3[i]+4*f4[i])/n[i]-u*u;
    f0[i]+=omega*(n[i]/4.*(4-5*T+5*T*T-(5-6*T)*u*u+u*u*u*u)-f0[i]);
    f1[i]+=omega*(n[i]/6.*(4*T-5*T*T+(4-3*T)*u+(4-6*T)*u*u-u*u*u-u*u*u*u)-f1[i]);
    f2[i]+=omega*(n[i]/6.*(4*T-5*T*T-(4-3*T)*u+(4-6*T)*u*u+u*u*u-u*u*u*u)-f2[i]);
    f3[i]+=omega*(n[i]/24.*(-T+5*T*T-(2-6*T)*u-(1-6*T)*u*u+2*u*u*u+u*u*u*u)-f3[i]);
    f4[i]+=omega*(n[i]/24.*(-T+5*T*T+(2-6*T)*u-(1-6*T)*u*u-2*u*u*u+u*u*u*u)-f4[i]);
  }
  tmp=f1[0];
  memmove(&f1[0],&f1[1],(xdim-1)*sizeof(double));
  f1[xdim-1]=tmp;
  tmp=f2[xdim-1];
  memmove(&f2[1],&f2[0],(xdim-1)*sizeof(double));
  f2[0]=tmp;
  tmp=f3[0];
  tmp2=f3[1];
  memmove(&f3[0],&f3[2],(xdim-2)*sizeof(double));
  f3[xdim-1]=tmp2;
  f3[xdim-2]=tmp;
  tmp=f4[xdim-1];
  tmp2=f4[xdim-2];
  memmove(&f4[2],&f4[0],(xdim-2)*sizeof(double));
  f4[0]=tmp2;
  f4[1]=tmp;
}

void iteration_forcing(){
  static double p[xdim],pp[xdim];
  double u,T,nu[xdim],TT[xdim],F,tmp,tmp2,phi,dn,ddn,dp,ddp;
  int i,ip,im,ipp,imm;

  iterations++;
  for (i=0;i<xdim;i++){
    n[i]=f0[i]+f1[i]+f2[i]+f3[i]+f4[i];
  }
  for (i=0;i<xdim;i++){
    ip=(i+1)%xdim;
    im=(i+xdim-1)%xdim;
    dn=0.5*(n[ip]-n[im]);
    ddn=(n[ip]-2*n[i]+n[im]);
    phi=(n[i]-nc)/nc;
    nu[i]=f1[i]-f2[i]+2*f3[i]-2*f4[i];
    TT[i]=0.333333333 /*(f1[i]+f2[i]+4*f3[i]+4*f4[i]-nu[i]*nu[i])/n[i]*/;
    p[i]=pc*(phi+1)*(phi+1)*(3*phi*phi-2*phi+1-2*theta)
      -kappa*(n[i]*ddn-0.5*dn*dn)-n[i]*TT[i];
  }
  for (i=0;i<xdim;i++){
    ip=(i+1)%xdim;
    im=(i+xdim-1)%xdim;
    pp[i]=1.5*p[i]-0.25*(p[ip]+p[im]);
  }
  for (i=0;i<xdim;i++){
    ip=(i+1)%xdim;ipp=(ip+1)%xdim;
    im=(i+xdim-1)%xdim; imm=(im+xdim-1)%xdim;
    T=TT[i];
    dp=(4.-3*T)/6.*(pp[ip]-pp[im])-(1-3*T)/12.*(pp[ipp]-pp[imm]);
    F=dp;
    u=nu[i]/n[i];
    f0[i]+=omega*(n[i]/4.*(4-5*T+5*T*T-(5-6*T)*u*u+u*u*u*u)-f0[i]);
    f1[i]+=omega*(n[i]/6.*(4*T-5*T*T+(4-3*T)*u+(4-6*T)*u*u-u*u*u-u*u*u*u)-f1[i])
      +(4-3*T)*F;
    f2[i]+=omega*(n[i]/6.*(4*T-5*T*T-(4-3*T)*u+(4-6*T)*u*u+u*u*u-u*u*u*u)-f2[i])
      -(4-3*T)*F;
    f3[i]+=omega*(n[i]/24.*(-T+5*T*T-(2-6*T)*u-(1-6*T)*u*u+2*u*u*u+u*u*u*u)-f3[i])
      -(2-6*T)*F;
    f4[i]+=omega*(n[i]/24.*(-T+5*T*T+(2-6*T)*u-(1-6*T)*u*u-2*u*u*u+u*u*u*u)-f4[i])
      +(2-6*T)*F;
  }
  tmp=f1[0];
  memmove(&f1[0],&f1[1],(xdim-1)*sizeof(double));
  f1[xdim-1]=tmp;
  tmp=f2[xdim-1];
  memmove(&f2[1],&f2[0],(xdim-1)*sizeof(double));
  f2[0]=tmp;
  tmp=f3[0];
  tmp2=f3[1];
  memmove(&f3[0],&f3[2],(xdim-2)*sizeof(double));
  f3[xdim-1]=tmp2;
  f3[xdim-2]=tmp;
  tmp=f4[xdim-1];
  tmp2=f4[xdim-2];
  memmove(&f4[2],&f4[0],(xdim-2)*sizeof(double));
  f4[0]=tmp2;
  f4[1]=tmp;
}

void SetIdeal(){
  iteration=iteration_ideal;
}

void SetForcing(){
  iteration=iteration_forcing;
}



void GUI(){
  static int xdimi=xdim;

  DefineGraphN_R("n",&n[0],&xdimi,NULL);
  DefineGraphN_R("u",&ug[0],&xdimi,&ugreq);
  DefineGraphN_R("p",&pg[0],&xdimi,&pgreq);
  DefineGraphN_R("T",&tg[0],&xdimi,&tgreq);
  StartMenu("D1Q3",1);
  DefineFunction("Set Ideal",&SetIdeal);
  DefineFunction("Set Forcing",&SetForcing);
  DefineInt("iterations",&iterations);
  DefineDouble("theta",&theta);
  DefineDouble("pc",&pc);
  DefineDouble("nc",&nc);
  DefineDouble("kappa",&kappa);
  DefineDouble("omega",&omega);
  DefineDouble("Amp",&Amp);
  DefineDouble("n0",&n0);
  DefineDouble("T0",&T0);
  DefineFunction("init",&init);
  DefineGraph(curve2d_,"Graphs");
  DefineInt("Repeat",&Repeat);
  DefineBool("next",&next);
  DefineBool("pause",&pause);
  DefineBool("done",&done);
  EndMenu();
}

void GetData(){
  int i;

  if (ugreq||tgreq) {
    for (i=0;i<xdim;i++) ug[i]=(f1[i]-f2[i]+2*f3[i]-2*f4[i])
			   /(f0[i]+f1[i]+f2[i]+f3[i]+f4[i]);
    ugreq=0;
  }
  if (pgreq) {
    for (i=0;i<xdim;i++) pg[i]=f1[i]+f2[i]+4*f3[i]+4*f4[i];
    pgreq=0;
  }
  if (tgreq) {
    for (i=0;i<xdim;i++) tg[i]=(f1[i]+f2[i]+4*f3[i]+4*f4[i])
			   /(f0[i]+f1[i]+f2[i]+f3[i]+f4[i])-ug[i]*ug[i];
    tgreq=0;
  }
}

int main(int argc, char *argv[]){
  int newdata=1;
  int i;

  SetIdeal();
  init();
  GUI();

  while (done==0){
    Events(newdata);
    GetData();
    DrawGraphs();
    if (next|| !pause){
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
