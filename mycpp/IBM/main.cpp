#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

double delta2(float xs, float ys, float x, float y) {

  double answer = 0.0;
  double tmpx, tmpy;

  if(x < xs) tmpx = fabs(xs - (x + dxValue));
  if(x >= xs) tmpx = fabs(xs - (x - dxValue));
  if(y < ys) tmpy = fabs(ys - (y + dyValue));
  if(y >= ys) tmpy = fabs(ys - (y - dyValue));
  if(fabs(x - xs) >= dxValue) tmpx = 0;
  if(fabs(y - ys) >= dyValue) tmpy = 0;
  answer = tmpx * tmpy / (dxValue * dyValue);
  return(answer);
}

main(void) {
  float xx[XD][YD], yy[XD][YD], xp[XD][YD], yp[XD][YD];
  float xu[XD][YD], xv[XD][YD], yu[XD][YD], yv[XD][YD];
  float x[CD], y[CD];
  int uc[XD][YD], vc[XD][YD], pc[XD][YD], buc[CD], bvc[CD];

// grid
  FILE *fpg;
  int i, j, h;

  if((fpg = fopen("result/grid.xyz", "w")) == NULL) {
    printf("cant open grid\n");
  }
  fprintf(fpg, "%d %d\n", imax, jmax);
  for (j = 1; j <= jmax; j++) {
    for (i = 1; i <= imax; i++) {
      xx[i][j] = 0+(i-1)*dxValue;
      fprintf(fpg, "%f ", xx[i][j]);
    }
    fprintf(fpg, "\n");
  }
  for (j = 1; j <= jmax; j++) {
    for (i = 1; i <= imax; i++) {
      yy[i][j] = 0+(j-1)*dyValue;
      fprintf(fpg, "%f ", yy[i][j]);
    }
    fprintf(fpg, "\n");
  }
  fclose(fpg);

//Circle(IB) 
  FILE *fpc;
  double ang;

  ang=(pi*2)/crclnm;
  if((fpc = fopen("result/circle.xy", "w")) == NULL) {
    printf("cant open circle\n");
  }
  for(i=0; i<crclnm; i++) {
    x[i+1] = cx+rr*cos(i*ang);
    y[i+1] = cy+rr*sin(i*ang);
    fprintf(fpc, "%f  %f\n", x[i+1], y[i+1]);
  }
  fclose(fpc);

//U,V,P grid
  for(i = 1; i <= imax; i++) {
    for(j = 2; j <= jmax; j++) { 
      yu[i][j] = (yy[i][j - 1] + yy[i][j]) / 2;
      xu[i][j] = xx[i][j];
    }
  }
  for(i = 1; i <= imax; i++) {
    xu[i][1] = xx[i][1];
    xu[i][jmax + 1] = xx[i][jmax];
    yu[i][1] = 2 * yy[i][1] - yu[i][2];
    yu[i][jmax + 1] = 2 * yy[i][jmax] - yu[i][jmax];
  }

  for(j = 1; j <= jmax; j++) {
    for(i = 2; i <= imax; i++) {
      xv[i][j] = (xx[i - 1][j] + xx[i][j]) / 2;
      yv[i][j] = yy[i][j];
    }
  }
  for(j = 1; j <= jmax; j++) {
    xv[1][j] = 2 * xx[1][j] - xv[2][j];
    yv[1][j] = yy[1][j];
    xv[imax + 1][j] = 2 * xx[imax][j] - xv[imax][j];
    yv[imax + 1][j] = yy[imax][j];
  }
  for(j = 1; j <= (jmax + 1); j++) {
    for(i = 1; i <= (imax + 1); i++) {
      xp[i][j] = xv[i][j];
      yp[i][j] = yu[i][j];
    }
  }
  xp[imax + 1][jmax + 1] = xp[imax + 1][jmax];
  yp[imax + 1][jmax + 1] = yp[imax][jmax + 1];
  xp[1][jmax + 1] = xp[1][jmax];
  yp[1][jmax + 1] = yp[2][jmax + 1];
  xp[imax + 1][1] = xp[imax + 1][2];
  yp[imax + 1][1] = yp[imax][1];

//select inner points
  double length;
  for(j = 1; j <= jmax; j++) {
    for(i = 1; i <= imax; i++) {
      length = pow((xu[i][j] - cx), 2) + pow((yu[i][j] - cy), 2);
      if(sqrt(length) < rr) uc[i][j] = 1;
      else uc[i][j] = 0;

      length = pow((xv[i][j] - cx), 2) + pow((yv[i][j] - cy), 2);
      if(sqrt(length) < rr) vc[i][j] = 1;
      else vc[i][j] = 0;

      length = pow((xp[i][j] - cx), 2) + pow((yp[i][j] - cy), 2);
      if(sqrt(length) < rr) pc[i][j] = 1;
      else pc[i][j] = 0;
    }
  }

//find internal layer
  double disx, disy;

  for(h = 1; h <= crclnm; h++) {
    for(j = 1; j <= jmax; j++) {
      for(i = 1; i <= imax; i++) {
        disx = fabs(x[h] - xu[i][j]);
        disy = fabs(y[h] - yu[i][j]);
        if(disx < (0.1 * dxValue) && disy < (0.1 * dyValue)) {
          uc[i][j] = 2;
          buc[h] = 1;
        }

        disx = fabs(x[h] - xv[i][j]);
        disy = fabs(y[h] - yv[i][j]);
        if(disx < (0.1 * dxValue) && disy < (0.1 * dyValue)) {
          vc[i][j] = 2;
          bvc[h] = 1;
        }
      }
    }
  }

//coincide points
  double d;

  for(j = 1; j <= jmax; j++) {
    for(i = 1; i <= imax; i++) {
      if(uc[i][j] == 1) {
        for(h=1; h<=crclnm; h++) {
          d = delta2(x[h], y[h], xu[i][j], yu[i][j]);
          if(d>0.01 && buc[h]==0) uc[i][j] = 3;
        }
      }
      if(vc[i][j] == 1) {
        for(h=1; h<=crclnm; h++) {
          d = delta2(x[h], y[h], xv[i][j], yv[i][j]);
          if(d>0.01 && bvc[h]==0) vc[i][j] = 3;
        }
      }
    }
  }

//output flag

 FILE *fpu1, *fpu2, *fpu3, *fpv1, *fpv2, *fpv3, *fpp1;

 fpu1 = fopen("result/Uflag1", "w");
 fpu2 = fopen("result/Uflag2", "w");
 fpu3 = fopen("result/Uflag3", "w");
 fpv1 = fopen("result/Vflag1", "w");
 fpv2 = fopen("result/Vflag2", "w");
 fpv3 = fopen("result/Vflag3", "w");
 fpp1 = fopen("result/Pflag1", "w");
 for(j = 1; j <= jmax; j++) {
    for(i = 1; i <= imax; i++) {
      if (uc[i][j] == 1) fprintf(fpu1, "%f %f\n", xu[i][j], yu[i][j]);
      if (uc[i][j] == 2) fprintf(fpu2, "%f %f\n", xu[i][j], yu[i][j]);
      if (uc[i][j] == 3) fprintf(fpu3, "%f %f\n", xu[i][j], yu[i][j]);
      if (vc[i][j] == 1) fprintf(fpv1, "%f %f\n", xv[i][j], yv[i][j]);
      if (vc[i][j] == 2) fprintf(fpv2, "%f %f\n", xv[i][j], yv[i][j]);
      if (vc[i][j] == 3) fprintf(fpv3, "%f %f\n", xv[i][j], yv[i][j]);
      if (pc[i][j] == 1) fprintf(fpp1, "%f %f\n", xp[i][j], yp[i][j]);
    }
  }
 fclose(fpu1);
 fclose(fpu2);
 fclose(fpu3);
 fclose(fpv1);
 fclose(fpv2);
 fclose(fpv3);
 fclose(fpp1);
}



