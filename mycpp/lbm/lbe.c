/*
 * lbe.c Translation of Sauro Succi's fortran program
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* parameters */
#define nx 64
#define ny 32
#define npop 9

/* constants */
double cs2,cs22,cssq,fpois,den,visc,rt0,rt1,rt2;

/* logic */
int iforce = 1;			/* is there an applied force? */
int iobst = 0;			/* is there an obstacle? */

/* phys */
double omega = 1.9;		/* relaxation frequency */
double rho = 1.0;		/* initial density for Poiseuille force */
double u0 = 0;			/* initial velocity for Poiseuille force */
double uf = 0.1;		/* final velocity for Poiseuille force */
double fom;

/* arrays */
double u[nx+2][ny+2],
     v[nx+2][ny+2],
     feq[npop][nx+2][ny+2],
     f[npop][nx+2][ny+2];

/* count */
int istep;			/* step count */
int nsteps = 1000;		/* Number of steps */
int nout = 100;			/* steps between printing profile */
int ndiag = 100;		/* steps between performing diagnostics */
int nobst;			/* length of the obstacle (multiple of 2) */

/* files */
char fileout[6] = "RUN01";

/* functions */
void input(void), pbc(void), mbc(void), move(void),
     hydrovar(void), equili(void), collis(void),
     force(int, double*), obst(void), diag0D(void), profil(int, double);

int main (int argc, char *argv[]) {
     double frce;
     
     input();
     inithydro();
     equili();
     initpop();

     for (istep = 1; istep <= nsteps; istep++) {
	  pbc();
	  mbc();
	  move();
	  hydrovar();
	  equili();
	  collis();
	  if (iforce)
	       force(istep, &frce);
	  if (iobst)
	       obst();
	  if (istep % ndiag == 0)
	       diag0D();
	  if (istep % nout == 0)
	       profil(istep, frce);
     }

     return 0;
}

void input (void) {
     FILE *inf;
     char ans[100];
     int prompt = 0;

     printf("Enter data file name or press RETURN: ");
     fgets(ans, 100, stdin);
     if ((inf = fopen(ans, "r")) == 0) {
	  prompt = 1;
	  inf = stdin;
     }
     if (prompt) printf("Number of steps: ");
     fgets(ans, 100, inf);
     sscanf(ans, "%d", &nsteps);
     if (prompt) printf("Number of steps between printing profile: ");
     fgets(ans, 100, inf);
     sscanf(ans, "%d", &ndiag);
     if (prompt) printf("Relaxation frequency omega: ");
     fgets(ans, 100, inf);
     sscanf(ans, "%lf", &omega);
     if (prompt) printf("Applied force? 1=YES, 0=NO: ");
     fgets(ans, 100, inf);
     sscanf(ans, "%d", &iforce);
     if (prompt) printf("Initial density and velocity for Poiseuille force: ");
     fgets(ans, 100, inf);
     sscanf(ans, "%lf %lf", &rho, &u0);
     if (prompt) printf("Final velocity for the Poise force: ");
     fgets(ans, 100, inf);
     sscanf(ans, "%lf", &uf);
     if (prompt) printf("Linear obstacle? 1=YES, 0=NO: ");
     fgets(ans, 100, inf);
     sscanf(ans, "%d", &iobst);
     if (iobst) {
	  if (prompt) printf("Length of obstacle (multiple of 2): ");
	  fgets(ans, 100, inf);
	  sscanf(ans, "%d", &nobst);
     }
     if (prompt) printf("File for output (5 chars): ");
     fgets(fileout, 6, inf);
}
