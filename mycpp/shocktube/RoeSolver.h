// PHY 411-506 Computational Physics II Spring 2003
// RoeSolver.h

// Translation of G. Mellema's Roe Solver {\tt roesol.f}

#ifndef ROESOLVER_H_INCLUDED
#define ROESOLVER_H_INCLUDED

#include <algorithm>

void RoeSolve(double dr,               // spatial step
              double dt,               // time step
              double gamma,            // adiabatic index
              double *vol,             // volume factor for 3-D problem
              double **state,          // (rho, rho*u, e)  -- input
              double **flux,           // flux at cell boundaries -- output
              int meshr,               // number of interior points
              int& icntl               // diagnostic -- bad if != 0
             ) 
{
    const double tiny = 1e-30;
    const double sbpar1 = 2.0;
    const double sbpar2 = 2.0;

    // allocate temporary arrays
    double **fludif = new double* [meshr+2];
    double *rsumr   = new double [meshr+2];
    double *utilde  = new double [meshr+2];
    double *htilde  = new double [meshr+2];
    double *absvt   = new double [meshr+2];
    double *uvdif   = new double [meshr+2];
    double *ssc     = new double [meshr+2];
    double *vsc     = new double [meshr+2];
    double **a      = new double* [meshr+2];
    double **ac1    = new double* [meshr+2];
    double **ac2    = new double* [meshr+2];
    double **w      = new double* [meshr+2];
    double **eiglam = new double* [meshr+2];
    double **sgn    = new double* [meshr+2];
    double **fluxc  = new double* [meshr+2];
    double **fluxl  = new double* [meshr+2];
    double **fluxr  = new double* [meshr+2];
    double *ptest   = new double [meshr+2];
    int **isb       = new int* [meshr+2];
    for (int i = 0; i < meshr + 2; i++) {
        fludif[i] = new double [3];
        a[i]      = new double [3];
        ac1[i]    = new double [3];
        ac2[i]    = new double [3];
        w[i]      = new double [4];
        eiglam[i] = new double [3];
        sgn[i]    = new double [3];
        fluxc[i]  = new double [3];
        fluxl[i]  = new double [3];
        fluxr[i]  = new double [3];
        isb[i]    = new int [3];
    }

    // initialize control variable to 0
    icntl = 0;

    // find parameter vector w
    for (int i = 0; i <= meshr + 1; i++) {
        w[i][0] = sqrt(vol[i] * state[i][0]);
        w[i][1] = w[i][0] * state[i][1] / state[i][0];
        w[i][3] = (gamma - 1) * (state[i][2] - 0.5 * state[i][1] 
                  * state[i][1] / state[i][0]);
        w[i][2] = w[i][0] * (state[i][2] + w[i][3]) / state[i][0];
    }

    // calculate the fluxes at the cell center
    for (int i = 0; i <= meshr + 1; i++) {
        fluxc[i][0] = w[i][0] * w[i][1];
        fluxc[i][1] = w[i][1] * w[i][1] + vol[i] * w[i][3];
        fluxc[i][2] = w[i][1] * w[i][2];
    }

    // calculate the fluxes at the cell walls 
    // assuming constant primitive variables
    for (int n = 0; n < 3; n++) {
        for (int i = 1; i <= meshr + 1; i++) {
            fluxl[i][n] = fluxc[i - 1][n];
            fluxr[i][n] = fluxc[i][n];
        }
    }

    // calculate the flux differences at the cell walls
    for (int n = 0; n < 3; n++)
        for (int i = 1; i <= meshr + 1; i++)
            fludif[i][n] = fluxr[i][n] - fluxl[i][n];

    // calculate the tilded state variables = mean values at the interfaces
    for (int i = 1; i <= meshr + 1; i++) {
        rsumr[i] = 1 / (w[i - 1][0] + w[i][0]);

        utilde[i] = (w[i - 1][1] + w[i][1]) * rsumr[i];
        htilde[i] = (w[i - 1][2] + w[i][2]) * rsumr[i];

        absvt[i] = 0.5 * utilde[i] * utilde[i];
        uvdif[i] = utilde[i] * fludif[i][1];

        ssc[i] = (gamma - 1) * (htilde[i] - absvt[i]);
        if (ssc[i] > 0.0)
            vsc[i] = sqrt(ssc[i]);
        else {
            vsc[i] = sqrt(abs(ssc[i]));
            ++icntl;
        }
    }

    // calculate the eigenvalues and projection coefficients for each 
    // eigenvector
    for (int i = 1; i <= meshr + 1; i++) {
        eiglam[i][0] = utilde[i] - vsc[i];
        eiglam[i][1] = utilde[i];
        eiglam[i][2] = utilde[i] + vsc[i];
        for (int n = 0; n < 3; n++)
            sgn[i][n] = eiglam[i][n] < 0.0 ? -1 : 1;
        a[i][0] = 0.5 * ((gamma - 1) * (absvt[i] * fludif[i][0] + fludif[i][2]
                  - uvdif[i]) - vsc[i] * (fludif[i][1] - utilde[i] 
                  * fludif[i][0])) / ssc[i];
        a[i][1] = (gamma - 1) * ((htilde[i] - 2 * absvt[i]) * fludif[i][0]
                  + uvdif[i] - fludif[i][2]) / ssc[i];
        a[i][2] = 0.5 * ((gamma - 1) * (absvt[i] * fludif[i][0] + fludif[i][2]
                  - uvdif[i]) + vsc[i] * (fludif[i][1] - utilde[i]
                  * fludif[i][0])) / ssc[i];
    }

    // divide the projection coefficients by the wave speeds
    // to evade expansion correction
    for (int n = 0; n < 3; n++)
        for (int i = 1; i <= meshr + 1; i++)
            a[i][n] /= eiglam[i][n] + tiny;

    // calculate the first order projection coefficients ac1
    for (int n = 0; n < 3; n++)
        for (int i = 1; i <= meshr + 1; i++)
            ac1[i][n] = - sgn[i][n] * a[i][n] * eiglam[i][n];

    // apply the 'superbee' flux correction to made 2nd order projection
    // coefficients ac2
    for (int n = 0; n < 3; n++) {
        ac2[1][n] = ac1[1][n];
        ac2[meshr + 1][n] = ac1[meshr + 1][n];
    }

    double dtdx = dt / dr;
    for (int n = 0; n < 3; n++) {
        for (int i = 2; i <= meshr; i++) {
            isb[i][n] = i - int(sgn[i][n]);
            ac2[i][n] = ac1[i][n] + eiglam[i][n] * 
                        ((max(0.0, min(sbpar1 * a[isb[i][n]][n], max(a[i][n],
                        min(a[isb[i][n]][n], sbpar2 * a[i][n])))) +
                        min(0.0, max(sbpar1 * a[isb[i][n]][n], min(a[i][n],
                        max(a[isb[i][n]][n], sbpar2 * a[i][n])))) ) *
                        (sgn[i][n] - dtdx * eiglam[i][n]));
        }
    }

    // calculate the final fluxes
    for (int i = 1; i <= meshr + 1; i++) {
        flux[i][0] = 0.5 * (fluxl[i][0] + fluxr[i][0] + ac2[i][0]
                     + ac2[i][1] + ac2[i][2]);
        flux[i][1] = 0.5 * (fluxl[i][1] + fluxr[i][1] + 
                     eiglam[i][0] * ac2[i][0] + eiglam[i][1] * ac2[i][1] +
                     eiglam[i][2] * ac2[i][2]);
        flux[i][2] = 0.5 * (fluxl[i][2] + fluxr[i][2] + 
                     (htilde[i] - utilde[i] * vsc[i]) * ac2[i][0] +
                     absvt[i] * ac2[i][1] +
                     (htilde[i] + utilde[i] * vsc[i]) * ac2[i][2]);
    }

    // calculate test variable for negative pressure check
    for (int i = 1; i <= meshr; i++) {
        ptest[i] = dr * vol[i] * state[i][1] + 
                   dt * (flux[i][1] - flux[i + 1][1]);
        ptest[i] = - ptest[i] * ptest[i] + 2 * (dr * vol[i] * state[i][0] +
                   dt * (flux[i][0] - flux[i + 1][0])) * (dr * vol[i] *
                   state[i][2] + dt * (flux[i][2] - flux[i + 1][2]));
    }

    // check for negative pressure/internal energy and set fluxes
    // left and  right to first order if detected
    for (int i = 1; i <= meshr; i++) {
        if (ptest[i] <= 0.0 || (dr * vol[i] * state[i][0] + dt * (flux[i][0] 
                                - flux[i + 1][0])) <= 0.0) {

            flux[i][0] = 0.5 * (fluxl[i][0] + fluxr[i][0] +
                ac1[i][0] + ac1[i][1] + ac1[i][2]);
            flux[i][1] = 0.5 * (fluxl[i][1] + fluxr[i][1] + 
                eiglam[i][0] * ac1[i][0] + eiglam[i][1] * ac1[i][1] + 
                eiglam[i][2] * ac1[i][2]);
            flux[i][2] = 0.5 * (fluxl[i][2] + fluxr[i][2] + 
                (htilde[i]-utilde[i] * vsc[i]) * ac1[i][0] + 
                absvt[i] * ac1[i][1] + 
                (htilde[i] + utilde[i] * vsc[i]) * ac1[i][2]);
            flux[i + 1][0] = 0.5 * (fluxl[i + 1][0] + fluxr[i + 1][0] + 
                 ac1[i + 1][0] + ac1[i + 1][1] + ac1[i + 1][2]);
            flux[i + 1][1] = 0.5 * (fluxl[i + 1][1] + fluxr[i + 1][1] + 
                 eiglam[i + 1][0] * ac1[i + 1][0] + eiglam[i + 1][1] * 
                 ac1[i + 1][1] + eiglam[i + 1][2] * ac1[i + 1][2]);
            flux[i + 1][2] = 0.5 * (fluxl[i + 1][2] + fluxr[i + 1][2] + 
                 (htilde[i + 1] - utilde[i + 1] * vsc[i + 1]) * ac1[i + 1][0] 
                 + absvt[i + 1] * ac1[i + 1][1] + 
                 (htilde[i + 1] + utilde[i + 1] * vsc[i + 1]) * ac1[i + 1][2]);

            // Check if it helped, set control variable if not

            ptest[i] = (dr * vol[i] * state[i][1] + 
                       dt * (flux[i][1] - flux[i + 1][1]));
            ptest[i] = 2.0 * (dr * vol[i] * state[i][0] 
                + dt * (flux[i][0]-flux[i + 1][0])) * (dr * vol[i] * 
                state[i][2] + dt * (flux[i][2] - flux[i + 1][2]))
                - ptest[i] * ptest[i];
            if (ptest[i] <= 0.0 || (dr * vol[i] * state[i][0] + 
                    dt * (flux[i][0] - flux[i + 1][0])) <= 0.0) 
                icntl = icntl + 1;
        }
    }                    

    // free temporary arrays
    for (int i = 0; i < meshr + 2; i++) {
        delete [] fludif[i];
        delete [] a[i];
        delete [] ac1[i];
        delete [] ac2[i];
        delete [] w[i];
        delete [] eiglam[i];
        delete [] sgn[i];
        delete [] fluxc[i];
        delete [] fluxl[i];
        delete [] fluxr[i];
        delete [] isb[i];
    }
    delete [] fludif;
    delete [] rsumr;
    delete [] utilde;
    delete [] htilde;
    delete [] absvt;
    delete [] uvdif;
    delete [] ssc;
    delete [] vsc;
    delete [] a;
    delete [] ac1;
    delete [] ac2;
    delete [] w;
    delete [] eiglam;
    delete [] sgn;
    delete [] fluxc;
    delete [] fluxl;
    delete [] fluxr;
    delete [] ptest;
    delete [] isb;
}

#endif /* ROESOLVER_H_INCLUDED */
