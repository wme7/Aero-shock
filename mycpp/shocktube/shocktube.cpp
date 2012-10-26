// PHY 411-506 Computational Physics II Spring 2003
// shocktube.cpp

// Program to solve Sod's shock tube problem

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cstring>

using namespace std;

#include <GL/glut.h>
#include "RoeSolver.h"           // Roe's Riemann solver for Euler equations
#include "Riemann.h"             // Laney's upwind Godunov Riemann solver

double L = 1;                    // length of shock tube
double gama = 1.4;               // ratio of specific heats
int N = 200;                     // number of grid points

double **U;                      // solution with 3 components
double **newU;                   // new solution
double **F;                      // flux with 3 components
double *vol;                     // for Roe solver

double h;                        // lattice spacing
double tau;                      // time step
double CFL = 0.9;                // Courant-Friedrichs-Lewy number
int step;

void allocate() {
    static int oldN = 0;
    if (N != oldN) {
        if (U != 0) {
            for (int j = 0; j < oldN; j++) {
                delete [] U[j]; delete newU[j]; delete [] F[j];
            }
            delete [] U; delete [] newU; delete [] F; delete [] vol;
        }
        oldN = N;
        U = new double* [N];
        newU = new double* [N];
        F = new double* [N];
        vol = new double [N];
        for (int j = 0; j < N; j++) {
            U[j] = new double [3];
            newU[j] = new double [3];
            F[j] = new double [3];
        }
    }
}

double cMax() {
    double uMax = 0;
    for (int i = 0; i < N; i++) {
        if (U[i][0] == 0)
            continue;
        double rho = U[i][0];
        double u = U[i][1] / rho;
        double p = (U[i][2] - rho * u * u / 2) * (gama - 1);
        double c = sqrt(gama * abs(p) / rho);
        if (uMax < c + abs(u))
            uMax = c + abs(u);
    }
    return uMax;
}

void initialize() {

    allocate();
    h = L / (N - 1);
    for (int j = 0; j < N; j++) {
        double rho = 1, p = 1, u = 0;
        if (j > N / 2)
            rho = 0.125, p = 0.1;
        double e = p / (gama - 1) + rho * u * u / 2;
        U[j][0] = rho;
        U[j][1] = rho * u;
        U[j][2] = e;
        vol[j] = 1;
    }
    tau = CFL * h / cMax();
    step = 0;
}

void boundaryConditions(double **U) {

    // reflection boundary conditions at the tube ends
    U[0][0] = U[1][0];
    U[0][1] = -U[1][1];
    U[0][2] = U[1][2];
    U[N - 1][0] = U[N - 2][0];
    U[N - 1][1] = -U[N - 2][1];
    U[N - 1][2] = U[N - 2][2];
}

void LaxWendroffStep() {

    // compute flux F from U
    for (int j = 0; j < N; j++) {
        double rho = U[j][0];
        double m = U[j][1];
        double e = U[j][2];
        double p = (gama - 1) * (e - m * m / rho / 2);
        F[j][0] = m;
        F[j][1] = m * m / rho + p;
        F[j][2] = m / rho * (e + p);
    }    

    // half step
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            newU[j][i] = (U[j + 1][i] + U[j][i]) / 2 -
                         tau / 2 / h * (F[j + 1][i] - F[j][i]);
    boundaryConditions(newU);

    // compute flux at half steps
    for (int j = 0; j < N; j++) {
        double rho = newU[j][0];
        double m = newU[j][1];
        double e = newU[j][2];
        double p = (gama - 1) * (e - m * m / rho / 2);
        F[j][0] = m;
        F[j][1] = m * m / rho + p;
        F[j][2] = m / rho * (e + p);
    }    

    // step using half step flux
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            newU[j][i] = U[j][i] - tau / h * (F[j][i] - F[j - 1][i]);

    // update U from newU
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            U[j][i] = newU[j][i];
}

void LaxFriedrichsStep() {

    // compute flux F from U
    for (int j = 0; j < N; j++) {
        double rho = U[j][0];
        double m = U[j][1];
        double e = U[j][2];
        double p = (gama - 1) * (e - m * m / rho / 2);
        F[j][0] = m;
        F[j][1] = m * m / rho + p;
        F[j][2] = m / rho * (e + p);
    }

    // Lax-Friedrichs step
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            newU[j][i] = (U[j + 1][i] + U[j - 1][i]) / 2 -
                         tau / h * (F[j + 1][i] - F[j - 1][i]);
    boundaryConditions(newU);

    // update U from newU
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            U[j][i] = newU[j][i];
}

void upwindGodunovStep() {

    // find fluxes using Riemann solver
    for (int j = 0; j < N - 1; j++)
        Riemann(U[j], U[j + 1], F[j]);

    // update U
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            U[j][i] -= tau / h * (F[j][i] - F[j - 1][i]);
}

void RoeStep() {

    // compute fluxes at cell boundaries
    int icntl;
    RoeSolve(h, tau, gama, vol, U, F, N - 2, icntl);

    // update U
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            U[j][i] -= tau / h * (F[j + 1][i] - F[j][i]);
}

double nu = 0.0;

void LapidusViscosity() {

    // store Delta_U values in newU
    for (int j = 1; j < N; j++)
        for (int i = 0; i < 3; i++)
            newU[j][i] = U[j][i] - U[j - 1][i];

    // multiply Delta_U by |Delta_U|
    for (int j = 1; j < N; j++)
        for (int i = 0; i < 3; i++)
            newU[j][i] *= abs(newU[j][i]);

    // add artificial viscosity
    for (int j = 2; j < N; j++)
        for (int i = 0; i < 3; i++) 
            U[j][i] += nu * tau / h * (newU[j][i] - newU[j - 1][i]);
}

void (*stepAlgorithm)() = RoeStep; 
void redraw();

void takeStep() {
    boundaryConditions(U);
    tau = CFL * h / cMax();
    stepAlgorithm();
    LapidusViscosity();
    redraw();
    ++step;
}

int mainWindow, controlWindow, plotWindow[4];
int margin = 10, controlHeight = 30;
int buttons = 4;
int algorithm = 0;
char algorithmName[][20] = {"Roe Solver", "Lax Wendroff", 
                            "Upwind Godunov", "Lax Friedrichs"};
double yMin[] = {-1, -1, -0.2, -0.2};
double yMax[] = {2, 1, 3, 1.2};

void redraw() {
    for (int i = 0; i < 4; i++) {
        glutSetWindow(plotWindow[i]);
        glutPostRedisplay();
    }
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (glutGetWindow() == plotWindow[0])
        gluOrtho2D(0, 1, yMin[0], yMax[0]);
    else if (glutGetWindow() == plotWindow[1])
        gluOrtho2D(0, 1, yMin[1], yMax[1]);
    else if (glutGetWindow() == plotWindow[2])
        gluOrtho2D(0, 1, yMin[2], yMax[2]);
    else if (glutGetWindow() == plotWindow[3])
        gluOrtho2D(0, 1, yMin[3], yMax[3]);
    else
        gluOrtho2D(0, w, 0, h);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void display() {

    glClear(GL_COLOR_BUFFER_BIT);

    glColor3ub(0, 0, 0);
    glBegin(GL_LINES);
        glVertex2d(0, 0);
        glVertex2d(1, 0);
    glEnd();

    int plot = glutGetWindow();
    if (plot == plotWindow[0]) 
        glColor3ub(255, 0, 0);
    if (plot == plotWindow[1]) 
        glColor3ub(0, 255, 0);
    if (plot == plotWindow[2]) 
        glColor3ub(0, 0, 255);
    if (plot == plotWindow[3]) 
        glColor3ub(255, 0, 255);

    double avg = 0;
    glBegin(GL_LINE_STRIP);
        for (int j = 0; j < N; j++) {
            double y;
            if (plot == plotWindow[0]) 
                y = U[j][0];
            if (plot == plotWindow[1]) 
                y = U[j][1] / U[j][0];
            if (plot == plotWindow[2]) 
                y = U[j][2];
            if (plot == plotWindow[3]) 
                y = (U[j][2] - U[j][1] * U[j][1] / U[j][0] / 2) * (gama - 1);
            glVertex2d(j * h, y);
            avg += y;
        }
    glEnd();

    if (avg != 0.0) 
        avg /= N;
    for (int i = 0; i < 4; i++) {
        if (plot == plotWindow[i]) {
            glRasterPos2d(0.05, yMin[i] + 0.92 * (yMax[i] - yMin[i]));
            char plotName[][20] = {"Density", "Velocity", 
                                   "Energy", "Pressure"};
            char str[50];
            sprintf(str, "<%s> = %.4g", plotName[i], avg);
            for (int j = 0; j < strlen(str); j++)
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[j]);
        }
    }
    glColor3ub(0, 0, 0);
    glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {

    static bool running = true;

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        if (running) {
            glutIdleFunc(NULL);
            running = false;
        } else {
            glutIdleFunc(takeStep);
            running = true;
        }
        redraw();
    }
}

void mouseControl(int button, int state, int x, int y) {

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        int w = glutGet(GLUT_WINDOW_WIDTH);
        algorithm = buttons * x / w;
        switch (algorithm) {
        case 0:
            stepAlgorithm = RoeStep;
            initialize();
            break;
        case 1:
            stepAlgorithm = LaxWendroffStep;
            initialize();
            break;
        case 2:
            stepAlgorithm = upwindGodunovStep;
            initialize();
            break;
        case 3:
            stepAlgorithm = LaxFriedrichsStep;
            initialize();
            break;
        default:
            break;
        }
        glutPostRedisplay();
        redraw();
    }
}

void displayControl() {
    glClear(GL_COLOR_BUFFER_BIT);
    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    double dx = w / buttons;
    for (int b = 0; b < buttons; b++) {
        if (b == algorithm)
            glColor3ub(255, 0, 0);
        else
            glColor3ub(0, 255, 0);
        glRectd(b * dx, 0, (b + 1) * dx, h);
        glColor3ub(0, 0, 0);
        glRasterPos2d((b + 0.2) * dx, 0.3 * h);
        char *str = algorithmName[b];
        for (int j = 0; j < strlen(str); j++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[j]);
    }
    glColor3ub(0, 0, 255);
    double d = 0.1 * h;
    glRectd(0, 0, w, d);
    glRectd(0, h - d, w, h);
    for (int b = 0; b <= buttons; b++) {
        double x = b * dx - d / 2;
        if (b == 0) x += d / 2;
        if (b == buttons) x -= d / 2;
        glRectd(x, 0, x + d, h);
    }
    glutSwapBuffers();
}

void makeSubWindows() {

    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    int dx = (w - 3 * margin) / 2;
    int dy = (h - 4 * margin - controlHeight) / 2;
    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
        int x0 = margin * (1 + i) + i * dx;
        int y0 = margin * (1 + j) + j * dy;
        int n = 2 * i + j;
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
        plotWindow[n] = glutCreateSubWindow(mainWindow, x0, y0, dx, dy);
        glClearColor(1.0, 1.0, 0, 0);
        glShadeModel(GL_FLAT);
        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
        glutMouseFunc(mouse);
    }
    controlWindow = glutCreateSubWindow(mainWindow, margin, 
                                        h - margin - controlHeight,
                                        2 * dx + margin, controlHeight);
    glClearColor(0.0, 0.0, 1.0, 0.0);
    glutDisplayFunc(displayControl);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseControl);
}

void displayMain() {
    glClear(GL_COLOR_BUFFER_BIT);

    glutSwapBuffers();
}

void reshapeMain(int w, int h) {
    reshape(w, h);
    int dx = (w - 3 * margin) / 2;
    int dy = (h - 4 * margin - controlHeight) / 2;
    for (int i = 0; i < 2; i++) 
    for (int j = 0; j < 2; j++) {
        glutSetWindow(plotWindow[2 * i + j]);
        glutPositionWindow(margin * (1 + i) + i * dx,
                           margin * (1 + j) + j * dy);
        glutReshapeWindow(dx, dy);
    }
    glutSetWindow(controlWindow);
    glutPositionWindow(margin, h - margin - controlHeight);
    glutReshapeWindow(w - 2 * margin, controlHeight);
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    initialize();
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(100, 100);
    mainWindow = glutCreateWindow("Sod's shock tube problem");
    glutDisplayFunc(displayMain);
    glutReshapeFunc(reshapeMain);
    glutIdleFunc(takeStep);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
    makeSubWindows();
    glutMainLoop();
}