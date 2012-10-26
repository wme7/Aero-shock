// PHY 411-506 Computational Physics II Spring 2003
// burgers.cpp

// Program to solve the 1-D Burgers' Equation

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

#include <GL/glut.h>

const double pi = 4 * atan(1.0); // value of pi

double L = 1;                    // size of periodic region
int N = 200;                     // number of grid points
double h;                        // lattice spacing
double tau;                      // time step
double CFLRatio = 1;             // Courant-Friedrichs-Lewy ratio tau/h
enum {SINE, STEP};
int initialWaveform = SINE;      // sine function, step, etc.

double nu = 1e-6;                // kinematic viscosity
double *u;                       // the solution
double *uNew;                    // for updating
double *F;                       // the flow
double *uPlus, *uMinus;          // for Godunov scheme
int step;                        // integration step number

void allocate() {
    static int oldN = 0;
    if (oldN != N) {
        if (u != 0)
            delete [] u; delete [] uNew; delete [] F;
            delete [] uPlus; delete [] uMinus;
    }
    oldN = N;
    u = new double [N]; 
    uNew = new double [N];
    F = new double [N];
    uPlus = new double [N];
    uMinus = new double [N];
}

void initialize() {

    allocate();
    h = L / N;

    double uMax = 0;
    for (int i = 0; i < N; i++) {
        double x = i * h;
        switch (initialWaveform) {
        case SINE:
            u[i] = sin(2 * pi * x) + 0.5 * sin(pi * x);
        break;
        case STEP:
            u[i] = 0;
            if (x > L / 4 && x < 3 * L / 4)
                u[i] = 1;
            break;
        default:
            u[i] = 1;
            break;
        }
        if (abs(u[i]) > uMax)
            uMax = abs(u[i]);
    } 

    tau = CFLRatio * h / uMax;
    step = 0;
}

void (*integrationAlgorithm)();
void redraw();

void takeStep() {
    integrationAlgorithm();
    double *swap = u; 
    u = uNew;
    uNew = swap;
    redraw();
    ++step;
}

void FTCS() {
    for (int j = 0; j < N; j++) {
        int jNext = j < N - 1 ? j + 1 : 0;
        int jPrev = j > 0 ? j - 1 : N - 1;
        uNew[j] = u[j] * (1 - tau / (2 * h) * (u[jNext] - u[jPrev])) +
                  nu * tau / h / h * (u[jNext] + u[jPrev] - 2 * u[j]);
    }
}

void Lax() {
    for (int j = 0; j < N; j++) {
        int jNext = j < N - 1 ? j + 1 : 0;
        int jPrev = j > 0 ? j - 1 : N - 1;
        uNew[j] = (u[jNext] + u[jPrev]) / 2
                  - u[j] * tau / (2 * h) * (u[jNext] - u[jPrev])
                  + nu * tau / h / h * (u[jNext] + u[jPrev] - 2 * u[j]);
    }
}

void LaxWendroff() {
    for (int j = 0; j < N; j++)
        F[j] = u[j] * u[j] / 2;
    for (int j = 0; j < N; j++) {
        int jMinus1 = j > 0 ? j - 1 : N - 1;
        int jPlus1 = j < N - 1 ? j + 1 : 0;
        int jPlus2 = jPlus1 < N - 1 ? jPlus1 + 1 : 0;
        uNew[j] = (u[j] + u[jPlus1]) / 2 - 
                  (tau / 2 / h) * (F[jPlus1] - F[j]) + 
                  (nu * tau / (2 * h * h)) * (
                      (u[jPlus1] + u[jMinus1] - 2 * u[j]) / 2 +
                      (u[jPlus2] + u[j] - 2 * u[jPlus1]) / 2 );
    }
    for (int j = 0; j < N; j++)
        F[j] = uNew[j] * uNew[j] / 2;
    for (int j = 0; j < N; j++) {
        int jMinus1 = j > 0 ? j - 1 : N - 1;
        int jPlus1 = j < N - 1 ? j + 1 : 0;
        uNew[j] = u[j] - (tau / h) * (F[j] - F[jMinus1]) +
                  (nu * tau / (h * h)) * (u[jPlus1] + u[jMinus1] - 2 * u[j]);
    }
}

void Godunov() {

    for (int j = 0; j < N; j++) {
        uPlus[j] = u[j] > 0 ? u[j] : 0;
        uMinus[j] = u[j] < 0 ? u[j] : 0;
    }
    for (int j = 0; j < N; j++) {
        int jNext = j < N - 1 ? j + 1 : 0;
        int jPrev = j > 0 ? j - 1 : N - 1;
        double f1 = uPlus[jPrev] * uPlus[jPrev] / 2;
        double f2 = uMinus[j] * uMinus[j] / 2;
        F[jPrev] = f1 > f2 ? f1 : f2;
        f1 = uPlus[j] * uPlus[j] / 2;
        f2 = uMinus[jNext] * uMinus[jNext] / 2;
        F[j] = f1 > f2 ? f1 : f2;
        uNew[j] = u[j]  + nu * tau / h / h * (u[jNext] + u[jPrev] - 2 * u[j]);
        uNew[j] -= (tau / h) * (F[j] - F[jPrev]);
    }
}

int mainWindow, solutionWindow, controlWindow;
int margin = 10;
int controlHeight = 30;

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, 0, h);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void redraw() {
    glutSetWindow(solutionWindow);
    glutPostRedisplay();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    glutSwapBuffers();
}

void displaySolution() {
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3ub(255, 255, 255);
    glBegin(GL_LINE_STRIP);
        for (int i = 0; i < N; i++) {
            int iNext = i < N - 1 ? i + 1 : 0;
            glVertex2d(i * h, u[i]);
            glVertex2d((i + 1) * h, u[iNext]);
        }
    glEnd();
    char str[100];
    sprintf(str, "CFL Ratio = %.4f     nu = %.4g     t = %.4f", 
            CFLRatio, nu, step * tau);
    glRasterPos2d(0.02, -0.95);
    for (int j = 0; j < strlen(str); j++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[j]);
    glutSwapBuffers();
}

void (*method[])() = {FTCS, Lax, LaxWendroff, Godunov};
char methodName[][20] = {"FTCS", "Lax", "Lax Wendroff", "Godunov"};

void displayControl() {
    glClear(GL_COLOR_BUFFER_BIT);
    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    for (int i = 0; i < 4; i++) {
        if (method[i] == integrationAlgorithm)
            glColor3ub(255, 0, 0);
        else
            glColor3ub(0, 0, 255);
        glRectd((i + 0.025) * w / 4, 0.1 * h, (i + 0.975) * w / 4, 0.9 * h);
        glColor3ub(255, 255, 255);
        glRasterPos2d((i + 0.2) * w / 4, 0.3 * h);
        for (int j = 0; j < strlen(methodName[i]); j++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, methodName[i][j]);
    }
    glutSwapBuffers();
}

void reshapeMain(int w, int h) {
    reshape(w, h);

    glutSetWindow(solutionWindow);
    glutPositionWindow(margin, margin);
    glutReshapeWindow(w - 2 * margin, h - 3 * margin - controlHeight);

    glutSetWindow(controlWindow);
    glutPositionWindow(margin, h - margin - controlHeight);
    glutReshapeWindow(w - 2 * margin, controlHeight);
}

void reshapeSolution(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, 1, -1, +1.5);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void mouseSolution(int button, int state, int x, int y) {
    static bool running = false;

    switch (button) {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            if (running) {
                glutIdleFunc(NULL);
                running = false;
            } else {
                glutIdleFunc(takeStep);
                running = true;
            }
        }
        break;
    default:
        break;
    }
}

void mouseControl(int button, int state, int x, int y) {

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        int w = glutGet(GLUT_WINDOW_WIDTH);
        int algorithm = int(x / double(w) * 4);
        if (algorithm >= 0 && algorithm < 4)
            integrationAlgorithm = method[algorithm];
        glutPostRedisplay();
    }
}

void makeMainWindow() {
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(600, 400);
    glutInitWindowPosition(100, 100);
    mainWindow = glutCreateWindow("One-dimensional Burgers' Equation");
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glShadeModel(GL_FLAT);
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeMain);
}

void solutionMenu(int menuItem) {
    switch (menuItem) {
    case 1:
        initialWaveform = SINE;
        break;
    case 2:
        initialWaveform = STEP;
        break;
    default:
        break;
    }
    initialize();
    glutPostRedisplay();
}

void makeSolutionWindow() {
    glutSetWindow(mainWindow);
    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    solutionWindow = glutCreateSubWindow(mainWindow, margin, margin,
                     w - 2 * margin, h - 3 * margin - controlHeight);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
    glutDisplayFunc(displaySolution);
    glutReshapeFunc(reshapeSolution);
    glutMouseFunc(mouseSolution);
    integrationAlgorithm = Lax;
    glutCreateMenu(solutionMenu);
    glutAddMenuEntry("Initial Sine Waveform", 1);
    glutAddMenuEntry("Initial Step Waveform", 2);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void makeControlWindow() {
    glutSetWindow(mainWindow);
    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    controlWindow = glutCreateSubWindow(mainWindow, 
                      margin, h - margin - controlHeight,
                      w - 2 * margin, controlHeight);
    glClearColor(0.0, 1.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
    glutDisplayFunc(displayControl);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseControl);
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    if (argc > 1)
        CFLRatio = atof(argv[1]);
    if (argc > 2)
        nu = atof(argv[2]);
    initialize();
    makeMainWindow();
    makeSolutionWindow();
    makeControlWindow();
    glutMainLoop();
}
