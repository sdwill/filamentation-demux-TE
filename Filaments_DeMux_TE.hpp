#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <complex.h>

using namespace std;

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

/* Arrays */
double complex **Gz, **Mz, **Dz, **Ez;
double complex **Tx, **Bx, **Hx;
double complex **Ty, **By, **Hy;

double complex ***t_Gz, ***t_Mz, ***t_Dz, ***t_Ez;
double complex ***t_Tx, ***t_Bx, ***t_Hx;
double complex ***t_Ty, ***t_By, ***t_Hy;

double **alpha_x, **beta_x;
double **alpha_y, **beta_y;

double **plasmaFreq; 
double **nu;

/* Parallellization Parameters */
int numThreads;
int chunk;

/* Mesh Size */
int nx;
int ny;

/* Constants */
double eps_d; 			// Dielectric permittivity
double eps_0; 	// Permittivity of free space
double mu_0;		// Permeability of free space
double e_c; 	// Electron charge
double m_e; 	// Electron mass
double c; 		// Speed of light

/* Simulation Parameters */
double ff; 						// Fill fraction of filaments
double lambda; 			// Wavelength
double ddx;					// Mesh spatial step size
double dt; 					// Mesh temporal step size
double freqWave;				// Input frequency
double delta_f;
double delta_x;
double nu_0; 				// Collision Frequency
double plasmaFreq_0; // Plasma Frequency
double numSteps; 					// Number of iterations simulation should execute

/* Initial Parameters */
double n_0;
double beta;
double h;
double b;

/* UPML Parameters */
double sigma_max;
double kappa_max;
int upmlSize_x;
int upmlSize_y;
int mp;
double *sigma_x;
double *kappa_x;
double *sigma_y;
double *kappa_y;

/* Functions */
void initFieldArrays(int nx, int ny);
void delFieldArrays(int nx);
void initUPML(int nx, int upmlSize_x, int ny, int upmlSize_y, double kappa_max, double sigma_max, int mp);
void delUPML();
void initAB(int nx, int ny);
void initFreqArrays(int nx, int ny, int start_x, int end_x, int start_y, int end_y, double plasmaFreq_0, int nu_0, int chunk);
void initMesh(int nx, int ny, int start_x, int end_x, int start_y, int end_y, double plasmaFreq_0, double nu_0, int upmlSize_x, int upmlSize_y, double kappa_max, double sigma_max, int mp, int chunk);
void delMesh(int nx);
double n (int T);
void updatePlasmaFreq(int nx, int ny, int start_x, int end_x, int start_y, int end_y, double plasmaFreq_0, int chunk);
double C1(int i, int j);
double C2(int i, int j);
double C3(int i, int j);
double C4(int i, int j);
double C5(int i, int j);
double C6(int i, int j);
double C7(int i, int j);
double C8(int i, int j);
double C9(int i, int j);
double C10(int i, int j);
double C11(int i, int j);
double C12(int i, int j);
double C13(int i, int j);
double C14(int i, int j);
double C15(int i, int j);
double AZ1(int i, int j);
double AZ2(int i, int j);
double AZ3(int i, int j);
double AZ4(int i, int j);
double AZ5(int i, int j);
double AZ6(int i, int j);
double AZ7(int i, int j);
double AZ8(int i, int j);
double AZ9(int i, int j);
double AZ10(int i, int j);
double AZ11(int i, int j);
double AZ12(int i, int j);
double AZ13(int i, int j);
double AZ14(int i, int j);
double AZ15(int i, int j);
double AX1(int i, int j);
double AX2(int i, int j);
double AX3(int i, int j);
double AX4(int i, int j);
double AX5(int i, int j);
double AY1(int i, int j);
double AY2(int i, int j);
double AY3(int i, int j);
double AY4(int i, int j);
double AY5(int i, int j);
void inputSource(int T);
void updateElectricField(int nx, int ny, int T, int chunk);
void updateMagneticField(int nx, int ny, int T, int chunk);
void resetTempArrays(int nx, int ny, int chunk);
void updateTempArrays(int nx, int ny, int chunk);
void mainLoop(int numSteps, int nx, int ny, int chunk);
void writeFields(int nx, int ny, int n);
