#include "Filaments_DeMux_TE.hpp"

void initFieldArrays(int nx, int ny)
{
/* Initializes the arrays for electric and magnetic fields */

/* Electric Components */
	Gz = new double complex *[nx];
	Mz = new double complex *[nx];
	Dz = new double complex *[nx];
	Ez = new double complex *[nx];
	
	
	for (int i = 0; i < nx; i++)
	{
		Gz[i] = new double complex [ny];
		Mz[i] = new double complex [ny];
		Dz[i] = new double complex [ny];
		Ez[i] = new double complex [ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Gz[i][j] = 0.0 + I*0.0;
			Mz[i][j] = 0.0 + I*0.0;
			Dz[i][j] = 0.0 + I*0.0;
			Ez[i][j] = 0.0 + I*0.0;
		}
	}
	
/* Magnetic Components */
	Tx = new double complex *[nx];
	Bx = new double complex *[nx];
	Hx = new double complex *[nx];
	
	for (int i = 0; i < nx; i++)
	{
		Tx[i] = new double complex [ny];
		Bx[i] = new double complex [ny];
		Hx[i] = new double complex [ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Tx[i][j] = 0.0;
			Bx[i][j] = 0.0;
			Hx[i][j] = 0.0;
		}
	}
	
	Ty = new double complex *[nx];
	By = new double complex *[nx];
	Hy = new double complex *[nx];
	
	for (int i = 0; i < nx; i++)
	{
		Ty[i] = new double complex [ny];
		By[i] = new double complex [ny];
		Hy[i] = new double complex [ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Ty[i][j] = 0.0;
			By[i][j] = 0.0;
			Hy[i][j] = 0.0;
		}
	}
	
/* Temporary arrays for finite difference approximations */
/* t_Gz */
	t_Gz = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Gz[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Gz[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Gz[i][j][0] = 0.0;
			t_Gz[i][j][1] = 0.0;
		}
	}
	
/* t_Mz */
	t_Mz = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Mz[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Mz[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Mz[i][j][0] = 0.0;
			t_Mz[i][j][1] = 0.0;
		}
	}

/* t_Dz */
	t_Dz = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Dz[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Dz[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Dz[i][j][0] = 0.0;
			t_Dz[i][j][1] = 0.0;
		}
	}

/* t_Ez */
	t_Ez = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Ez[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Ez[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Ez[i][j][0] = 0.0;
			t_Ez[i][j][1] = 0.0;
		}
	}

/* t_Tx */
	t_Tx = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Tx[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Tx[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Tx[i][j][0] = 0.0;
			t_Tx[i][j][1] = 0.0;
		}
	}

/* t_Bx */	
	t_Bx = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Bx[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Bx[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Bx[i][j][0] = 0.0;
			t_Bx[i][j][1] = 0.0;
		}
	}

/* t_Hx */	
	t_Hx = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Hx[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Hx[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Hx[i][j][0] = 0.0;
			t_Hx[i][j][1] = 0.0;
		}
	}

/* t_Ty */	
	t_Ty = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Ty[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Ty[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Ty[i][j][0] = 0.0;
			t_Ty[i][j][1] = 0.0;
		}
	}

/* t_By */	
	t_By = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_By[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_By[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_By[i][j][0] = 0.0;
			t_By[i][j][1] = 0.0;
		}
	}

/* t_Hy */	
	t_Hy = new double complex **[nx];
	
	for (int i = 0; i < nx; i++)
	{
		t_Hy[i] = new double complex *[ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Hy[i][j] = new double complex [2];
		}
	}
	
	for (int i  = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			t_Hy[i][j][0] = 0.0;
			t_Hy[i][j][1] = 0.0;
		}
	}
	
}

void delFieldArrays(int nx)
{
/* Deallocates memory space from field arrays */
	for (int i = 0; i < nx; i++)
	{
		delete[] Gz[i];
		delete[] Mz[i];
		delete[] Dz[i];
		delete[] Ez[i];
		
		delete[] Tx[i];
		delete[] Bx[i];
		delete[] Hx[i];
		
		delete[] Ty[i];
		delete[] By[i];
		delete[] Hy[i];
	}
	
	delete[] Gz;
	delete[] Mz;
	delete[] Dz;
	delete[] Ez;
	
	delete[] Tx;
	delete[] Bx;
	delete[] Hx;
	
	delete[] Ty;
	delete[] By;
	delete[] Hy;
}

void initUPML(int nx, int upmlSize_x, int ny, int upmlSize_y, double kappa_max, double sigma_max, int mp)
{
/* Initializes the UPML boundary conditions */
	sigma_x = new double [nx];
	kappa_x = new double [nx];
	sigma_y = new double [ny];
	kappa_y = new double [ny];
	
	for (int i = 0; i < nx; i++)
	{
		sigma_x[i] = 0.0; // Check these values!!!
		kappa_x[i] = 1.0; //
	}
	
	for (int j = 0; j < ny; j++)
	{
		sigma_y[j] = 0.0; //
		kappa_y[j] = 1.0; //
	}
	
// Boundary functions *************************************************
	for (int i = 0; i <= upmlSize_x; i++)
	{
			sigma_x[upmlSize_x - i] = double(sigma_max*pow(double((i*ddx)/(upmlSize_x*ddx)), mp));
			sigma_x[(nx - 1) - upmlSize_x + i] = double(sigma_max*pow(double((i*ddx)/(upmlSize_x*ddx)), mp));
	}
	
	for (int j = 0; j <= upmlSize_y; j++)
	{
			sigma_y[upmlSize_y - j] = double(sigma_max*pow(double((j*ddx)/(upmlSize_y*ddx)), mp));
			sigma_y[(ny - 1) - upmlSize_y + j] = double(sigma_max*pow(double((j*ddx)/(upmlSize_y*ddx)), mp));
	}
}

void delUPML()
{
	delete[] sigma_x;
	delete[] sigma_y;
	delete[] kappa_x;
	delete[] kappa_y;
}

void initAB(int nx, int ny)
{
/* Initializes the alpha and beta arrays */
	alpha_x = new double *[nx];
	alpha_y = new double *[nx];
	beta_x = new double *[nx];
	beta_y = new double *[nx];
	
	for (int i = 0; i < nx; i++)
	{
		alpha_x[i] = new double [ny];
		beta_x[i] = new double [ny];
		alpha_y[i] = new double [ny];
		beta_y[i] = new double [ny];
	}
}

void initFreqArrays(int nx, int ny, int start_x, int end_x, int start_y, int end_y, double plasmaFreq_0, double nu_0, int chunk)
{
/* Initializes plasma frequency array */
	plasmaFreq = new double *[nx];
	
	for (int i = 0; i < nx; i++)
	{
		plasmaFreq[i] = new double [ny];
	}

	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			plasmaFreq[i][j] = 0.0;
		} 
	}

	for (int j = start_y; j < end_y; j++)
	{
		for (int i = start_x; i < end_x; i++)
		{
			plasmaFreq[i][j] = plasmaFreq_0*(tanh((i - 50)/10) - tanh((i - 750)/10))/2.0;
		}
	}
	
/* Initializes collision frequency array */
	nu = new double *[nx];
	
	for (int i = 0; i < nx; i++)
	{
		nu[i] = new double [ny];
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			nu[i][j] = nu_0;
		} 
	}
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			alpha_x[i][j] = nu_0 + sigma_x[i]/kappa_x[i];
			alpha_y[i][j] = nu[i][j] + sigma_y[j]/kappa_y[j];
			beta_x[i][j] = pow(plasmaFreq[i][j], 2) + nu[i][j]*sigma_x[i]/kappa_x[i];
			beta_y[i][j] = pow(plasmaFreq[i][j], 2) + nu[i][j]*sigma_y[j]/kappa_y[j];
		}
	}
}

void initMesh(int nx, int ny, int start_x, int end_x, int start_y, int end_y, double plasmaFreq_0, double nu_0, int upmlSize_x, int upmlSize_y, double kappa_max, double sigma_max, int mp, int chunk)
{
/* Initializes the simulation mesh */
	initFieldArrays(nx, ny);
	initUPML(nx, upmlSize_x, ny, upmlSize_y, kappa_max, sigma_max, mp);
	initAB(nx, ny);
	initFreqArrays(nx, ny, start_x, end_x, start_y, end_y, plasmaFreq_0, nu_0, chunk);
}

void delMesh(int nx)
{
/* Deallocates the memory space for the mesh arrays */
	delFieldArrays(nx);
	delUPML();
	
	for (int i = 0; i < nx; i++)
	{
		delete[] plasmaFreq[i];
		delete[] nu[i];
		delete[] alpha_x[i];
		delete[] alpha_y[i];
		delete[] beta_x[i];
		delete[] beta_y[i];
	}
	
	delete[] plasmaFreq;
	delete[] nu;
	delete[] alpha_x;
	delete[] alpha_y;
	delete[] beta_x;
	delete[] beta_y;
}

double n (int T)
{
/* Free electron density function */
	return n_0*exp(-1.0*h*double(T)*dt)/(1.0 + b*double(T)*dt);
}

void updatePlasmaFreq(int nx, int ny, int start_x, int end_x, int start_y, int end_y, double plasmaFreq_0, int chunk)
{
/* Updates the plasma frequency in the mesh as the simulation iterates */
	#pragma omp parallel for schedule (static, chunk)
	for (int j = start_y; j < end_y; j++)
	{
		for (int i = start_x; i < end_x; i++)
		{
			plasmaFreq[i][j] = plasmaFreq_0*(tanh((i - 50)/10) - tanh((i - 750)/10))/2.0;
		}
	}
}
double C1(int i, int j)
{
	return (kappa_x[i]*(-1 + nu[i][j]*dt - 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) - sigma_x[i]*(dt - 0.5*nu[i][j]*pow(dt, 2)));
}

double C2(int i, int j)
{
	return (2*kappa_x[i]);
}

double C3(int i, int j)
{
	return (kappa_x[i]*(-1 - nu[i][j]*dt + 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) + sigma_x[i]*(dt - 0.5*nu[i][j]*pow(dt, 2)));
}

double C4(int i, int j)
{
	return (-1 + nu[i][j]*dt - 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2));
}

double C5(int i, int j)
{ 
	return 2.0;
}

double C6(int i, int j)
{
	return (-1 - nu[i][j]*dt + 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2));
}

double C7(int i, int j)
{
	return (kappa_y[j]*(-1 + nu[i][j]*dt - 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) - sigma_y[j]*(dt - 0.5*nu[i][j]*pow(dt, 2)));
}

double C8(int i, int j)
{
	return (2*kappa_y[j]);
}

double C9(int i, int j)
{
	return (kappa_y[j]*(-1 - nu[i][j]*dt + 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) + sigma_y[j]*(dt - 0.5*nu[i][j]*pow(dt, 2)));
}

double C10(int i, int j)
{
	return (eps_d*(1 + ff)*(-1 + nu[i][j]*dt - 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) - pow(eps_d, 2)*(1 - ff)*(1 - nu[i][j]*dt));
}

double C11(int i, int j)
{ 
	return (2*eps_d*(1 + ff + eps_d*(1 - ff)));
}

double C12(int i, int j) 
{
	return (eps_d*(1 + ff)*(-1 - nu[i][j]*dt + 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) - pow(eps_d, 2)*(1 - ff)*(1 + nu[i][j]*dt));
}

double C13(int i, int j)
{
	return ((1 - ff)*(-1 + nu[i][j]*dt - 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) - eps_d*(1 + ff)*(1 - nu[i][j]*dt));
}

double C14(int i, int j)
{
	return (2*(1 - ff + eps_d*(1 + ff)));
}

double C15(int i, int j)
{
	return ((1 - ff)*(-1 - nu[i][j]*dt + 0.5*pow(plasmaFreq[i][j], 2)*pow(dt, 2)) - eps_d*(1 + ff)*(1 + nu[i][j]*dt));
}

double AZ1(int i, int j)
{
	return C4(i, j)/C1(i, j);
}

double AZ2(int i, int j)
{
	return C5(i, j)/C1(i, j);
}

double AZ3(int i, int j)
{
	return C6(i, j)/C1(i, j);
}

double AZ4(int i, int j)
{ 
	return C2(i, j)/C1(i, j);
}

double AZ5(int i, int j)
{ 
	return C3(i, j)/C1(i, j);
}

double AZ6(int i, int j) 
{
	return C4(i, j)/C7(i, j);
}

double AZ7(int i, int j)
{ 
	return C5(i, j)/C7(i, j);
}

double AZ8(int i, int j)
{ 
	return C6(i, j)/C7(i, j);
}

double AZ9(int i, int j)
{ 
	return C8(i, j)/C7(i, j);
}

double AZ10(int i, int j)
{ 
	return C9(i, j)/C7(i, j);
}

double AZ11(int i, int j)
{ 
	return C13(i, j)/C10(i, j);
}

double AZ12(int i, int j)
{ 
	return C14(i, j)/C10(i, j);
}

double AZ13(int i, int j) 
{
	return C15(i, j)/C10(i, j);
}

double AZ14(int i, int j)
{
	return C11(i, j)/C10(i, j);
}

double AZ15(int i, int j)
{
	return C12(i, j)/C10(i, j);
}

double AX1(int i, int j)
{
	return C1(i, j)/C7(i, j);
}

double AX2(int i, int j)
{
	return C2(i, j)/C7(i, j);
}

double AX3(int i, int j)
{
	return C3(i, j)/C7(i, j);
}

double AX4(int i, int j) 
{
	return C8(i, j)/C7(i, j);
}

double AX5(int i, int j)
{
	return C9(i, j)/C7(i, j);
}

double AY1(int i, int j)
{
	return C7(i, j)/C1(i, j);
}

double AY2(int i, int j)
{
	return C8(i, j)/C1(i, j);
}

double AY3(int i, int j)
{
	return C9(i, j)/C1(i, j);
}

double AY4(int i, int j)
{
	return C2(i, j)/C1(i, j);
}

double AY5(int i, int j)
{
	return C3(i, j)/C1(i, j);
}

void inputSource(int T)
{
	for (int i = 400; i <= 600; i++)
	{
		Gz[i][700] = Gz[i][700] + sin(2*3.14*freqWave*dt*T);
	}
}

void updateElectricField(int nx, int ny, int T, int chunk)
{
/* Gz */
	for (int i = 1; i < nx - 1; i++)
	{
		for (int j = 1; j < ny - 1; j++)
		{
			Gz[i][j] = t_Gz[i][j][1] + (dt/ddx)*(Hy[i + 1][j] - Hy[i - 1][j] - Hx[i][j + 1] + Hx[i][j - 1]);
		}
	}
/* Inserting the input source */
	inputSource(T);

/* Mz */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			Mz[i][j] = AZ1(i, j)*Gz[i][j] + AZ2(i, j)*t_Gz[i][j][1] + AZ3(i, j)*t_Gz[i][j][0] - AZ4(i, j)*t_Mz[i][j][1] - AZ5(i, j)*t_Mz[i][j][0];
		}
	}
	
/* Dz */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			Dz[i][j] = AZ6(i, j)*Mz[i][j] + AZ7(i, j)*t_Mz[i][j][1] + AZ8(i, j)*t_Mz[i][j][0] - AZ9(i, j)*t_Dz[i][j][1] - AZ10(i, j)*t_Dz[i][j][0];
		}
	}

/* Ez */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			Ez[i][j] = AZ11(i, j)*Dz[i][j] + AZ12(i, j)*t_Dz[i][j][1] + AZ13(i, j)*t_Dz[i][j][0] - AZ14(i, j)*t_Ez[i][j][1] - AZ15(i, j)*t_Ez[i][j][0];
		}
	}
}

void updateMagneticField(int nx, int ny, int T, int chunk)
{
/* x components */
/* Tx */
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny - 1; j++)
		{
			Tx[i][j] = t_Tx[i][j][1] - (dt/ddx)*(Ez[i][j + 1] - Ez[i][j - 1]);
		}
	}
	
/* Bx */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			Bx[i][j] = AX1(i, j)*Tx[i][j] + AX2(i, j)*t_Tx[i][j][1] + AX3(i, j)*t_Tx[i][j][0] - AX4(i, j)*t_Bx[i][j][1] - AX5(i, j)*t_Bx[i][j][0];
		}
	}

/* Hx */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			Hx[i][j] = (1/mu_0)*Bx[i][j];
		}
	}
	
/* y components */
/* Ty */
	for (int i = 1; i < nx - 1; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			Tx[i][j] = t_Ty[i][j][1] + (dt/ddx)*(Ez[i + 1][j] - Ez[i - 1][j]);
		}
	}
/* By */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			By[i][j] = AY1(i, j)*Ty[i][j] + AY2(i, j)*t_Ty[i][j][1] + AY3(i, j)*t_Ty[i][j][0] - AY4(i, j)*t_By[i][j][1] - AY5(i, j)*t_By[i][j][0];
		}
	}

/* Hy */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			Hy[i][j] = (1/mu_0)*By[i][j];
		}
	}
}

void resetTempArrays(int nx, int ny, int chunk)
{
/* Resets the values of temporary arrays to 0 */
	#pragma omp parallel for schedule (static, chunk)
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			t_Gz[i][j][0] = t_Gz[i][j][1] = 0.0 + I*0.0;
			t_Mz[i][j][0] = t_Mz[i][j][1] = 0.0 + I*0.0;
			t_Dz[i][j][0] = t_Dz[i][j][1] = 0.0 + I*0.0;
			t_Ez[i][j][0] = t_Ez[i][j][1] = 0.0 + I*0.0;
			t_Tx[i][j][0] = t_Tx[i][j][1] = 0.0 + I*0.0;
			t_Bx[i][j][0] = t_Bx[i][j][1] = 0.0 + I*0.0;
			t_Hx[i][j][0] = t_Hx[i][j][1] = 0.0 + I*0.0;
			t_Ty[i][j][0] = t_Ty[i][j][1] = 0.0 + I*0.0;
			t_By[i][j][0] = t_By[i][j][1] = 0.0 + I*0.0;
			t_Hy[i][j][0] = t_Hy[i][j][1] = 0.0 + I*0.0;
		}
	}
}

void updateTempArrays(int nx, int ny, int chunk)
{
/* Moves older field values from index 1 to index 0, and newer values to index 1 */
	#pragma omp parallel for schedule (static, chunk)
	for (int i = 1; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			t_Gz[i][j][0] = t_Gz[i][j][1];
			t_Gz[i][j][1] = Gz[i][j];
			t_Mz[i][j][0] = t_Mz[i][j][1];
			t_Mz[i][j][1] = Mz[i][j];
			t_Dz[i][j][0] = t_Dz[i][j][1];
			t_Dz[i][j][1] = Dz[i][j];
			t_Ez[i][j][0] = t_Ez[i][j][1];
			t_Ez[i][j][1] = Ez[i][j];
			
			t_Tx[i][j][0] = t_Tx[i][j][1];
			t_Tx[i][j][1] = Tx[i][j];
			t_Bx[i][j][0] = t_Bx[i][j][1];
			t_Bx[i][j][1] = Bx[i][j];
			t_Hx[i][j][0] = t_Hx[i][j][1];
			t_Hx[i][j][1] = Hx[i][j];
			
			t_Ty[i][j][0] = t_Ty[i][j][1];
			t_Ty[i][j][1] = Ty[i][j];
			t_By[i][j][0] = t_By[i][j][1];
			t_By[i][j][1] = By[i][j];
			t_Hy[i][j][0] = t_Hy[i][j][1];
			t_Hy[i][j][1] = Hy[i][j];
		}
	}
}

void mainLoop(int numSteps, int nx, int ny, int chunk)
{
/* The primary calculation and i/o loop */
	double plasmaFreq_current;  // Initial value of plasma frequency
	int T = 0; // Loop counter
	resetTempArrays(nx, ny, chunk);
	
	for (int t = 1; t <= numSteps; t++)
	{
		T = T + 1;
		if ((T > 1000) && (T % 10 == 0))
		{
			plasmaFreq_current = sqrt(pow(e_c, 2)*n(T)/(eps_0*m_e))/(2.0*3.14);
			updatePlasmaFreq(nx, ny, 0, nx, 600, 900, plasmaFreq_current, chunk);
			cout << "!" << T << "!";
		}
		
		updateTempArrays(nx, ny, chunk);
		updateElectricField(nx, ny, T, chunk);
		updateMagneticField(nx, ny, T, chunk);
		
		if (T % 10 == 0)
		{
			writeFields(nx, ny, t);
		}
		cout << T << endl;
	}
}

void writeFields(int nx, int ny, int n)
{
	stringstream ss;
	string filename_Ez, filename_Hx, filename_Hy;
	ss << n;
	filename_Ez = "Out/E_field/Ez" + ss.str() + ".dat";
	
	const char* fname_Ez = filename_Ez.c_str();
	const char* fname_Hx = filename_Hx.c_str();
	const char* fname_Hy = filename_Hy.c_str();

	ofstream out_Ez, out_Hx, out_Hy;
	out_Ez.open(fname_Ez);
	out_Hx.open(fname_Hx);
	out_Hy.open(fname_Hy);
	
	for (int i = 0; i < nx; i = i + 3)
	{
		for (int j = 0; j < ny; j = j + 3)
		{
			out_Ez << i << " " << j << " " << creal(Ez[i][j]) << " " << cimag(Ez[i][j]) << " " << pow(cabs(Ez[i][j]), 2) << endl;
//			out_Hx << i << " " << j << " " << Hx[i][j] << endl;
//			out_Hy << i << " " << j << " " << Hy[i][j] << endl;
		}
		
		out_Ez << endl;
//		out_Hx<<endl;
//		out_Hy<<endl;
	}
	
	out_Ez.close();
	out_Hx.close();
	out_Hy.close();
}

int main()
{
	/* Parallellization Parameters */
	numThreads = 24;
	chunk = 100;
	omp_set_num_threads(numThreads);

	/* Mesh Size */
	nx = 800;
	ny = 1400;

	/* Constants */
	eps_d = 1.0; 						// Dielectric permittivity
	eps_0 = 8.8*pow(10, -12); 			// Permittivity of free space
	mu_0 = 1.2566*pow(10, -7);			// Permeability of free space
	e_c = -1.602176565*pow(10, -19); 	// Electron charge
	m_e = 9.10938291*pow(10, -31); 		// Electron mass
	c = 2.99792458*pow(10, 8); 			// Speed of light

	/* Simulation Parameters */
	ff = 0.00; 							// Fill fraction of filaments
	lambda = 6*pow(10, -4); 			// Wavelength
	ddx = lambda/30.0;					// Mesh spatial step size
	dt = ddx/(2.0*c); 					// Mesh temporal step size
	freqWave = c/lambda; 				// Input frequency
	delta_f = freqWave/1.5;
	delta_x = 50*ddx;
	nu_0 = 0.8*pow(10, 12); 				// Collision Frequency
	plasmaFreq_0 = sqrt(pow(e_c, 2)*n_0/(eps_0*m_e))/(2.0*3.14); // Initial Plasma Frequency
	numSteps = 3500; 					// Number of iterations simulation should execute

	/* Initial Parameters */
	n_0 = 7*pow(10, 16 + 6);
	beta = 0.12*pow(10, -13);
	h = 2.5*pow(10, 7);
	b = beta*n_0;

	/* UPML Parameters */
	upmlSize_x = 40;
	upmlSize_y = 40;
	mp = 4;

	initMesh(nx, ny, 0, nx, 600, 900, plasmaFreq_0, nu_0, upmlSize_x, upmlSize_y, kappa_max, sigma_max, mp, chunk);
	mainLoop(numSteps, nx, ny, chunk);

	return 0;
}
