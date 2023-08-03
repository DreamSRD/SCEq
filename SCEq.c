#include <stdio.h>
#include "math.h"
#include <complex.h>
#include <fftw3.h>
#include "SCEq.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

double Random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}

void CheckFile(FILE* file_name)
{
	if (file_name == NULL)
	{
		printf("Cannot open or find the file \n");
		exit(0);
	}
}

void Set_new_data()
{
	// C (k, t) 

	for (int j = 0; j < N; j++)
	{
		//Ck[j] = Random(1e-11, 1e-10) + I * Random(1e-11, 1e-10);
		Ck[j] = 0.0 + I * 0.0;
	}

	Ck[100] = 1.0 + I * 0.0;

	printf("%s\n", "New data was successfully set");
}

void Load_data()
{
	sprintf(Name, "./RestartFiles//Single_Breather_mu01_k0_100");
	pfile = fopen (Name,"rb");
	CheckFile(pfile);
	number_of_numbers = fread(&Ck[0], sizeof(fftw_complex), N, pfile);
	fclose(pfile);

	printf("Restart File %s has been succesfully read \n", Name);

/*	pfile = fopen ("C_k_data","rb");
	number_of_numbers = fread(&Ck[0], sizeof(fftw_complex), N, pfile);
	fclose(pfile); */
}

void Set_data()
{
	Bx = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Bk = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	Cx = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Ck = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	CCx = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	CCk = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	Cdx= (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Cdk = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	CkOld = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	CkNew = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	RHSk1 = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	RHSk2 = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	RHSk3 = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	RHSk4 = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	CCCdx = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	CCCdk = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	K = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	KCCx = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	KCCk = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	KCCCx = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	KCCCk = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));


	Eta_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Eta_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Eta_2_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Eta_2_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	Psi_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Psi_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	Psi_2_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Psi_2_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Psi_2m_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Psi_2m_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Psi_2p_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Psi_2p_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	Mu_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Mu_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
	Mu_2_x = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));	
	Mu_2_k = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

	t = 0.0;
	h = (b - a)/N;
	Dk = 2 * M_PI/ (b - a);

	k0dim = Dk * k0;
	w0 = sqrt(g * k0dim);
	T0 = 2.0 * M_PI/w0;



	// C (k, t) 

	for ( int j = 0; j < N; j++)
	{
		Ck[j] = 0.0 + I * 0.0;
	} 

	//Set_new_data();

	Load_data();


	//  ik massive initializing

	for ( int j = N/2 + 1; j < N; j++)
	{
		K[j] = cabs(j - N) * Dk;
	} 

	for ( int j = 0; j < N/2; j++)
	{
		K[j] = j * Dk;
		Eta_k_average[j] = 0.0; // Для усреднения спектра поверхности
	}

	for (int j = 0; j < N; j++)
	{
		CkOld[j] = Ck[j];
	}

	K[N/2] = 0.0 + I * 0.0;
}


void Closing()
{

	fftw_free(Cx);
	fftw_free(Ck);

	fftw_free(CCx);
	fftw_free(CCk);

	fftw_free(Cdx);
	fftw_free(Cdk);
	fftw_free(CkOld);
	fftw_free(CkNew);

	fftw_free(RHSk1);
	fftw_free(RHSk2);
	fftw_free(RHSk3);
	fftw_free(RHSk4);

	fftw_free(K);
	fftw_free(CCCdx);
	fftw_free(CCCdk);

	fftw_free(KCCx);
	fftw_free(KCCk);
	fftw_free(KCCCx);
	fftw_free(KCCCk);

	fftw_free(Eta_x);
	fftw_free(Eta_k);
	fftw_free(Eta_2_x);
	fftw_free(Eta_2_k);

	fftw_free(Psi_x);
	fftw_free(Psi_k);

	fftw_free(Psi_2_x);
	fftw_free(Psi_2_k);
	fftw_free(Psi_2m_k);
	fftw_free(Psi_2m_x);
	fftw_free(Psi_2p_k);
	fftw_free(Psi_2p_x);

	fftw_free(Mu_x);
	fftw_free(Mu_k);
	fftw_free(Mu_2_x);
	fftw_free(Mu_2_k);

}

void CreatePlans()
{
	Cx_to_Ck = fftw_plan_dft_1d(N,&Cx[0],&Ck[0],FFTW_FORWARD,FFTW_ESTIMATE);
 	Ck_to_Cx = fftw_plan_dft_1d(N,&Ck[0],&Cx[0],FFTW_BACKWARD,FFTW_ESTIMATE);

 	Bk_to_Bx = fftw_plan_dft_1d(N,&Bk[0],&Bx[0],FFTW_BACKWARD,FFTW_ESTIMATE);

 	CCx_2_CCk = fftw_plan_dft_1d(N,&CCx[0],&CCk[0],FFTW_FORWARD,FFTW_ESTIMATE);
 	Cdk_to_Cdx = fftw_plan_dft_1d(N,&Cdk[0],&Cdx[0],FFTW_BACKWARD,FFTW_ESTIMATE);

 	KCCk_to_KCCx = fftw_plan_dft_1d(N,&KCCk[0],&KCCx[0],FFTW_BACKWARD,FFTW_ESTIMATE);
 	KCCCx_2_KCCCk = fftw_plan_dft_1d(N,&KCCCx[0],&KCCCk[0],FFTW_FORWARD,FFTW_ESTIMATE);
 	CCCdx_2_CCCdk = fftw_plan_dft_1d(N,&CCCdx[0],&CCCdk[0],FFTW_FORWARD,FFTW_ESTIMATE);

 	Eta_k_to_Eta_x = fftw_plan_dft_1d(N, &Eta_k[0], &Eta_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);
 	Psi_k_to_Psi_x = fftw_plan_dft_1d(N, &Psi_k[0], &Psi_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);

 	Eta_2_x_to_Eta_2_k = fftw_plan_dft_1d(N, &Eta_2_x[0], &Eta_2_k[0], FFTW_FORWARD,FFTW_ESTIMATE);
 	Eta_2_k_to_Eta_2_x = fftw_plan_dft_1d(N, &Eta_2_k[0], &Eta_2_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);

 	Psi_2_x_to_Psi_2_k = fftw_plan_dft_1d(N, &Psi_2_x[0], &Psi_2_k[0], FFTW_FORWARD,FFTW_ESTIMATE);
 	Psi_2_k_to_Psi_2_x = fftw_plan_dft_1d(N, &Psi_2_k[0], &Psi_2_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);

 	Psi_2m_k_to_Psi_2m_x = fftw_plan_dft_1d(N, &Psi_2m_k[0], &Psi_2m_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);
 	Psi_2p_k_to_Psi_2p_x = fftw_plan_dft_1d(N, &Psi_2p_k[0], &Psi_2p_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);

 	Mu_k_to_Mu_x = fftw_plan_dft_1d(N, &Mu_k[0], &Mu_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);
 	Mu_2_k_to_Mu_2_x = fftw_plan_dft_1d(N, &Mu_2_k[0], &Mu_2_x[0], FFTW_BACKWARD,FFTW_ESTIMATE);
}

void DestroyPlans()
{
	fftw_destroy_plan(Cx_to_Ck);
	fftw_destroy_plan(Ck_to_Cx);

	fftw_destroy_plan(Bk_to_Bx);

	fftw_destroy_plan(CCx_2_CCk);
	fftw_destroy_plan(Cdk_to_Cdx);

	fftw_destroy_plan(CCCdx_2_CCCdk);
	fftw_destroy_plan(KCCCx_2_KCCCk);
	fftw_destroy_plan(KCCk_to_KCCx);

	fftw_destroy_plan(Eta_k_to_Eta_x);
	fftw_destroy_plan(Psi_k_to_Psi_x);
	
	fftw_destroy_plan(Eta_2_x_to_Eta_2_k);
	fftw_destroy_plan(Eta_2_k_to_Eta_2_x);

	fftw_destroy_plan(Psi_2_x_to_Psi_2_k);
	fftw_destroy_plan(Psi_2_k_to_Psi_2_x);

	fftw_destroy_plan(Psi_2m_k_to_Psi_2m_x);
	fftw_destroy_plan(Psi_2p_k_to_Psi_2p_x);

	fftw_destroy_plan(Mu_k_to_Mu_x);
	fftw_destroy_plan(Mu_2_k_to_Mu_2_x);

}

void Normalize(fftw_complex *Arr_x)
{
	for (int j = 0; j < N; j++)
	{
		Arr_x[j] = Arr_x[j]/N;
	}

	Arr_x[0] = 0.0 + I * 0.0;
	Arr_x[N/2] = 0.0 + I * 0.0;
}

double MaxFinder(fftw_complex *Ck)
{
	double Max = 0.0;
	for (int j = 0; j < N; j++)
	{
		if (cabs(Ck[j]) > Max)
		{
			Max = cabs(Ck[j]);
			Max_int = j;
		}
	}

	return Max;
}

void Damping()
{
	double Gamma = 0.001;

/*	double Gamma_exp = 0.5;
	double Sigma = 0.075;
	double x;
	double L = b - a;*/

	for ( int j = 0; j < N; j++)
	{
		Cx[j] = Cx[j] * exp( - Gamma * pow(cos(M_PI * (a + j * h)/(b - a)), 6.0) * tau); ///// Затухание с косинусом
	} 

/*	for ( int j = 0; j < N; j++)
	{
		x = a + j * h;
		Cx[j] = Cx[j] * exp( -Gamma_exp * exp(x * (x - L)/(Sigma * L)) * tau); /////// Затухание с экспонентой
	} */

	fftw_execute(Cx_to_Ck);
	Normalize(&Ck[0]);

	for (int j = 0; j < N; j++)
	{
		CkOld[j] = Ck[j];
	}

	for (int j = N/2; j < N; j++)
	{
		Ck[j] = 0.0 + I * 0.0;
		CkOld[j] = 0.0 + I * 0.0;
	}

}


void CalculateStatC(fftw_complex* C, unsigned long long int* HR, unsigned long long int* Habs, unsigned long long int* Habs2, double h)
{
	double N_curr = 0.0;
	int N_curr_PDF = 0;

	for (int j = 0; j < N; j++)
	{	
		N_curr = creal(C[j])/h;
		N_curr_PDF = (int)(N_curr + (N_curr<0.0?-0.5:+0.5)) + ((Nstat/2.0)-1);

		if(N_curr_PDF < 0)
		{
			HR[0] = HR[0] + 1;
		}

		else if(N_curr_PDF >= Nstat)
		{
			HR[Nstat - 1] = HR[Nstat - 1] + 1;	
		}

		else
		{
			HR[N_curr_PDF] = HR[N_curr_PDF] + 1;	
		}
	

		N_curr = sqrt(pow(creal(C[j]), 2.0) + pow(cimag(C[j]), 2.0) )/h;
		N_curr_PDF = (int)(N_curr + (N_curr<0.0?-0.5:+0.5)) + ((Nstat/2.0)-1);


		if(N_curr_PDF < Nstat)
		{
			Habs[N_curr_PDF] = Habs[N_curr_PDF] + 1;
		}

		else
		{

			Habs[Nstat - 1] = Habs[Nstat - 1] + 1;	
		}

		N_curr = (pow(creal(C[j]), 2.0) + pow(cimag(C[j]), 2.0))/h;
		N_curr_PDF = (int)(N_curr + (N_curr<0.0?-0.5:+0.5)) + ((Nstat/2.0)-1);


		if(N_curr_PDF < Nstat)
		{
			Habs2[N_curr_PDF] = Habs[N_curr_PDF] + 1;
		}

		else
		{

			Habs2[Nstat - 1] = Habs[Nstat - 1] + 1;	
		}

	}
	

}

void CalculateStatSurf(fftw_complex* Eta_x, unsigned long long int* H_surf, double h)
{
	double N_curr = 0.0;
	int N_curr_PDF = 0;
	for (int j = 0; j < N; j++)
	{

		N_curr = creal(Eta_x[j])/h;
		N_curr_PDF = (int)(N_curr + (N_curr<0.0?-0.5:+0.5)) + ((Nstat/2.0) - 1);

		if(N_curr_PDF < 0)
		{
			H_surf[0] = H_surf[0] + 1;
		}

		else if(N_curr_PDF >= Nstat - 1)
		{
			H_surf[Nstat - 1] = H_surf[Nstat - 1] + 1;	
		}

		else
		{
			H_surf[N_curr_PDF] = H_surf[N_curr_PDF] + 1;	
		}
	}	

}

void RHSNonlin()
{
	for (int j = 0; j < N/2; j++)
	{
		Cdk[j] = I * j * Dk * Ck[j];
	}

	fftw_execute(Cdk_to_Cdx);

	for (int j = 0; j < N; j++)
	{
			CCCdx[j] = CCx[j] * Cdx[j];
	}

	fftw_execute(CCCdx_2_CCCdk);

	Normalize(&CCCdk[0]);
}

void RHSAdvect()
{
	fftw_execute(CCx_2_CCk);

	Normalize(&CCk[0]);

	for (int j = 0; j < N; j++)
	{
		KCCk[j] = K[j] * CCk[j];
	}

	fftw_execute(KCCk_to_KCCx);

	for (int j = 0; j < N; j++)
	{
		KCCCx[j] = KCCx[j] * Cx[j];
	}

	fftw_execute(KCCCx_2_KCCCk);

	Normalize(&KCCCk[0]);


}

void RHS(fftw_complex *RHS)
{
	fftw_execute(Ck_to_Cx);

	for (int j = 0; j < N; j++)
	{
		CCx[j] = Cx[j] * conj(Cx[j]);
	}

	RHSNonlin();

	RHSAdvect();
	
	for (int j = 0; j < N/2; j++)
	{
		RHS[j] = -I * sqrt(g * j * Dk) * Ck[j] + V * I * j * Dk * Ck[j] - j * Dk * CCCdk[j] + I * j * Dk * KCCCk[j];
	}

}

void Shift(fftw_complex *Shift, double step)
{
	for (int j = 0; j < N; j++)
	{
		Ck[j] = CkOld[j] + step * Shift[j];
	}

}

void Single_RK4_step()
{
	for ( int j = 0; j < N; j++)
	{
		CkNew[j] = CkOld[j] + (tau/6.0) * (RHSk1[j] + 2.0 * RHSk2[j] + 2.0 * RHSk3[j] + RHSk4[j]);
	}

	for ( int j = 0; j < N; j++)
	{
		CkOld[j] = CkNew[j];
		Ck[j] = CkNew[j];
	}
}


double NumberOfWaves()
{
	double Num = 0.0;

	for (int j = 1; j < N/2; j++)
	{
		Num = Num + creal((Ck[j] * conj(Ck[j]))/(j * Dk));
	}

	return Num;
}

double Momentum()
{
	double Mom = 0.0;

	for (int j = 0; j < N/2; j++)
	{
		Mom = Mom + creal(conj(Ck[j]) * Ck[j]);
	}

	return Mom;
}

double TotalEnergy()
{
	double E1 = 0.0;
	double E2 = 0.0;
	double E3 = 0.0;

	for (int j = 1; j < N/2; j++)
	{
		E1 = E1 + creal(conj(Ck[j]) * sqrt(g * j * Dk)/(j * Dk) * Ck[j]);
	}

	for (int j = 0; j < N; j++)
	{
		E2 = E2 + cimag(CCx[j] * conj(Cx[j]) * Cdx[j]);
	}

	for (int j = 0; j < N; j++)
	{
		E3 = E3 + creal(CCx[j] * KCCx[j]);
	}
	
	return E1 + (E2 - E3)/(2.0 * N);
} 

void Print_Cx()
{
	h = (b - a)/N;
	sprintf(Name, "./Cx//Cx_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	fprintf(pfile, "\t%s\t\t\t  %s\t\t\t N_av = %d\n", "[x]", "|c(x)|", N_av);

	for (int j = 0; j < N; j++)
	{
		fprintf(pfile, "%e\t%e\n", a + j * h, cabs(Cx[j]));
	}

	fclose(pfile);
}

void Print_Ck()
{
	sprintf(Name, "./Ck//Ck_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	fprintf(pfile, "%s\t\t%s\t\t N_av = %d\n", "[k]", "|c(k)|", N_av);

	for (int j = 0; j < N/2; j++)
	{
		fprintf(pfile, "%d\t%e\n", j, cabs(Ck[j]));
	}

	fclose(pfile);
}

void Increment()
{
	pfile = fopen("Gamma.txt", "a");
	CheckFile(pfile);

	fprintf (pfile, "%e\t%e\t%e\n", t, cabs(Ck[90]),cabs(Ck[110]));

	fclose(pfile);
}

void Print_NPE()
{
	pfile = fopen("NPE.txt", "a");
	CheckFile(pfile);

	if ( nstep == 0)
	{
		fprintf (pfile, "%s\t\t%s\t%s\t%s\n", " Time [t] ", " Number of Waves [N] ", " Total Momentum [P] ", " Total Energy [E] ");
	}
	
	fprintf (pfile, "%e\t%.15e\t%.15e\t%.15e\n", t,  NumberOfWaves(), Momentum(), TotalEnergy());
	fclose(pfile);
}

void Eta_Psi_Mu_recover_1()
{	

	for ( int j = 1; j < N/2; j ++)
	{
		Eta_k[j] = 1.0/(sqrt(2.0) * pow(g,0.25)) * pow(K[j], -0.25) * Ck[j];
		Psi_k[j] = -I * pow(g,0.25)/sqrt(2.0) * pow(K[j], -0.75) * Ck[j];
		Eta_k[N - j] = 1.0/(sqrt(2.0) * pow(g,0.25)) * pow(K[j], -0.25) * conj(Ck[j]); //// %%%% Первый порядок правильный, проверялся с FFTW на волне
		Psi_k[N - j] = I * pow(g,0.25)/sqrt(2.0) * pow(K[j], -0.75) * conj(Ck[j]); ////// eta_k = c_k + conj(c_(-k)) !!!
	}

	Eta_k[0] = 0.0 + I * 0.0;
	Eta_k[N/2] = 0.0 + I * 0.0;

	Psi_k[0] = 0.0 + I * 0.0;
	Psi_k[N/2] = 0.0 + I * 0.0;

	// int binnum;
	// pfile = fopen("./eta_k", "w+b"); 
	// CheckFile(pfile);
	// binnum = fwrite(&Eta_k[0], sizeof(fftw_complex), Nstat, pfile);
	// fclose(pfile);

	// pfile = fopen("./psi_k", "w+b"); 
	// CheckFile(pfile);
	// binnum = fwrite(&Psi_k[0], sizeof(fftw_complex), Nstat, pfile);
	// fclose(pfile);

	fftw_execute(Eta_k_to_Eta_x);
	fftw_execute(Psi_k_to_Psi_x);

	////////////////////////////////////////////////////////////////////////////////////////////// Steepness (Mu)

	for ( int j = N/2 + 1; j < N; j++)
	{
		Mu_k[j] = I * (j - N) * Dk * Eta_k[j];
	}
	for ( int j = 0; j < N/2; j++)
	{
		Mu_k[j] = I * j * Dk * Eta_k[j];
	}

	fftw_execute(Mu_k_to_Mu_x);

}

void Eta_Psi_Mu_recover_2()
{

	for ( int j = 1; j < N/2; j ++)
	{
		Eta_2_k[j] = pow(K[j], -0.25) * Ck[j];
		Eta_2_k[N - j] = - pow(K[j], -0.25) * conj(Ck[j]);
		Psi_2p_k[j] = pow(K[j], 0.25) * Ck[j];
		Psi_2m_k[j] = pow(K[j], -0.25) * Ck[j];
	}

 	Eta_2_k[0] = 0.0 + I * 0.0;
	Eta_2_k[N/2] = 0.0 + I * 0.0;

	Psi_2p_k[0] = 0.0 + I * 0.0;
	Psi_2p_k[N/2] = 0.0 + I * 0.0;
	Psi_2m_k[0] = 0.0 + I * 0.0;
	Psi_2m_k[N/2] = 0.0 + I * 0.0;

	fftw_execute(Eta_2_k_to_Eta_2_x);

	fftw_execute(Psi_2m_k_to_Psi_2m_x);
	fftw_execute(Psi_2p_k_to_Psi_2p_x);

	for ( int j = 0; j < N; j ++)
	{
		Eta_2_x[j] = Eta_2_x[j] * Eta_2_x[j];
		Psi_2_x[j] = 2.0 * creal(Psi_2m_x[j] * conj(Psi_2p_x[j]));
	}

	fftw_execute(Eta_2_x_to_Eta_2_k);
	Normalize(&Eta_2_k[0]);

	fftw_execute(Psi_2_x_to_Psi_2_k);
	Normalize(&Psi_2_k[0]);

	for ( int j = 0; j < N; j ++)
	{
		Eta_2_k[j] = (K[j] * Eta_2_k[j])/(4.0 * sqrt(g));
	}

	for ( int j = 0; j < N/2; j ++)
	{
		Psi_2_k[j] = I/2.0 * Psi_2_k[j];
	}

	for ( int j = N/2 + 1; j < N; j ++)
	{
		Psi_2_k[j] = -I/2.0 * Psi_2_k[j];
	}

	Psi_2_k[0] = 0.0 + I * 0.0;
	Psi_2_k[N/2] = 0.0 + I * 0.0;

	fftw_execute(Eta_2_k_to_Eta_2_x);
	fftw_execute(Psi_2_k_to_Psi_2_x);

	for ( int j = 0; j < N; j ++)
	{
		Psi_2_x[j] = cimag(Psi_2m_x[j] * Psi_2p_x[j]) + Psi_2_x[j];
	}

	fftw_execute(Psi_2_x_to_Psi_2_k); /// Для итогового eta_k второго порядка
	Normalize(&Psi_2_k[0]);


	////////////////////////////////////////////////////////////////////////////////////////////// Steepness (Mu)

	for ( int j = N/2 + 1; j < N; j++)
	{
		Mu_2_k[j] = I * (j - N) * Dk * Eta_2_k[j];
	}

	for ( int j = 0; j < N/2; j++)
	{
		Mu_2_k[j] = I * j * Dk * Eta_2_k[j];
	}

	fftw_execute(Mu_2_k_to_Mu_2_x);

	//////////////////////////////////////////////////////////////////////////////////////////////


	for (int j = 0; j < N; j++)
	{
		Eta_2_x[j] = Eta_x[j] + Eta_2_x[j];
		Eta_2_k[j] = Eta_k[j] + Eta_2_k[j];

		Psi_2_x[j] = Psi_x[j] + Psi_2_x[j];
		Psi_2_k[j] = Psi_k[j] + Psi_2_k[j];

		Mu_2_x[j] = Mu_x[j] + Mu_2_x[j];
	}

}

double AverageSteep()
{
	double Sum = 0.0;

	for (int j = 0; j < N/2; j++)
	{
		Sum = Sum + (j * Dk * j * Dk) * Eta_k[j] * conj(Eta_k[j]);
	}

	return sqrt(Sum);
}

void AverageSurface_k()
{
	for (int j = 0; j < N/2; j++)
	{
		Eta_k_average[j] = Eta_k_average[j] + cabs(Eta_k[j]);
	}
}

void Print_Eta_Psi_Mu()
{
	Eta_Psi_Mu_recover_1();
	Eta_Psi_Mu_recover_2();

	//AverageSurface_k();	

	h = (b - a)/N;
	sprintf(Name, "./Eta_Psi_x//Eta_Psi_x_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	for (int j = 0; j < N; j++)
	{
		fprintf(pfile, "%e\t%e\t%e\t%e\t%e\n", a + j * h, creal(Eta_x[j]), creal(Eta_2_x[j]), creal(Psi_x[j]), creal(Psi_2_x[j]));
	}

	fclose(pfile);

	sprintf(Name, "./Eta_Psi_k//Eta_Psi_k_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);


	fprintf(pfile, "%s\t\t%s\t\t\t%s\t\t\t N_av = %d\n", "[k]", "|Eta(k)|", "|Psi(k)|", N_av);


	for (int j = N/2; j < N; j++)
	{
		fprintf(pfile, "%d\t%e\t%e\n", j - N, cabs(Eta_2_k[j]), cabs(Psi_2_k[j]));
	}

	for (int j = 0; j < N/2; j++)
	{
		fprintf(pfile, "%d\t%e\t%e\n", j, cabs(Eta_2_k[j]), cabs(Psi_2_k[j]));
	}

	fclose(pfile);

	sprintf(Name, "./Mu//Mu_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	fprintf(pfile, "\t%s\t\t\t%s\t\t\t N_av = %d\n", "x", "Mu(x)", N_av);

	for (int j = 0; j < N; j++)
	{
		fprintf(pfile, "%e\t%e\n", a + j * h, creal(Mu_x[j]));
	}

	fclose(pfile);
}

void RV_RestartFiles()
{
	sprintf(Name, "./RV_restart//data_etafin_%d", nstep/N_av);
	pfile = fopen(Name, "w+b"); 
	CheckFile(pfile);

	number_of_numbers = fwrite(&Eta_2_k[0], sizeof(fftw_complex), N, pfile);
	fclose(pfile);

	sprintf(Name, "./RV_restart//data_psifin_%d", nstep/N_av);
	pfile = fopen(Name, "w+b"); 
	CheckFile(pfile);

	number_of_numbers = fwrite(&Psi_2_k[0], sizeof(fftw_complex), N, pfile);
	fclose(pfile);
}



void StatNormalize(double* Arr, double h)
{
	double Sum = 0.0;
	for (int j = 0; j < Nstat; j++)
	{
		Sum = Sum + Arr[j];
	}

	Sum = Sum * h;

	for (int j = 0; j < Nstat; j++)
	{
		Arr[j] = Arr[j]/Sum;
	}
}

/*void PrintStatC(double Amp0)
{
	sprintf(Name, "./Statistics//StatC_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	fprintf(pfile, "\t%s\t\t%s\t\t\t%s\t\t%s\t\t\t N_av = %d\n", "c/c0", "PDF(c)", "PDF(|c|)", "PDF(c * c)", N_av);

	for (int j = 0; j < Nstat; j++)
	{
		HC_double[j] = HC_int[j]/(double)(N_av);
		HCabs_double[j] = HCabs_int[j]/(double)(N_av);
		HCabs2_double[j] = HCabs2_int[j]/(double)(N_av);
	}

	StatNormalize(&HC_double[0], hstat1);
	StatNormalize(&HCabs_double[0], hstat1);
	StatNormalize(&HCabs2_double[0], hstat1);

	for (int j = 0; j < Nstat; j++)
	{
		fprintf(pfile, "%e\t%e\t%e\t%e\n", (l_lim1 + j * hstat1)/Amp0, HC_double[j], HCabs_double[j], HCabs2_double[j]);
	}

	fclose(pfile);

}*/

/*
void PrintStatSurf(double Surf0)
{
	sprintf(Name, "./Statistics//StatSurf_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	fprintf(pfile, "\t%s\t%s\t\t\t N_av = %d\n", "Eta/Eta0", "PDF(Eta)", N_av);

	for (int j = 0; j < Nstat; j++)
	{
		HSurf_double[j] = HSurf_int[j]/(double)(N_av);
	}

	StatNormalize(&HSurf_double[0], hstat2);

	for (int j = 0; j < Nstat; j++)
	{
		fprintf(pfile, "%e\t%e\n", (l_lim2 + j * hstat2)/Surf0, HSurf_double[j]);
	}

	fclose(pfile);

}*/

/*void ResetStat()
{
	for (int j = 0; j < Nstat; j++)
	{
		HC_int[j] = 0;
		HCabs_int[j] = 0;
		HCabs2_int[j] = 0;
		HC_double[j] = 0.0;
		HCabs_double[j] = 0.0;
		HCabs2_double[j] = 0.0;

		HSurf_int[j] = 0;
		HSurf_double[j] = 0.0;
	}
}*/


void PrintAverageSurf_k()
{

	for (int j = 0; j < N/2; j++)
	{
		Eta_k_average[j] = Eta_k_average[j]/(double)(N_av);
	}

	sprintf(Name, "./BitStat//Eta_k_%d", nstep/N_av);
	pfile = fopen(Name, "w+b"); 
	CheckFile(pfile);
	number_of_numbers = fwrite(&Eta_k_average[0], sizeof(fftw_complex), Nstat, pfile);
	fclose(pfile);

	sprintf(Name, "./Statistics//Eta_k_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	fprintf(pfile, "%s\t%s\t\t\t N_av = %d\n", "[k]", "|Eta(k)|", N_av);

	for (int j = 0; j < N/2; j++)
	{
		fprintf(pfile, "%d\t%e\n", j, Eta_k_average[j]);
	}

	fclose(pfile);

	for (int j = 0; j < N/2; j++)
	{
		Eta_k_average[j] = 0.0;
	}
}


void CreateDirectories()
{
	mkdir("Bx", 0777);
	mkdir("Cx", 0777);
	mkdir("Ck", 0777);
	mkdir("Eta_Psi_x", 0777);
	mkdir("Eta_Psi_k", 0777);
	mkdir("CheckPoints", 0777);
	mkdir("Mu", 0777);
	mkdir("Statistics", 0777);
}

void Clean()
{
	remove("NPE.txt");
	remove("StatMom.txt");
	system("exec rm -r Cx/*");
	system("exec rm -r Bx/*");
	system("exec rm -r Ck/*");
	system("exec rm -r Eta_Psi_x/*");
	system("exec rm -r Eta_Psi_k/*");
	system("exec rm -r RV_restart/*");
	system("exec rm -r Mu/*");
	system("exec rm -r Statistics/*");
	system("exec rm -r BitStat/*");
}

void CheckPoint()
{
	sprintf(Name, "./CheckPoints//C_k_data_%d", nstep/N_av);
	pfile = fopen (Name,"w+b");
	CheckFile(pfile);
	number_of_numbers = fwrite(&Ck[0], sizeof(fftw_complex), N, pfile);
	fclose(pfile);
}

/*void BitStat()
{
	sprintf(Name, "./BitStat//BitStatReC_%d", nstep/N_av);
	pfile = fopen (Name,"w+b");
	CheckFile(pfile);
	int number_of_numbers = fwrite(&HC_double[0], sizeof(fftw_complex), Nstat, pfile);
	fclose(pfile);

	sprintf(Name, "./BitStat//BitStatAbsC_%d", nstep/N_av);
	pfile = fopen (Name,"w+b");
	CheckFile(pfile);
	number_of_numbers = fwrite(&HCabs_double[0], sizeof(fftw_complex), Nstat, pfile);
	fclose(pfile);

	sprintf(Name, "./BitStat//BitStatSurf_%d", nstep/N_av);
	pfile = fopen (Name,"w+b");
	CheckFile(pfile);
	number_of_numbers = fwrite(&HSurf_double[0], sizeof(fftw_complex), Nstat, pfile);
	fclose(pfile);
}*/


/*void StatMake()
{
	double k0dim = g/(4.0 * V * V);
	// ############################################ Гистограмма
	u_lim1 = 4.0 * MaxFinder(&Cx[0]);
	u_lim2 = 4.0 * sqrt(2.0) * pow(k0dim, -0.25) * MaxFinder(&Cx[0])/pow(g,0.25);
	l_lim1 = -4.0 * MaxFinder(&Cx[0]);
	l_lim2 = -4.0 * sqrt(2.0) * pow(k0dim, -0.25) * MaxFinder(&Cx[0])/pow(g,0.25);
	hstat1 = (u_lim1 - l_lim1)/Nstat;
	hstat2 = (u_lim2 - l_lim2)/Nstat;

	for (int j = 0; j < Nstat; j++)
	{
		CR[j] = (l_lim1 + j * hstat1);
		SurfR[j] = (l_lim2 + j * hstat2);
		HC_int[j] = 0;
		HCabs_int[j] = 0;
		HCabs2_int[j] = 0;
		HSurf_int[j] = 0;
		HC_double[j] = 0.0;
		HSurf_double[j] = 0.0;
	}
	// ############################################
}*/

double Check_Break(fftw_complex *Ck, double A)
{
	double Average_Ck = 0.0;
	for (int j = N/2 - 100; j < N/2; j++)
	{
		Average_Ck = Average_Ck + cabs(Ck[j]);
	}

	Average_Ck = Average_Ck/100.0;

	if (Average_Ck > 1e-4 * A)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void b_Calculate()
{
	double h = (b - a)/N;
	for (int j = 0; j < N; j++)
	{
		Bk[j] = 0.0 + I * 0.0;
	}

	for (int j = 1; j < N/2; j++)
	{
		Bk[j] = Ck[j]/(sqrt(j * Dk));
	}

	fftw_execute(Bk_to_Bx);
	
	sprintf(Name, "Bx//bx_%d.txt", nstep/N_av);
	pfile = fopen(Name, "w");
	CheckFile(pfile);

	if (nstep == 0)
	{
		fprintf(pfile, "\t%s\t\t\t%s\t\t\t%s\t\t\t N_av = %d\n", "x", "Re[b(x)]", "n(x)", N_av);
	}

	for (int j = 0; j < N; j++)
	{	
		fprintf(pfile, "%e\t%e\t%e\n", a + j * h, creal(Bx[j]), creal(Bx[j] * conj(Bx[j])));
	}

	fclose(pfile);
}
////////////////////////////////////////////////////////////////////////////////////////////////////// STAT MOM CALCULATIONS

double Expectation(double* Arr)
{
	double Sum = 0.0;

	for (int j = 0; j < Nstat; j++)
	{
		Sum = Sum + Arr[j];
	}

	return Sum/Nstat;
}

double CentralMom(double* Arr, int power)
{
	double Sum = 0.0;
	for (int j = 0; j < Nstat; j++)
	{
		Sum = Sum + pow(Arr[j] - Expectation(&Arr[0]), power);
	}

	return Sum/Nstat;
}

double RMS(double* Arr)
{
	double D = CentralMom(&Arr[0], 2.0);
	return sqrt(D);
}

double SkewNess(double* Arr)
{
	double mu3 = CentralMom(&Arr[0], 3.0);

	double sigma3 = pow(RMS(&Arr[0]),3.0);

	return mu3/sigma3;
}

double Kurtosis(double* Arr)
{
	double mu4 = CentralMom(&Arr[0], 4.0);
	double sigma4 = pow(RMS(&Arr[0]),4.0);

	return (mu4/sigma4) - 3.0;
}

/*void Print_StatMoments()
{
	pfile = fopen("StatMom.txt", "a");
	CheckFile(pfile);

	if (nstep == 0)
	{
		fprintf(pfile, "%s\t\t%s\t%s\t%s\t%s\n", "Time, [t]", "Expectation, E[X]", "RMS S[X]", "Skewness, G[X]", "Kurtosis, K[X]");
	}

	else
	{
		fprintf(pfile, "%e\t%e\t%e\t%e\t%e\n", t, Expectation(&HSurf_double[0]), RMS(&HSurf_double[0]), SkewNess(&HSurf_double[0]), Kurtosis(&HSurf_double[0]));
	}
	
	fclose(pfile);
}*/



void Step_Check()
{
		double taumax_Ny = M_PI/(4.0 * sqrt(g * Dk * N/2.0)); // Критерий по частоте Найквиста
		double taumax_C = 2.0 * (b - a) * sqrt(Dk)/(N * sqrt(g)); // Критерий Куранта

		//printf(" Tau = %e\t  Tau_max = %e\t Tau_max_1 =  %e\n", tau, taumax, taumax1);

		if (tau > taumax_Ny)
		{
			printf("Warning! The time step is too high. Please decrease the step to avoid possible calculation instability. \n Tau = %e\t  Tau_Nyquist = %e\t  Tau_Courant = %e\n", tau, taumax_Ny, taumax_C);
			exit(0);
		}
		else
		{
			printf("Time step is adequate. \n Tau = %e\t Tau_Nq = %e\t Tau_Crnt = %e\n\n", tau, taumax_Ny, taumax_C);
		}
}


double Av_int_finder(fftw_complex* Cx)
{
	double x_av = 0.0;
	int Av_int = 0;

	double Sum = 0.0;
	double Norm = 0.0;

	for (int j = 0; j < N - 1; j++)
	{
		Sum = Sum + (a + j * h) * cabs(Cx[j]) * h;
		Norm = Norm + (cabs(Cx[j]) * h);
	}

	x_av = Sum/Norm;

	Av_int = (int)(x_av/h);
	
	return Av_int;

}

double Velocity_Correction(int Calculation_number, double delta_t, int Number_of_measurements)
{
		if (Calculation_number < Number_of_measurements)
		{			
			Coordinate[Calculation_number] = a + Av_int_finder(&Cx[0]) * h;
		}

		if (Calculation_number == Number_of_measurements)
		{
			for (int j = 0; j < Number_of_measurements - 1; j++)
			{
				Velocity[j] = (Coordinate[j + 1] - Coordinate[j])/delta_t;
				V_br = V_br + Velocity[j];
			}

			V_br = V_br/(Number_of_measurements - 1);
			//V = V + V_br; // Velocity correction when active and velocity measurement when off
			printf(" V = %.8e\n", V_br);
		}
	
}

void Calculating_Evolution(double Time)
{
	int N_msr = 0; // Счётчик измерения скорости

	double V_av = 0.0;

	double Ak = MaxFinder(&Ck[0]); // Максимум амплитуды в k-пространстве для определения опрокидывания

	T = 1000.0; // Характерный период

	while (t < Time) // 
		{
			RHS(&RHSk1[0]);
			Shift(&RHSk1[0], tau/2.0);

			RHS(&RHSk2[0]);
			Shift(&RHSk2[0], tau/2.0);

			RHS(&RHSk3[0]);
			Shift(&RHSk3[0], tau);

			RHS(&RHSk4[0]);

			Single_RK4_step();


			////////////////////////////////////////////////////////////////////////////////////////////////////// Damping procedure

	/*		if (t < 90.0 * T + tau)
			{
				fftw_execute(Ck_to_Cx);
				Damping();
			}*/

			
				///////////////////////////////////////////////////////////////////////////////////////////////////////// Velocity shift

/*				if (!(nstep%(100)))
				{	

					Velocity_Correction(N_msr, tau * 100, 10);

					if (N_msr < 11)
					{
						N_msr++;
					}

					if (N_msr == 11)
					{
						N_msr = 0;
					}

				}*/

				//////////////////////////////////////////////////////////////////////////////////////////////////////// Print Data


/*				if (t > 1000.0 * T + tau)
				{
					if (!(nstep%(1 * N_av)))
					{	
					
						//Print_Ck();
						Print_Cx();
						//Print_Eta_Psi_Mu();
					
						//b_Calculate();
					}
				}*/

				if (!(nstep%(1 * N_av)))
				{
					Print_Cx();
					CheckPoint();
					Print_NPE();
					Print_Ck();
					Print_Eta_Psi_Mu();
					RV_RestartFiles();
				
				}

				if (!(nstep%(1 * N_av)))
				{

					printf("t = %e\t N = %d\n", t, nstep);
					//CheckPoint();


					//PrintAverageSurf_k();

					//PrintStatC(Amp0_Ini);
					//PrintStatSurf(Surf0_Ini);
					//BitStat();
					
					//Print_StatMoments();
				
					//ResetStat();

				}

				if (Check_Break(&Ck[0], Ak) == 1)
				{
					printf("%s\n", "Extreme amplitude wave occured, calculations stopped.");
					break;
				}

				t = t + tau;
				nstep++;
				
		}
}


















///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////// 				MAIN
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////













int main(void)
{
		clock_t t_time = clock(); // запуск таймера
		srand(time(NULL)); // Обнуление генератора функции задания случайных чисел
	
		Set_data(); // Задание всех необходимых массивов и параметров

		CreatePlans(); // Создание планов Фурье

		Clean(); // Очистка старых файлов

		fftw_execute(Ck_to_Cx);

		N_av = 100000; // Период усреднения и записи

		t = 0.0; // начальное время
		nstep = 0; // число шагов
		tau = 1e-2; // шаг по времени

		Step_Check(); // Проверка шага по времени на устойчивость

		Eta_Psi_Mu_recover_1();
		Eta_Psi_Mu_recover_2(); // Вычисление поверхности и потенциала на поверхности

		 

		double Amp0_Ini = MaxFinder(&Cx[0]); // Начальная наибольшая амплитуда 
		double Surf0_Ini = MaxFinder(&Eta_x[0]); // Начальная наибольшая амплитуда поверхности

		CheckPoint(); // Создание файла рестарта

		Print_Cx(); ///////////////
		Print_Ck(); /////////////// Печать в файлы начальных условий (используют N_av)
		Print_Eta_Psi_Mu(); //////
		RV_RestartFiles(); // запись файлов рестарта Eta(k) и Psi(k) для расчёта в RV - уравнениях
		b_Calculate(); ///////

		//Print_StatMoments();

		//B_Calculate(nstep);

		//StatMake(); // Создание функции счёта статистики

		printf(" cx_max_0 = %e\n Eta_max_0 = %e\n Mu_max_0 =  %e\n k0_dim = %e\n", Amp0_Ini, Surf0_Ini, sqrt(2) * pow(k0dim, 0.75) * Amp0_Ini/ pow(g,0.25), k0dim);

		

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Calculating_Evolution(10000.0); // Расчёт эволюции

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
		fftw_execute(Ck_to_Cx);

		Print_Cx(); /////////////////////
		Print_Ck(); ///////////////////// Печать в файл конечного состояния 
		Print_NPE(); ////////////////////
		Print_Eta_Psi_Mu(); ///////////////
		b_Calculate();
		RV_RestartFiles(); 


		DestroyPlans(); // Освобождение памяти
		Closing(); ///////

		t_time = clock() - t_time; // остановка таймера

		printf("It required %f seconds for the program to complete \n", ((float)t_time)/CLOCKS_PER_SEC);


	return 0;
}
