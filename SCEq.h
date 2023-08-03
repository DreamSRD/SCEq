
	#define M_PI 3.1415926535

	#define N 8192
	#define Nstat 2000
	
	int k0 = 100;
	double a = 0.0;
	double b = 10000;
	double x0 = 0 * M_PI;

	double V = 6.24920; //k0 = 100
	//double V = 9.878; //k0 = 50
	//double V = 1.9757;

	double V_br = 0.0; // Скорость бризера
	double Coordinate[21];
	double Velocity[21];

	double g = 9.81;
	double Dk,h,t,tau,w0, k0dim, T0, u_lim1, u_lim2, l_lim1, l_lim2, hstat1, hstat2,x1,x2, T;
	int nstep = 0;
	int N_av = 0;
	char Name[100];
/*	double CR[Nstat + 1];
	double SurfR[Nstat + 1];
	unsigned long long  HSurf_int[Nstat + 2];
	unsigned long long  HC_int[Nstat + 2];
	unsigned long long int HCabs_int[Nstat + 2];
	unsigned long long int HCabs2_int[Nstat + 2];
	double HSurf_double[Nstat + 2];
	double HC_double[Nstat + 1];
	double HCabs_double[Nstat + 2];
	double HCabs2_double[Nstat + 2];*/
	double Eta_k_average[N/2];

	int Max_int = 0;

	int number_of_numbers;


	fftw_complex *Bx, *Bk, *Cx, *Ck, *CkOld, *CkNew, *CCx, *CCk, *Cdx, *Cdk, *CCCdx, *CCCdk, *KCCk, *KCCx, *KCCCx, *KCCCk;
	fftw_complex *RHSk1, *RHSk2, *RHSk3, *RHSk4, *K, *Eta_x, *Psi_x, *Psi_k, *Eta_2_x, *Eta_2_k, *Eta_k, *Mu_x, *Mu_k, *Mu_2_k, *Mu_2_x;
	fftw_complex *Psi_2m_k, *Psi_2p_k, *Psi_2m_x, *Psi_2p_x, *Psi_2_x, *Psi_2_k;
	
	fftw_plan Cx_to_Ck, Ck_to_Cx, CCx_2_CCk, Cdk_to_Cdx, KCCk_to_KCCx, CCCdx_2_CCCdk, KCCCx_2_KCCCk, Bk_to_Bx;
	fftw_plan Eta_k_to_Eta_x, Eta_2_x_to_Eta_2_k, Eta_2_k_to_Eta_2_x, Mu_k_to_Mu_x, Mu_2_k_to_Mu_2_x;
	fftw_plan Psi_k_to_Psi_x, Psi_2m_k_to_Psi_2m_x, Psi_2p_k_to_Psi_2p_x, Psi_2_x_to_Psi_2_k,  Psi_2_k_to_Psi_2_x;

	FILE *pfile;