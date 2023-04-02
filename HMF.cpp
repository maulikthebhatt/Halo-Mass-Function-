#include "Parameters.H"

// USE GSL FOR THE INTEGTATION 
// TAKE SPHERICAL TOP HAT TO DEFINE SIGMA 8
// USE MORE SOPHISTICATED THAN TRAPOZODIAL 

// COMPARE THE GALAXY DENSITY WITH HMF_CALC ONLINE TOOL
// TRY TO ANSWER SOME SCIENTIFIC QUESTION USING OUR CODE

int isequal (double a, double b){
	if ( fabs(a-b) <= fabs(a*FRACT_FLOAT_ERR)){
		return 1;
	}
	else{
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////////////
// EQ. 7 from Liddle et al 1996
double OMm_z(double z){	
	return OMm*pow(1+z, 3) / (OMl + OMm*pow(1+z, 3));
}

///////////////////////////////////////////////////////////////////////////////////
// EQ. 6 from Liddle et al 1996
double g_z(double z){
	return (5/2.0) * OMm_z(z) / ((1/70.0) + (209*OMm_z(z)/140.0) - (OMm_z(z)*OMm_z(z)/140.0) + pow(OMm_z(z),0.57143));
}

///////////////////////////////////////////////////////////////////////////////////
// Function GROWTH_FACTOR returns the linear growth of matter perturbation, normalized to 1 at z=0 (from Liddle et al. 1996)

double growth_factor(double z){
	
	double OMm_z;
	
	// Einstein de Sitter 
	if (isequal(OMtot,1) && isequal(OMm,1)){
		return 1/(1.0+z);
	}
	
	// Chnage this function
	//if(isequal(OMtot,OMb+OMc+OMl) && isequal(OMtot,1)){
	if(isequal(OMtot,1)){
	// Flat LCDM eq. 7 in Liddle 
		//printf("We are in LCDM! yay!! \n");
		OMm_z = OMm * pow(1+z,3.0) / (OMl + pow(1+z,3.0)*OMm); 
		return g_z(OMm_z)/g_z(OMm)*1/(1+z);
	}
	
	fprintf(stderr, "Undefined cosmology ; -(\n");
	
	return -1;
}

///////////////////////////////////////////////////////////////////////////////////
// WINDOW FUNCTION IN FOURIER SPACE
double window_Wk(double k, double R)
{	
	if(window_function_flag == 0){
		// TOP HAT
		return 4 * M_PI * R * R * R *(sin(k*R)/k/k/k/R/R/R - cos(k*R)/k/k/R/R);
	}
	else if(window_function_flag == 1){
		// SHARP k
		if (k*R <= 1){
			return 1;
		}
		else{	
			return 0;
		}
	}
	else if(window_function_flag == 2){
		// Gaussian
		return pow(2*M_PI,1.5)* R * R * R * exp(-k*k*R*R/2);
	}
	else{
		printf("ERROR!  Window function is not defined!!");
	}
	return 0;
}

double window_volume_Vw(double R)
{	
	if(window_function_flag == 0){
		// TOP HAT 
		return 4 * M_PI * R * R * R / 3;
	}
	else if(window_function_flag == 1){
		// SHARP k
		return 6 * M_PI * M_PI * R * R * R;
	}
	else if(window_function_flag == 2){
		// Gaussian
		return pow(2*M_PI,1.5)* R * R * R;
	}
	else{
		printf("ERROR!  Window function is not defined!!");
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////
// Eq 2 in Daniel J. Eisenstein and Wayne Hu 1997 
double z_eq()
{
	return 2.5*1e4 * OMm * hlittle * hlittle / Theta_T / Theta_T / Theta_T / Theta_T ;
}

/////////////////////FOR THE COMPLEX TRANSFER FUNCTION ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// Eq 2 in Daniel J. Eisenstein and Wayne Hu 1997 
double b1()
{
	return 0.313 * pow((OMm*hlittle*hlittle),-0.419) * (1 + 0.607* pow((OMm*hlittle*hlittle),0.674)) ;
}

double b2()
{
	return 0.238 * pow((OMm*hlittle*hlittle),0.223);
}

double z_d()
{
	return 1291 * pow((OMm*hlittle*hlittle),0.251) * (1 + b1() * pow((OMm*hlittle*hlittle),b2())) / (1 + 0.689 * pow((OMm*hlittle*hlittle),0.828))   ;
}

//////////////////////////////////////////////////////////////////////////////////
// Eq 3 in Daniel J. Eisenstein and Wayne Hu 1997 
double y_d()
{
	return (1 + z_eq()) / (1 + z_d());
}

//////////////////////////////////////////////////////////////////////////////////
// Eq 4 in Daniel J. Eisenstein and Wayne Hu 1997 
double sound_horizon()
{
	// Sound horizone 
	return (44.5 * log(9.83/(OMm*hlittle*hlittle))) / sqrt(1+10*pow((OMb*hlittle*hlittle),0.75));
}

//////////////////////////////////////////////////////////////////////////////////
// Eq 5 in Daniel J. Eisenstein and Wayne Hu 1997 
double q_(double k)
{
	//return (k / 19.0) * pow((OMm*100*hlittle*100*hlittle),-0.5) * pow((1+z_eq()),-0.5); //
	return k * Theta_T * Theta_T / (OMm*hlittle*hlittle); 
}

//////////////////////////////////////////////////////////////////////////////////
// Eq 9 in Daniel J. Eisenstein and Wayne Hu 1997 
double g2_z(double z)
{
	return OMm * pow((1+z),3) + (1-OMm-OMl)*pow((1+z),2) + OMl;
}

/////////////////////////////////////////////////////////////////////////////////
// Eq 10 in Daniel J. Eisenstein and Wayne Hu 1997 

// Here OmegaLamda(z) is defind as 1-OmegaMatter(z) (flat LCDM) 
double OM_z(double z)
{
	return OMm * pow((1+z),3) /g2_z(z);
}
double OMl_z(double z){
	return OMl /g2_z(z);
}
double D1_z(double z){
	return	((1 + z_eq())/(1+z)) * (5 * OM_z(z)/2.0) / ((pow(OM_z(z),0.57143) - OMl_z(z)) + ((1 + OM_z(z)/2.0)*(1 + OMl_z(z)/70.0))) ;
} 

////////////////////////////////////////////////////////////////////////////////
// Eq 11 in Daniel J. Eisenstein and Wayne Hu 1997 
double p_factor(double f)
{
	return (5 - sqrt(1+24*f))/4;
}

////////////////////////////////////////////////////////////////////////////////
// Eq 14 in Daniel J. Eisenstein and Wayne Hu 1997 
double epoch_fs(double q)
{
	return 17.2 * f_nu * ( 1 + 0.488 * pow(f_nu,-1.16667)) * (pow(N_nu*q/f_nu,2));
}

// Eq 12 in Daniel J. Eisenstein and Wayne Hu 1997 
double D_cb(double q, double z)
{
	return pow((1 + pow((D1_z(z)/(1+epoch_fs(q))),0.7)),(p_factor(f_cb)/0.7)) * pow(D1_z(z),(1-p_factor(f_cb)));
}

// Eq 13 in Daniel J. Eisenstein and Wayne Hu 1997 
double D_cbnu(double q,double z)
{
	return pow((pow(f_cb,(0.7/p_factor(f_cb)))+ pow(D1_z(z)/(1+epoch_fs(q)),0.7)),(p_factor(f_cb)/0.7)) * pow(D1_z(z),(1-p_factor(f_cb)));
}

// Eq 15 in Daniel J. Eisenstein and Wayne Hu 1997 
double alpha_nu()
{
	double term1, term2, term3;
	term1 = f_c * (5 - 2 * (p_factor(f_c) + p_factor(f_cb))) / f_cb / (5 - 4 *p_factor(f_cb));
	term2 = (1 - 0.553 * f_nub + 0.126 * f_nub * f_nub * f_nub ) * pow((1 + y_d()),(p_factor(f_cb)-p_factor(f_c))) / (1 - 0.193 * sqrt(f_nu*N_nu) + 0.169 * f_nu * pow(N_nu,0.2));
	term3 = 1 + ((p_factor(f_c) - p_factor(f_cb))/2) * (1 + (1/(3-4*p_factor(f_c))/(7-4*p_factor(f_cb))))/(1+y_d()) ;
	return  term1 * term2 * term3;
}

// Eq 16 in Daniel J. Eisenstein and Wayne Hu 1997 
double Gamma_eff(double k)
{
	return OMm * hlittle * hlittle * (sqrt(alpha_nu()) + (1 - sqrt(alpha_nu()))/pow((1+(0.43*k*sound_horizon()/hlittle)),4));
}

// Eq 17 in Daniel J. Eisenstein and Wayne Hu 1997 
double q_eff(double k)
{
	return  k * hlittle * pow(Theta_T,2) / Gamma_eff(k);   // we have extra hlittle compared to the exact Eq17 : this is because the "k" already had "h" term.
}

// Eq 21 in Daniel J. Eisenstein and Wayne Hu 1997 
double beta_c(double k)
{
	return 1 / (1 - (0.949 * f_nub));
}

// Eq 20 in Daniel J. Eisenstein and Wayne Hu 1997 
double C_(double k)
{	
	return 14.4 + ( 325 / (1+60.5*pow(q_eff(k),1.08)) );
}

// Eq 19 in Daniel J. Eisenstein and Wayne Hu 1997 
double L_(double k)
{
	return log(2.718 + 1.84 * beta_c(k) * sqrt(alpha_nu()) * q_eff(k));
}

// Eq 18 in Daniel J. Eisenstein and Wayne Hu 1997 
double Tk_sup(double k)
{
	return L_(k) / (L_(k) + (C_(k) * q_eff(k) * q_eff(k)));
}

// Eq 23 in Daniel J. Eisenstein and Wayne Hu 1997 
double q_nu(double k)
{	
	//return k / (3.42 * sqrt(f_nu/N_nu)) / 10000;
	return  3.92 * q_(k) * sqrt(N_nu/f_nu) ;
}

// Eq 22 in Daniel J. Eisenstein and Wayne Hu 1997 
double correction_B(double k)
{	
	return 1 + ( 1.24 * pow(f_nu,0.64) * pow(N_nu,(0.3+0.6*f_nu)) / (pow(q_nu(k),-1.6) + pow(q_nu(k),0.8)) );
}

// Eq 24 in Daniel J. Eisenstein and Wayne Hu 1997 
double Tk_master(double k)
{
	return Tk_sup(k) * correction_B(k);
}

// Eq 6 in Daniel J. Eisenstein and Wayne Hu 1997 
double Tk_cb(double k, double z)
{	
	//double q_ = k / 19 / pow(OMm*hlittle*hlittle,0.5) / pow(1+z_eq(),0.5) ;
	return Tk_master(k) * D_cb(k,z) / D1_z(z) ;
}

// Eq 7 in Daniel J. Eisenstein and Wayne Hu 1997 
double Tk_cbnu(double k, double z)
{	
	double q_ = k * Theta_T * Theta_T / (OMm*hlittle*hlittle); // * 19 * pow(OMm*100*100*hlittle*hlittle,0.5) * pow(1+z_eq(),0.5) ;
	//double q_ = k / 19 / pow(OMm*100*100*hlittle*hlittle,0.5) / pow(1+z_eq(),0.5) ;
	return Tk_master(k) * D_cbnu(q_,z) / D1_z(z);
}

//////////////////////////////////////////////////////////////////////////////////
// Transfer function
double Tk(double k, double z)
{	
	if (Tk_flag == 0)
	{
		double q_bardeen = k / Gamma_bardeen;
		return  log(1.0+2.34*q_bardeen) / (2.34 * q_bardeen * pow((1.0+3.89*q_bardeen+pow(16.1*q_bardeen,2.0)+pow(5.46*q_bardeen,3.0)+pow(6.71*q_bardeen,4.0)),0.25));
	}
	else if (Tk_flag == 1)
	{
		return Tk_cbnu(k,z) ;
	}
	else{
		cout << "ERROR! Tk flag is not corretly set up!" << endl;
		return 0;
	}	
}

//////////////////////////////////////////////////////////////////////////////////
// Eq 109 from class notes
double Pk(double k, double z, double norm_)
{	
	if (Tk_flag == 0)
	{	
		return norm_ * pow(k,POWER_INDEX) * Tk(k,z) * Tk(k,z) * growth_factor(z) * growth_factor(z) / growth_factor(0) / growth_factor(0); //* g_z(z) * g_z(z) / g_z(0) / g_z(0) ;
	}
	else if (Tk_flag == 1)
	{
		return norm_ * pow(k,POWER_INDEX) * Tk_cbnu(k,z) * Tk_cbnu(k,z) * D1_z(z) * D1_z(z) / D1_z(0) / D1_z(0) ;
	}
	else{
		cout << "ERROR! Tk flag is not corretly set up!" << endl;
		return 0;
	}	
}


////////////////////////////////////////////////////////////////////////////////////
////// Sigma M   eq 104 from class notes   ( Function for integratoin )
double sigma_M_sq (double k, void * params) {
	
	double R = ((double *)params)[0];
	double z = ((double *)params)[1];
	double norm_ = ((double *)params)[2];
	
	double sigma_M_sq = window_Wk(k,R) * window_Wk(k,R) * k * k * Pk(k,z,norm_) / window_volume_Vw(R)  / window_volume_Vw(R) / 2 / M_PI / M_PI;
	return sigma_M_sq;
}

int main(int argc, char ** argv){
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000000);
	
	double z, result, error, norm_, k, k_end, R_start, R_end, R, R_old, R_old2, fcoll, fcoll_old, fcoll_old2, Mass, Mass_old, dndlogm;
	
	//k_end = 100/hlittle;
	//k = 0.0001/hlittle;	
	k_end = 99.485;   
	k = 0.0001234098;     
	z = 0;
	R = 8.0;
	norm_ = 1.0;
	
	cout << "COSMOLOGICAL PARAMETERS USED! "<< endl;
	cout << "Omega_matter = " << OMm << endl;
	cout << "Omega_baryon = " << OMb << endl;
	cout << "Omega_lamda  = " << OMl << endl;
	cout << "Omega_nu     = " << OMnu << endl;
	cout << "Power index  = " << POWER_INDEX << endl;
	cout << "Sigma 8      = " << SIGMA8 << endl;
	cout << "T_CMB        = " << T_CMB << endl;
	cout << "N_nu         = " << N_nu << endl;
	cout << "h            = " << hlittle << endl << endl;
	
	cout << "Baryons are released from the Compton drag of the photons at red shift z_d = " << z_d() << endl;
	cout << "Matter Radiation equality red shift z_eq = " << z_eq() << endl << endl;
	
	double params[] = {R,z,norm_};	
	
	gsl_function G;
	G.function = &sigma_M_sq;
	G.params = &params;
	
	gsl_integration_qagiu (&G, 0, 0, 1e-8, 100000000, w, &result, &error);
	
	norm_ = SIGMA8*SIGMA8/result;	
	cout << "Default top hat window function is used to get the normalization !!" << endl;
	cout << "NORMALIZATION constant = " << norm_ << endl;
	
	ofstream Tkfile;
	ofstream Pkfile;
	ofstream sigmaM2file;
	ofstream fcolfile;
	ofstream dndlogmfile;
	
	z = 0;
	
	cout << "\nReading Flags!! "<<endl;
	if (window_function_flag == 0 ){
		cout << "Window function flag is set to 0 : Top hat window function will be used!" << endl; 
	}
	else if ( window_function_flag == 1 ){
		cout << "Window function flag is set to 1 : Sharp k window function will be used!" << endl;
	}
	else if ( window_function_flag == 2){
		cout << "Window function flag is set to 2 : Gaussian window function will be used!" << endl;
	}
	else{
		cout << "ERROR!  Window function is not defined!!" << endl;
	}
	
	if (Tk_flag == 0){
		Tkfile.open ("Tk0_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		Pkfile.open ("Pk0_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		sigmaM2file.open ("SigmaM0_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		fcolfile.open ("fcol0_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		dndlogmfile.open ("dndlogm0_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		cout << "Tk flag is set to 0 : Liddle 1995 transfer function will be used for the HMF!" << endl;
	}
	else if (Tk_flag == 1){
		Tkfile.open ("Tk1_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		Pkfile.open ("Pk1_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		sigmaM2file.open ("SigmaM1_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		fcolfile.open ("fcol1_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		dndlogmfile.open ("dndlogm1_z"+to_string(((int)(z*100))*0.01)+".txt",ios::trunc);
		cout << "Tk flag is set to 1 : Daniel J. Eisenstein and Wayne Hu 1997 transfer function will be used for the HMF!" << endl;
	}
	else{
		cout << "ERROR! Tk flag is not corretly set up!" << endl;
	}
	
	cout << "\nRedshift z = " << z << endl;
	Tkfile << "#k" << " " << "#Tk" << " "<< "#Tk_nume" << endl;
	Pkfile << "#k" << " " << "#Pk" << endl;
	sigmaM2file << "#M" << " " << "#SigmaM" << endl;
	fcolfile << "#R" << " " << "#fcol" << endl;
	dndlogmfile << "#M" << " " << "#dn/dlogm" << endl;
	
	R_start = 0.001;
	R_end   = 100;
	R       = R_start;
	R_old   = 0;
	R_old2  = 0;
	
	//cout << "growth factor z = " << growth_factor(z) << endl;
	//cout << " growth factor z0 = " << growth_factor(0) << endl;
	
	while(R<R_end){
		double params[] = {R,z,norm_};	
		G.params = &params;
		gsl_integration_qagiu (&G, 0, 0, 1e-8, 100000000, w, &result, &error);
		
		fcoll = erfc(1.686/growth_factor(z)/sqrt(2.0*result));
		
		fcolfile << R << " " << fcoll << endl;
		Mass = 4*M_PI*R*R*R*rho0/3;
		
		sigmaM2file << Mass << " " << result << endl;
		
		//if (R>R_start){
		//	dndlogm = abs((fcoll - fcoll_old)/ (R - R_old)) / R / R / 4 / M_PI;
		//	dndlogmfile << (Mass+Mass_old)/2 << " " << dndlogm << endl;
		//}
		if (R_old2>0){
			dndlogm = abs((fcoll - fcoll_old2)/ (R - R_old2)) / R_old / R_old / 4 / M_PI;
			dndlogmfile << Mass_old << " " << dndlogm << endl;
		}
		R_old2 = R_old;
		R_old  = R;
		
		fcoll_old2 = fcoll_old;
		fcoll_old = fcoll;
		
		Mass_old = Mass;
		
		R = R * 1.0005;
	}
	while(k<k_end){
		Tkfile << k << " " << Tk(k,z) << endl;
		Pkfile << k << " " << Pk(k,z,norm_) << endl;
		//k = k * 1.05;
		//k = exp(log(k)+ 0.05);
		k = exp(log(k)+ 0.01);
	}
	
	Tkfile.close();
	Pkfile.close();
	sigmaM2file.close();
	fcolfile.close();
	dndlogmfile.close();
	cout << "f nu = " << f_nu << endl;
	cout << "Total Omega = " << OMb + OMc + OMl + OMnu << endl;
	return 0;
}
