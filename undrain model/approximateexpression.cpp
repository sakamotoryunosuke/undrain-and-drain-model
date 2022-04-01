# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <cmath>
# include <iostream>
# include <fstream>
# include <complex>
using namespace std ;

double tn = 9.*pow(10.,8.); //totla normal stress[Pa]
double p_0 = 8.995*pow(10.,8.); //initial pore fluid pressure[Pa]

double eff = 1.;//nondimensional normal stress (tn-p0)

double T = 12.4*60.*60.; //tidal period[s]

double ad1 = 0.02;
double ad = ad1*pow(10.,3.)/(tn-p_0); //tidal stress amplitude[Pa]

double steady = -9.;
double V_0 = pow(10.,steady); //reference slip velocity [m/s]

double plate = -8.;
double V_pl = pow(10.,plate); //pulling spring velocity [m/s]

double yuna = plate - steady;

double ak = 0.003;// direct effect coefficient
double bk = 0.002;// evolution effect coefficient

double myu = 0.7; //sliding friction coeffecient (reference slip velocity)

double stemyu = myu + (ak-bk)*log(V_pl/V_0);// friction coefficient (spring pulling velocity)

double k = pow(10.,4.);// spring stiffness

double U = 1.;//dilatancy parameter

int main ()
{
	int n,sum;
	char filename[128];
	ofstream file;
	sprintf(filename,"approximatesolution U_%.2f plate_%.0f.dat",U,plate);
	file.open(filename);
for(double T_d= pow(10.,-2.) ;T_d<=pow(10.,-1.);T_d+=pow(10.,-3.)){//T_d=T_θ/T
	
double kd = k/(tn-p_0)*T_d*T/(2*M_PI)*V_pl;//nondimensional spring stiffness
double V_pld = 2*M_PI/T_d;// nondimensional pulling spring velocity

	std::complex<double> D(1.,T_d);//1+i*T_θ/T
	std::complex<double> C(1./sqrt(2.),1./sqrt(2.)); //√i
	std::complex<double> i(0.,1.);//i
	
	double Breal = std::real(ak-1./D*(bk-stemyu*U));
	double Bimag = std::imag(ak-1./D*(bk-stemyu*U));
	
	std::complex<double> BB(Breal,Bimag);//eq(14)
	
	double HDreal   = std::real(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	double HDimag   = std::imag(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	
	std::complex<double> HD(HDreal,HDimag);
	
	double result1 = std::arg(HD); //phase difference
	double result2 = std::real(HD)/((tn-p_0)*0.001); //tidal sensitivity

 	file  << T_d << " " << result1 << " " << result2 << endl ;
 
	}
	for(double T_d= pow(10.,-1.) ;T_d<=pow(10.,0.);T_d+=pow(10.,-2.)){
	
double kd = k/(tn-p_0)*T_d*T/(2*M_PI)*V_pl;//nondimensional spring stiffness
double V_pld = 2*M_PI/T_d;// nondimensional pulling spring velocity

	std::complex<double> D(1.,T_d);//1+i*T_θ/T
	std::complex<double> C(1./sqrt(2.),1./sqrt(2.)); //√i
	std::complex<double> i(0.,1.);//i
	
	double Breal = std::real(ak-1./D*(bk-stemyu*U));
	double Bimag = std::imag(ak-1./D*(bk-stemyu*U));
	
	std::complex<double> BB(Breal,Bimag);//eq(14)
	
	double HDreal   = std::real(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	double HDimag   = std::imag(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	
	std::complex<double> HD(HDreal,HDimag);
	
	double result1 = std::arg(HD); //phase difference
	double result2 = std::real(HD)/((tn-p_0)*0.001); //tidal sensitivity

 	file  << T_d << " " << result1 << " " << result2 << endl ;
 
	}
	
		for(double T_d= pow(10.,0.) ;T_d<=pow(10.,1.);T_d+=pow(10.,-1.)){

double kd = k/(tn-p_0)*T_d*T/(2*M_PI)*V_pl;//nondimensional spring stiffness
double V_pld = 2*M_PI/T_d;// nondimensional pulling spring velocity

	std::complex<double> D(1.,T_d);//1+i*T_θ/T
	std::complex<double> C(1./sqrt(2.),1./sqrt(2.)); //√i
	std::complex<double> i(0.,1.);//i
	
	double Breal = std::real(ak-1./D*(bk-stemyu*U));
	double Bimag = std::imag(ak-1./D*(bk-stemyu*U));
	
	std::complex<double> BB(Breal,Bimag);//eq(14)
	
	double HDreal   = std::real(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	double HDimag   = std::imag(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	
	std::complex<double> HD(HDreal,HDimag);
	
	double result1 = std::arg(HD); //phase difference
	double result2 = std::real(HD)/((tn-p_0)*0.001); //tidal sensitivity

 	file  << T_d << " " << result1 << " " << result2 << endl ;
 
	}
	
		for(double T_d= pow(10.,1.) ;T_d<=pow(10.,2.);T_d+=pow(10.,0.)){
	
double kd = k/(tn-p_0)*T_d*T/(2*M_PI)*V_pl;//nondimensional spring stiffness
double V_pld = 2*M_PI/T_d;// nondimensional pulling spring velocity

	std::complex<double> D(1.,T_d);//1+i*T_θ/T
	std::complex<double> C(1./sqrt(2.),1./sqrt(2.)); //√i
	std::complex<double> i(0.,1.);//i
	
	double Breal = std::real(ak-1./D*(bk-stemyu*U));
	double Bimag = std::imag(ak-1./D*(bk-stemyu*U));
	
	std::complex<double> BB(Breal,Bimag);//eq(14)
	
	double HDreal   = std::real(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	double HDimag   = std::imag(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	
	std::complex<double> HD(HDreal,HDimag);
	
	double result1 = std::arg(HD); //phase difference
	double result2 = std::real(HD)/((tn-p_0)*0.001); //tidal sensitivity

 	file  << T_d << " " << result1 << " " << result2 << endl ;
 
	}
	
		for(double T_d= pow(10.,2.) ;T_d<=pow(10.,3.);T_d+=pow(10.,1.)){
	
double kd = k/(tn-p_0)*T_d*T/(2*M_PI)*V_pl;//nondimensional spring stiffness
double V_pld = 2*M_PI/T_d;// nondimensional pulling spring velocity

	std::complex<double> D(1.,T_d);//1+i*T_θ/T
	std::complex<double> C(1./sqrt(2.),1./sqrt(2.)); //√i
	std::complex<double> i(0.,1.);//i
	
	double Breal = std::real(ak-1./D*(bk-stemyu*U));
	double Bimag = std::imag(ak-1./D*(bk-stemyu*U));
	
	std::complex<double> BB(Breal,Bimag);//eq(14)
	
	double HDreal   = std::real(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	double HDimag   = std::imag(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	
	std::complex<double> HD(HDreal,HDimag);
	
	double result1 = std::arg(HD); //phase difference
	double result2 = std::real(HD)/((tn-p_0)*0.001); //tidal sensitivity

 	file  << T_d << " " << result1 << " " << result2 << endl ;
 
	}

			for(double T_d= pow(10.,3.) ;T_d<=pow(10.,4);T_d+=pow(10.,2.)){
	
double kd = k/(tn-p_0)*T_d*T/(2*M_PI)*V_pl;//nondimensional spring stiffness
double V_pld = 2*M_PI/T_d;// nondimensional pulling spring velocity

	std::complex<double> D(1.,T_d);//1+i*T_θ/T
	std::complex<double> C(1./sqrt(2.),1./sqrt(2.)); //√i
	std::complex<double> i(0.,1.);//i
	
	double Breal = std::real(ak-1./D*(bk-stemyu*U));
	double Bimag = std::imag(ak-1./D*(bk-stemyu*U));
	
	std::complex<double> BB(Breal,Bimag);//eq(14)
	
	double HDreal   = std::real(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	double HDimag   = std::imag(2.*M_PI*i/(kd*V_pld+2.*M_PI*i*BB));
	
	std::complex<double> HD(HDreal,HDimag);
	
	double result1 = std::arg(HD); //phase difference
	double result2 = std::real(HD)/((tn-p_0)*0.001); //tidal sensitivity

 	file  << T_d << " " << result1 << " " << result2 << endl ;
 
	}
	
	return 0;
}