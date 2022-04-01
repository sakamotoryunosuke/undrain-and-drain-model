# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <cmath>
# include <time.h>
# include <iostream>
# include <fstream>
# include <iomanip>
using namespace std ;

const int Np1 = 35; //grid number

double tn = 9.*pow(10.,8.); //totla normal stress[Pa]
double p_0 = 8.995*pow(10.,8.); //initial pore fluid pressure[Pa]

double eff = 1.;//nondimensional normal stress (tn-p0)

double T = 12.4*60.*60.; //tidal period[s]

double time_max = 1000.; //simulation time
int time_step = 10;

double ad1 = 1.;
double ad = ad1*pow(10.,3.)/(tn-p_0); //tidal stress amplitude[Pa]

double omegad = 2.*M_PI;

double Sk = 0.9;// Skempton coefficient
double dd = 100000.;
double d_c = dd*pow(10.,-6.);//critical slip distance

double V_0 = pow(10.,-9.); //reference slip velocity [m/s]
double V_0d = V_0*T/d_c; //nondimensional reference slip velocity

double plate = -8.;
double V_pl = pow(10.,plate); //spring pulling velocity [m/s]
double V_pld = V_pl*T/d_c; //nondimensional spring pulling velocity [m/s]

double myu = 0.7; //sliding friction coeffecient (reference slip velocity)

double c = pow(10.,-4.);//
double dzd = 0.3;//space step size
double dtd = 0.0001;//time step size

double ak = 0.003;// direct effect coefficient
double bk = 0.002;// evolution effect coefficient

double k = pow(10.,4.);// spring stiffness
double kd = k*d_c/(tn-p_0);// nondimension spring stiffness

double sita0 = d_c/(V_pl*T);// state variable

double stemyu = myu + (ak-bk)*log(V_pld/V_0d);// friction coefficient (spring pulling velocity)

double E_p = 0.001;//dilatancy parameter
double w = 1.; //acceleration parameter (successive overrelaxation method)
double EPS=1e-8; //threshold of eeror

double velocity(double sita,double mu1)//rate and state friction law
{
	return V_pld*exp((mu1 - stemyu - bk*log(sita/sita0))/ak);
}

double theta(double Vd,double sita)// slip law
{
	return -sita*Vd*log(sita*Vd);
}

double mu(double mu1,double Vd,double t,double preP,double difP)//quasi static equation of motion
{
		double sigma = eff + (1-Sk)*ad*sin(omegad*t) - preP;//effective normal stress
	return 1./sigma*(ad*omegad*cos(omegad*t)+kd*(V_pld-Vd)+(-(1-Sk)*ad*omegad*cos(omegad*t)+difP)*mu1);
}

double SOR(double *nextP,double *P,int Z,double sita,double Vd,double *k1,double *pasP,double *difP1,double t)//pore fluid pressure change due to dilatancy/compaction
{
	double gamma = dtd/(dzd*dzd);
	double z;
	pasP[0]=P[1];
	double d,n,A,B,C,D;
	for(;;) {
	d =0.0;
	 z =log(c);
				nextP[0] = nextP[2] - 2*E_p*Vd*log(Vd*sita)*c*dzd;
		for(Z=1;Z<Np1-1;Z++){

		A = 1. - gamma/2.*exp(-z)*(exp(-(z-dzd/2.))+exp(-(z+dzd/2.)));
		B = gamma/2.*exp(-z)*exp(-(z-dzd/2.));
		C = gamma/2.*exp(-z)*exp(-(z+dzd/2.));
		D = 1. + gamma/2.*exp(-z)*(exp(-(z-dzd/2.))+exp(-(z+dzd/2.)));
		n = 1./D*(A*P[Z] + B*(nextP[Z-1]+P[Z-1])+C*(nextP[Z+1]+P[Z+1])); //Crank-Nicolson method

		n = (1.-w)*nextP[Z] + w*n;
		d+= abs(n-nextP[Z]);

		nextP[Z]=n;
		z += dzd; 

	}
	if(d<EPS){
		break;
		}
	}
		for(Z=1;Z<Np1-1;Z++){
			P[Z] = nextP[Z];
		}
				P[0] = P[2] - 2*E_p*Vd*log(Vd*sita)*c*dzd;
			k1[2]=P[1];
			difP1[0]=(P[1]-pasP[0])/dtd; //pore fluid pressure change
	return 0;
}

double adamsbashforth(double t,double dtd,double sita,double mu1,double *k1,double preP,double Vd,double difP,double *sita1,double *mu11)
{
	double kt1,km1;
	if(t==0){
	kt1 = theta(Vd,sita);
	km1 = mu(mu1,Vd,t,preP,difP);

	sita1[0] += kt1*dtd;
	mu11[0]	+= km1*dtd;

	k1[0] = kt1;
	k1[1] = km1;	
	}else if(t==dtd){
	kt1 = theta(Vd,sita);
	km1 = mu(mu1,Vd,t,preP,difP);

	sita1[0] += 0.5*dtd*(3.*kt1-k1[0]);
	mu11[0]	+= 0.5*dtd*(3.*km1-k1[1]);

	k1[3] = k1[0];
	k1[4] = k1[1];

	k1[0] = kt1;
	k1[1] = km1;
	}else{
	kt1 = theta(Vd,sita);
	km1 = mu(mu1,Vd,t,preP,difP);

	sita1[0] += dtd/12.*(23.*kt1-16.*k1[0]+5.*k1[3]);
	mu11[0]	+= dtd/12.*(23.*km1-16.*k1[1]+5.*k1[4]);

	k1[3] = k1[0];
	k1[4] = k1[1];

	k1[0] = kt1;
	k1[1] = km1;
	}
	return 0;
}


int main ()
{
	int N,Z,sum=0;
	double kt1,km1,P[Np1],Vd=V_pld,nextP[Np1],sita,k1[10],preP=0.,sita1[0],mu11[0],fai[0],mu1=stemyu,halfP[Np1],nexthalfP[Np1],difP=0.,difP1[0],pasP[0];
	sita = d_c/(V_pl*T);
	sita1[0] = d_c/(V_pl*T);
	mu11[0] = stemyu;
	difP1[0]=0.;
	pasP[0]=0.;
			for(Z=0;Z<Np1;Z++){
			P[Z] = 0.;//pore fluid pressure
			nextP[Z] = 0.;
			halfP[Z] = 0.;//pore fluid pressure
			nexthalfP[Z] = 0.;
	}
	
	char filename[128];
	ofstream file;
	sprintf(filename,"Ep_%.4f dc_%.1f plate_%.0f.dat",E_p,dd,plate);
	file.open(filename);
	clock_t start, end;
	for (double t=0.;t<=time_max+1.;t+=dtd){
			   	if(t>=time_max){
	if(sum%time_step==0){
	   	file << setprecision(3) << t-time_max  << " " << setprecision(16) << Vd/V_pld << " " << setprecision(6) << -stemyu*E_p*log(sita/sita0) << " " << setprecision(6) << ak*log(Vd/V_pld) << " " << setprecision(6) << bk*log(sita/sita0) << " " << setprecision(12) << ad*sin(omegad*t)*(1-(1-Sk)*stemyu) << endl ;
		}
sum++;
	}

		adamsbashforth(t,dtd,sita,mu1,k1,preP,Vd,difP,sita1,mu11);
		sita = sita1[0];
		mu1 = mu11[0];
		Vd = velocity(sita,mu1);
		SOR(nextP,P,Z,sita,Vd,k1,pasP,difP1,t);
		preP = k1[2];
		difP=difP1[0];
		}
	return 0;
}

