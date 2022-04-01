# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
using namespace std ;

double tn = 9.*pow(10.,8.); //totla normal stress[Pa]
double p_0 = 8.995*pow(10.,8.); //initial pore fluid pressure[Pa]

double eff = 1.;//nondimensional normal stress (tn-p0)

double T = 12.4*60.*60.; //tidal period[s]

double time_max = 1000.; //simulation time
int time_step = 1;

double ad1 = 1.;
double ad = ad1*pow(10.,3.)/(tn-p_0); //tidal stress amplitude[Pa]

double omegad = 2.*M_PI;

double Sk = 0.9;// Skmpton coefficient
double dd = 100000.;
double d_c = dd*pow(10.,-6.);//critical slip distance

double V_0 = pow(10.,-9.); //reference slip velocity [m/s]
double V_0d = V_0*T/d_c; //nondimensional reference slip velocity

double plate = -8.;
double V_pl = pow(10.,plate); //spring pulling velocity [m/s]
double V_pld = V_pl*T/d_c; //nondimensional spring pulling velocity [m/s]

double myu = 0.7; //sliding friction coeffecient (reference slip velocity)

double dtd = 0.001;// time step size

double ak = 0.003;// direct effect coefficient
double bk = 0.002;// evolution effect coefficient

double k = pow(10.,4.);// spring stiffness
double kd = k*d_c/(tn-p_0);// nondimension spring stiffness

double sita0 = d_c/(V_pl*T);// state variable

double stemyu = myu + (ak-bk)*log(V_pld/V_0d);// friction coefficient (spring pulling velocity)

double U = 1.;//dilarancy parameter

double velocity(double sita,double mu1)//rate and state friction law
{
	return V_pld*exp((mu1 - stemyu - bk*log(sita/sita0))/ak);
}

double theta(double Vd,double sita)//slip law
{
	return -sita*Vd*log(sita*Vd);
}

double mu(double mu1,double Vd,double t,double preP,double difP,double sita)//quasi static equation of motion
{
		double sigma = eff + (1.-Sk)*ad*sin(omegad*t) - preP; //effective normal stress
	return 1./sigma*(ad*omegad*cos(omegad*t)+kd*(V_pld-Vd)+(-(1.-Sk)*ad*omegad*cos(omegad*t)-U*Vd*log(sita*Vd))*mu1);
}

double dilatancy(double sita,double Vd)//pore fluid pressure change due to dilatancy/compaction
{
	return -U*Vd*log(sita*Vd);
}

double RungeKutta(double t,double dtd,double sita,double mu1,double preP,double Vd,double difP,double *sita1,double *mu11,double *P)
{
	double kt1,km1,kp1,kt2,km2,kp2,kt3,km3,kp3,kt4,km4,kp4,difnew;
	kt1 = theta(Vd,sita)*dtd;
	km1 = mu(mu1,Vd,t,preP,difP,sita)*dtd;
	kp1 = dilatancy(sita,Vd)*dtd;

	Vd = velocity(sita+0.5*kt1,mu1+0.5*km1);

	kt2 = theta(Vd,sita+0.5*kt1)*dtd;
	km2 = mu(mu1+0.5*km1,Vd,t+0.5*dtd,preP+kp1*0.5,difP,sita+kt1*0.5)*dtd;
	kp2 = dilatancy(sita+kt1*0.5,Vd)*dtd;

	Vd = velocity(sita+0.5*kt2,mu1+0.5*km2);

	kt3 = theta(Vd,sita+0.5*kt2)*dtd;
	km3 = mu(mu1+0.5*km2,Vd,t+0.5*dtd,preP+kp2*0.5,difP,sita+kt2*0.5)*dtd;
	kp3 = dilatancy(sita+0.5*kt2,Vd)*dtd;

	Vd = velocity(sita+0.5*kt3,mu1+0.5*km3);

	kt4 = theta(Vd,sita+kt3)*dtd;
	km4 = mu(mu1+km3,Vd,t+dtd,preP+kp3,difP,sita+kt3)*dtd;
	kp4 = dilatancy(sita,Vd)*dtd;

	sita1[0] += (kt1+2.*kt2+2.*kt3+kt4)/6.;
	mu11[0]	+=  (km1+2.*km2+2.*km3+km4)/6.;
	P[0] += (kp1+2.*kp2+2.*kp3+kp4)/6.;
	return 0;
}

int main ()
{
	int N,Z,sum=0;
	double P[0],Vd=V_pld,sita,preP=0.,sita1[0],mu11[0],mu1=stemyu,difP=0;
	sita = d_c/(V_pl*T);
	sita1[0] = d_c/(V_pl*T);
	mu11[0] = stemyu;
	P[0]=0.;
	
	char filename[128];
	ofstream file;
	sprintf(filename,"U_%.2f dc_%.1f plate_%.0f.dat",U,dd,plate);
	file.open(filename);
	for (double t=0.;t<=time_max+1;t+=dtd){
			   	if(t>=time_max){
	if(sum%time_step==0){
	   	file << setprecision(3) << t-time_max  << " " << setprecision(16) << Vd/V_pld << " " << setprecision(6) << -stemyu*U*log(sita/sita0) << " " << setprecision(6) << ak*log(Vd/V_pld) << " " << setprecision(6) << bk*log(sita/sita0) << " " << setprecision(12) << ad*sin(omegad*t)*(1-(1-Sk)*stemyu) << " " << preP << endl ;
		}
sum++;
	}
		RungeKutta(t,dtd,sita,mu1,preP,Vd,difP,sita1,mu11,P);
		sita = sita1[0];
		mu1 = mu11[0];
		preP = P[0];
		Vd = velocity(sita,mu1);
		}
	return 0;
}

