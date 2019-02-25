/* 
One of 3 scripts used to model electroporation as discussed in the article: 
Curvature-Driven Pore Growth in Charged Membranes during Charge-Pulse and Voltage-Clamp Experiments
Biophys J. 2009 Feb 4; 96(3): 907–916.
doi: 10.1016/j.bpj.2008.10.035
Article authors: Jens Kroeger, Dan Vernon, Martin Grant
Script author: Jens Kroeger
Revision date: February 21, 2019
Free for use, copy and modify with proper citation.
*/




//		Import libraries

// get gaussian random number generator (random_n)
#include "random.i"

//		Define functions


//		Declare constants and variables with corresponding units


//Declaration of main variables
palette, "gray.gp"
random_seed,.1234567
k_BT=0.5   // Tc is around 0.5 apparently
MCS=1000
time=20000;
radius=100;
fixedradius=40.0;
DEL_X=1.0/sqrt(2.0);
center=int(radius/2);
rad=array(0.5,[1,radius]);
densityofholes=array(0.0,[2,radius+20,radius+20]);
K=array(0.0,[2,radius+20,radius+20]);
equilibriumdensity=array(0.0,[2,radius+20,radius+20]);
holesatrimevolution=array(0.0,[1,time]);
voltageevolution=array(0.0,[1,time]);
equilibriumdensityevolution=array(0.0,[1,time]);
voidevolution=array(0.0,[1,time]);
domain=array(0.0,[2,radius+20,radius+20]);
holesaverageevolution=array(0.0,[1,time]);
equilibriumdensity=array(1.0,[1,radius]);


for(i=1;i<=radius;i++){
	rad(i)=abs(i-center)*DEL_X;
	}
kappas=10.0;//de facon a ce que f(0)\sim 10 dyne/cm en accord avec Kuzmin et al.//les unites de kappas sont en kT/nm^2.
epsilon1=0.02;
epsilon2=0.009;
deltat=0.001;
n0=1.0; //d'apres Debruin
hsquared=5.0;//5.0 nm^2 d'apres Kuzmin et al.
lambdasquared=8.0;//en accord avec Kuzmin et al. (2005) si \lambda^2<h_s^2
splay=-0.1;// units are 1/nm from Kozlovsky (2002)
beta=0.502027//in terms of cm/dynes at room temperature.
conductanceperdensity=3.77e-10;//in pScm^2;
q=2.46; // d'apres Debruin
alpha=100.0;//units in cm^{-1}ms^{-1} d'apres Debruin
N0=1.5e5;// d'apres Debruin
beta=0.502027//in terms of cm/dynes at room temperature.
equilibriumelectroporationvoltage=270.0;
equilibriumcurrentvoltage=45.0; //en mV
calciumout=0.55; // in mu M;
//C_o min=C_i*exp(equilibriumcurrentvoltage/25.8)=0.0575 for
// C_i=0.01 and equilibriumcurrentvoltage=45.0;
calciumin=0.01; //in \mu M;
voltage=25.8*log(calciumout/calciumin);

// 		Load data set

flr=open("domain3.txt");
read,flr,domain;
center=int((radius+20)/2)
sxy=radius+10;
bias=1.1;

m=array(0.0,[2,radius+20,radius+20]);
for(l=1;l<=radius+10;l++)
	for(k=1;k<=radius+10;k++){
		if((1.0*l-center)^2.0+(1.0*k-center)^2.0<=9.0){m(l,k)=-1.0;}};
ftm=fft(roll(m))/10.0/0.77;//works best for ftg=fft(m)/cutoff/0.5
psi=1.0+2.0*roll(m,[0,38])+2.0*roll(m,[1,38]); 
DEL_T=.1 
k_BT=0.000002   // Tc is around 0.5 apparently
MCS=1000
noise=sqrt(2.0*k_BT/(DEL_X*DEL_X*DEL_T))

summ=0.0;
for(i=1;i<=radius+20;i++) 
	for(j=1;j<=radius+20;j++) {
		summ=summ+(1.0-psi(i,j));}

	fma;pli,psi;


//		Loop over time for a given set of parameters
// Beginning of time loop/iteration over time using simple Newton algorithmn.

for(t=1;t<=time;t++) {

	calciumin=calciumin+10.0*DEL_T*conductanceperdensity*summ*(voltage);
	//calciumin=calciumin+10.0*DEL_T*conductanceperdensity*summ*(calciumout-calciumin);


	voltage=25.8*log(calciumout/calciumin);

	summ=0.0;
	for(i=1;i<=radius+20;i++) 
		for(j=1;j<=radius+20;j++) {
			summ=summ+(1.0-psi(i,j))/2.0;}


	K=alpha*exp((voltage/equilibriumelectroporationvoltage)^2.0+beta*domain);

	equilibriumdensity=N0*exp(q*(voltage/equilibriumelectroporationvoltage)^2.0+beta*domain);

	densityofholes=densityofholes+deltat*(K-K*densityofholes/equilibriumdensity);

	voltageevolution(t-1)=voltage;
	voidevolution(t-1)=summ;

	psi += DEL_T*(-0.5*((psi+1.0)/2.0)*(voltage/equilibriumelectroporationvoltage)^2.0-0.035*((psi+1.0)/2.0)*(domain)+psi-psi*psi*psi+0.3*laplacian(psi))

		if(k_BT>1e-6) { // speed up for T=0
			psi += DEL_T*noise*(random_n([2,radius+20,radius+20])+1.0-0.0045*exp(1.05*domain+1.5*(voltage/equilibriumelectroporationvoltage)^2.0))
		}
	

	write,format="%10d\r",t
	//pause,100

	if(i%10==0 && i>100000){
		ftz=fft(random_n([2,radius+20,radius+20])+0.15);
		prod=ftm*ftz/sxy/sxy;
		convolution=fft(prod,-1);
		psi += deltat*domain*re_part(convolution)/100.0;
	}	
	
	fma;pli,psi;
	write,format="%10d\r",t
	
}
write,format="%s\n","done"
animate,0


//		Plot results


window, 0;
palette, "yarg.gp";
limits;
fma;pli, psi;

window, 2;
palette, "yarg.gp";
limits;
fma;plg, marks=0, voidevolution;


window, 1;
palette, "yarg.gp";
limits;
fma;plg, marks=0, domain(,100);

window, 4; 
palette, "yarg.gp";
limits;
fma;plg, marks=0, voltageevolution;
