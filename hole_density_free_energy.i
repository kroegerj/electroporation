/* DOCUMENT modela.i

  You can adjust temperature, initial conditions,
and so on,  below.

*/
// initialize T = 0
// get gaussian random number generator (random_n)

//		Import libraries

#include "random.i"

//		 Declare constants and variables with corresponding units

random_seed,.1234567
k_BT=0.5   // Tc is around 0.5 apparently
MCS=1000
radius=100;
fixedradius=40.0;
DEL_X=0.75;
time=40000; 
height=array(0.0,[1,radius]);
n=array(0.0,[1,radius]);
freenergy=array(0.0,[1,radius]);
center=int(radius/2);
rad=array(0.5,[1,radius]);
densityofholes=array(0.0,[1,radius]);
K=array(0.0,[1,radius]);
holesatrimevolution=array(0.0,[1,time]);
voltageevolution=array(0.0,[1,time]);
equilibriumdensityevolution=array(0.0,[1,time]);
energyevolution=array(0.0,[1,time]);
holesaverageevolution=array(0.0,[1,time]);
equilibriumdensity=array(1.0,[1,radius]);


for(i=1;i<=radius;i++){
	rad(i)=abs(i-center)*DEL_X;}

n0=1.0;
N0=1.5e5;
kappas=10.0;//de facon a ce que f(0)\sim 10 dyne/cm en accord avec Kuzmin et al.//les unites de kappas sont en kT/nm^2.
epsilon1=0.02;
epsilon2=0.009;
deltat=0.001;
hsquared=5.0;//5.0 nm^2 d'apres Kuzmin et al.
lambdasquared=8.0;//en accord avec Kuzmin et al. (2005) si \lambda^2<h_s^2
splay=-0.1;// units are 1/nm from Kozlovsky (2002)
beta=0.502027//in terms of cm/dynes at room temperature.
conductanceperdensity=3.77e-10;//in pScm^2;
q=2.46; // d'apres Debruin
alpha=100.0;//units in cm^{-1}ms^{-1}
equilibriumelectroporationvoltage=70.0;
equilibriumcurrentvoltage=1.0; //en mV
calciumout=0.75; 
//C_o min=C_i*exp(equilibriumcurrentvoltage/25.8)=0.0575 for
// C_i=0.01 and equilibriumcurrentvoltage=45.0;
calciumin=0.01;
/* Initial voltage across membrane. */
voltage=25.8*log(calciumout/calciumin);
/* Boundary conditions for surface tension calculations. */
n(center)=n0;
n(1)=0.0;
n(radius)=0.0;
equilibriumdensity=N0*exp(q*(voltage/equilibriumelectroporationvoltage)^2.0+beta*freenergy);
/* Evolution of monomers axis director n. */
n=n-deltat*(hsquared^2.0*lap2(n)+4.0*(hsquared-lambdasquared)*lap(n)+4.0*n);
densityofholes=densityofholes+deltat*(K-K*densityofholes/equilibriumdensity);

//		Loop over time for a given set of parameters

for(t=1;t<=time;t++){

holesdensityaverage=0.0;

for(i=1;i<=radius;i++){holesdensityaverage=holesdensityaverage+densityofholes(i);}
holesdensityaverage=holesdensityaverage/radius;

/* Update of voltage across the membrane. */

voltage=25.8*log(calciumout/calciumin);

/* Various scenarios for calcium flow and conduction through the membrane. */

calciumin=calciumin+deltat*conductanceperdensity*holesdensityaverage*(voltage-equilibriumcurrentvoltage);

voltage=25.8*log(calciumout/calciumin);

K=alpha*exp((voltage/equilibriumelectroporationvoltage)^2.0+beta*freenergy);

equilibriumdensity=N0*exp(q*(voltage/equilibriumelectroporationvoltage)^2.0+beta*freenergy);

densityofholes=densityofholes+deltat*(K-K*densityofholes/equilibriumdensity);

/* Boundary conditions and evolution for monomers axis directors. */

n(center)=n0;
n(1)=0.0;
n(radius)=0.0;

n=n-deltat*(hsquared^2.0*lap2(n)+4.0*(hsquared-lambdasquared)*lap(n)+4.0*n);

densityofholes=densityofholes+deltat*(K-K*densityofholes/equilibriumdensity);


/* Record of time evolution of various quantities. */


holesaverageevolution(t-1)=holesdensityaverage;
holesatrimevolution(t-1)=densityofholes(center-1);
voltageevolution(t-1)=voltage;
equilibriumdensityevolution(t-1)=equilibriumdensity(center-1);

freenergy=kappas/2.0*(lambdasquared*((grad(n)+splay)^2.0)+(hsquared*lap(n)/2.0+n)^2.0-lambdasquared*splay^2.0);
}

holesatrimevolution(t-1)=holesatrimevolution(t-2);
voltageevolution(t-1)=voltageevolution(t-2);
equilibriumdensityevolution(t-1)=equilibriumdensityevolution(t-2);



for(i=1;i<=radius;i++){
height(i)=height(i-1)+hsquared/2.0*grad(n)(i);}



window, 1;
palette, "yarg.gp";
limits;
fma;plg, marks=0, freenergy;

window, 2;
palette, "yarg.gp";
limits;
fma;plg, marks=0, voltageevolution;

window, 4; 
palette, "yarg.gp";
limits;
fma;plg, marks=0, holesaverageevolution;

