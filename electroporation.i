
%Very simple code for langevin equation for membrane electroporation. Published in:
%Curvature-Driven Pore Growth in Charged Membranes during Charge-Pulse and Voltage-Clamp Experiments
%Biophys J. 2009 Feb 4; 96(3): 907–916.
%doi:  [10.1016/j.bpj.2008.10.035]
%Author: Jens Kroeger (2007)
%Free for reproduction with proper citation



%Declaration of main variables of the model.



#include "lmfit.i"



k_BT=0.5   // Tc is around 0.5 apparently
MCS=1000



LENGTHX=LENGTHY=100
CMIN=-1 //pli limits
CMAX=1
random_seed,.1234567
psi=array(1,[2,LENGTHX,LENGTHY])  // initialize T = 0
// get gaussian random number generator (random_n)
#include "random.i"

func grad(mat)
{	// near neighbours
	de=roll(mat,1)-roll(mat,-1);
	return(de/DEL_X/2.0)
}

func lap(mat)
{	// near neighbours
	de=roll(mat,1)+roll(mat,-1)-2.0*mat;
	return(de/DEL_X/DEL_X)
}

result=array(0.0,[1,12]);
//for(s=1;s<=12;s++){
s=1;


time=100000;

evol=array(0.0,[1,time]);
voltagetime=array(0.0,[1,time]);
voltagederivative=array(0.0,[1,time]);
criticalradius=array(0.0,[1,time]);
calctime=array(0.0,[1,time]);
timearray=array(0.0,[1,time]);
potentialofr=array(0.0,[1,1000]);

for(i=1;i<=time;i++){
timearray(i)=i;
}



surfacetension=1.0E-21;//J/nm^2 d'apres Chizmadzhev et al. biophys. j. 69:2489-2500.
conductanceperarea=1.9E-18;//S/nm^2 Debruin and Krassowska Biophys. J. 77:1213

//radius=1.0E-9;
radius=1.0//nm
C=9.67E-6//J^.25 nm Neu and Krassowska
deltat=100.1;
ap=6.9E-20//F/nm^2 if h=1.0 nm instead of 5.0 as given by Neu and Krassowska;

capacitance=3.0E-9;//of membrane patch in Farads from Wilhelm et al.
//linetension=1.8E-11 //J/m //from Neu and Krassowska PRE 59:3471
linetension=1.0E-20 //J/nm from Neu and Krassowska PRE 59:3471
restpotential=0.0; //in mV
voltage=0.200//V


for(i=1;i<=1000;i++){
potentialofr(i)=-(0.01*(i+50))^2.0*3.1416*(surfacetension+ap*voltage^2.0)+2.0*3.1416*linetension*(0.01*(i+50))+(C/(0.01*(i+50)))^4.0;
}

for(t=2;t<=time;t++){

voltage=voltage-deltat*((conductanceperarea*3.1416*radius^2.0/capacitance)*(voltage-restpotential));

voltagederivative(t)=((conductanceperarea*3.1416*radius^2.0/capacitance)*(voltage-restpotential));
if(radius>0.0){
radius=radius+deltat*0.8E14*(2.0*3.1416*(surfacetension+ap*voltage^2.0)-2.0*3.1416*linetension/radius+4*C^4.0/(radius^6.0));
indice=t;
}
else{radius=0.0;}

criticalradius(t)=3.1416*linetension/(3.1416*surfacetension+ap*voltage^2.0);
evol(t)=radius;

voltagetime(t)=voltage;

}

experimentalconductance=voltagederivative/voltage;


window, 0;
palette, "yarg.gp";
limits;
fma;plg, marks=0, voltagetime;

window, 1;
palette, "yarg.gp";
limits;
fma;plg, marks=0, evol;

window, 2;
palette, "yarg.gp";
limits;
fma;plg, marks=0, experimentalconductance;

window, 3;
palette, "yarg.gp";
limits;
fma;plg, marks=0, potentialofr;

/*
flw=create("voltagetime.txt");
write,flw,timearray,voltagetime;
close,flw;

flw=create("evol.txt");
write,flw,timearray,evol;
close,flw;

flw=create("experimentalconductance.txt");
write,flw,timearray,experimentalconductance;
close,flw;
*/