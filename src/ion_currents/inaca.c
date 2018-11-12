#include "syspara.h"

// Sodium-Calcium Exchanger V-S
// c1 = .00025;   // Scaling factor for inaca (uA/uF)
// c2 = 0.0001;   // Half-saturation concentration of NaCa exhanger (mM)
// gammas = 0.15;  // Position of energy barrier controlling voltage dependance of inaca
// x[16] = [Na]i
// x[17] = [K]i
// x[18] = [Ca]i
// x[21] = [Na]o
// x[22] = [K]o
// x[23] = [Ca]o

void comp_inaca (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double exp2,exp3; 

	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp2=var.Texp2[iV]*d2 + var.Texp2[iV+1]*d1;
	exp3=var.Texp3[iV]*d2 + var.Texp3[iV+1]*d1;

    var.inaca = var.c1*exp2*((exp3*x[16]*x[16]*x[16]*x[23]-x[21]*x[21]*x[21]*x[18])/(1.0+var.c2*exp2*(exp3*x[16]*x[16]*x[16]*x[23]+x[21]*x[21]*x[21]*x[18])));
    
}
