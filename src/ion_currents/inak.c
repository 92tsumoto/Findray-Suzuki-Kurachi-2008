#include "syspara.h"

// Sodium-Potassium Pump
// fnak;    // Voltage-dependance parameter of inak
// sigma;   // [Na]o dependance factor of fnak
// ibarnak = 2.25;   // Max. current through Na-K pump (uA/uF)
// kmnai = 10;    // Half-saturation concentration of NaK pump (mM)
// kmko = 1.5;    // Half-saturation concentration of NaK pump (mM)
// x[16] = [Na]i
// x[21] = [Na]o
// x[22] = [K]o

void comp_inak (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
    double exp0,exp1; 

	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp0=var.Texp0[iV]*d2 + var.Texp0[iV+1]*d1;
	exp1=var.Texp1[iV]*d2 + var.Texp1[iV+1]*d1;

	var.sigma = (exp(x[21]/67.3)-1.0)/7.0;

	var.fnak = 1.0/(exp0+0.0365*var.sigma*exp1);

	var.inak = var.ibarnak*var.fnak*(1.0/(1.0+pow(var.kmnai/x[16],2.0)))*(x[22]/(x[22]+var.kmko));

}
