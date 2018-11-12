#include "syspara.h"

// Current through T-type Ca Channel */
// gcat;    // Max. Conductance of the T-type Ca channel (mS/uF)
// eca;     // Reversal Potential of the T-type Ca channel (mV)
// bss;     // Steady-state value of activation gate b 
// taub;    // Time constant of gate b (ms^-1)
// gss;     // Steady-state value of inactivation gate g
// taug;    // Time constant of gate g (ms^-1)
// x[18] = [Ca]i (cai)
// x[23] = [Ca]o (cao)

void comp_icat(double x[])
{

MKL_INT iV=0;
double V1,V2,d1,d2;
     
V1 = (x[0]+200)*dvm;
V2 = (int)V1;
d1 = V1-V2;
d2 = 1.0-d1;
iV = (int)V2;

	var.bss = var.Tbss[iV]*d2 + var.Tbss[iV+1]*d1;
	var.taub = var.Ttaub[iV]*d2 + var.Ttaub[iV+1]*d1;
	var.gss = var.Tgss[iV]*d2 + var.Tgss[iV+1]*d1;
	var.taug = var.Ttaug[iV]*d2 + var.Ttaug[iV+1]*d1;
	
	var.Eca = (R*T/(var.zca*F))*log(x[23]/x[18]); // [Ca]i = x[18],[Ca]o = x[23]

	var.icat = var.gcat*x[9]*x[9]*x[10]*(x[0]-var.Eca);

}
