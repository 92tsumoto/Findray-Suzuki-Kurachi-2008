#include "syspara.h"

// Ca Background Current 
// gcab;  // Max. conductance of Ca background (mS/uF)
// Ecan;  // Nernst potential for Ca (mV)
// x[18] = [Ca]i (cai)
// x[23] = [Ca]o (cao)

void comp_icab (double x[])
{

	// var.gcab = 0.003016;
	var.Ecan = ((R*T)/(var.zca*F))*log(x[23]/x[18]); // = Eca;

	var.icab = var.gcab*(x[0]-var.Ecan);

}
