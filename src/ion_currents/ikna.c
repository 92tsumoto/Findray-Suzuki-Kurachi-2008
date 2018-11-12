#include "syspara.h"

// Na-Activated K Channel 
// pona;   // Open probability dependant on Nai
// pov;    // Open probability dependant on Voltage
// ekna;   // Reversal potential
// gkna = 0.12848;   // Maximum conductance (mS/uF)
// nkna = 2.8;       // Hill coefficient for Na dependance
// kdkna = 66;       // Dissociation constant for Na dependance(mM)
// x[16] = [Na]i (nai)

void comp_ikna (double x[])
{

	double pona,pov;

	var.Ekna = var.Ekr;

	pona = 0.85/(1.0 + pow((var.kdkna/x[16]),2.8));
	pov = 0.8-0.65/(1.0 + exp((x[0]+125.0)/15.0));
	
	var.ikna = var.gkna*pona*pov*(x[0] - var.Ekna);

}
