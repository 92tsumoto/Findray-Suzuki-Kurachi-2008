#include "syspara.h"

// Ito Transient Outward Current
// Dumaine et al. Circ Res 1999;85:803-809

// gitodv;    // Maximum conductance of Ito
// ekdv;      // Reversal Potential of Ito
// rvdv;      // Time independant voltage dependence of Ito
// azdv;      // Ito alpha-z rate constant
// bzdv;      // Ito beta-z rate constant
// tauzdv;    // Time constant of z gate
// zssdv;     // Steady-state value of z gate
// aydv;      // Ito alpha-y rate constant
// bydv;      // Ito beta-y rate constant
// tauydv;    // Time constant of y gate
// yssdv;     // Steady-state value of y gate

void comp_ito (double x[])
{

	
	var.Ekdv = var.Ekr;

	var.rvdv = exp(x[0]/100.0);

	var.azdv = (10.0*exp((x[0]-40.0)/25.0))/(1.0 + exp((x[0]-40.0)/25.0));
	var.bzdv = (10.0*exp(-(x[0]+90.0)/25.0))/(1.0 + exp(-(x[0]+90.0)/25.0));

	var.aydv = 0.015/(1.0+exp((x[0]+60.0)/5.0));
	var.bydv = (0.1*exp((x[0]+25.0)/5.0))/(1.0+exp((x[0]+25.0)/5.0));

	var.ito = var.itof*0.5*x[14]*x[14]*x[14]*x[15]*var.rvdv*(x[0]-var.Ekdv);

}
