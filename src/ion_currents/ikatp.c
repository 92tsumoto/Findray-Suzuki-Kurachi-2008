#include "syspara.h"

// ATP-Sensitive K Channel
// ekatp;    // K reversal potential (mV)
// gkbaratp; // Conductance of the ATP-sensitive K channel (mS/uF)
// gkatp;    // Maximum conductance of the ATP-sensitive K channel (mS/uF)
// patp;     // Percentage availibility of open channels
// natp = 0.24;          // K dependence of ATP-sensitive K current
// nicholsarea = 0.00005; // Nichol's ares (cm^2)
// atpi = 3;            // Intracellular ATP concentraion (mM)
// hatp = 2;     // Hill coefficient
// katp = 0.250; // Half-maximal saturation point of ATP-sensitive K current (mM)

// Note: If you wish to use this current in your simulations,
// there are additional changes which must be made to the code 
// as detailed in Cardiovasc Res 1997;35:256-272 
// [K]o = x[22] extracellular K ion concentration

void comp_ikatp (double x[])
{

	var.Ekatp = var.Ekr;

	var.gkatp = 0.000195/var.nicholsarea;
	var.patp = 1.0/(1.0+(pow((var.atpi/var.katp),var.hatp)));
	var.gkbaratp = var.gkatp*var.patp*(pow((x[22]/4.0),var.natp));

	var.ikatp = var.gkbaratp*(x[0]-var.Ekatp);

}
