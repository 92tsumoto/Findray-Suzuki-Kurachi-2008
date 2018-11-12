#include "syspara.h"

// Potassium Current (time-independant)
// gki;    // Channel Conductance of Time Independant K Current (mS/uF)
// aki;    // K alpha-ki rate constant (ms^-1)
// bki;    // K beta-ki rate constant (ms^-1)
// kin;    // K inactivation
// x[22] = [K]o (ko)

void comp_iki (double x[])
{

	double gki,aki,bki,kin;


	var.Eki = var.Ekr;
	gki = 0.75*sqrt(x[22]/5.4);

	aki = 1.02/(1.0+exp(0.2385*(x[0]-var.Eki-59.215)));
	bki = (0.49124*exp(0.08032*(x[0]-var.Eki+5.476))+exp(0.06175*(x[0]-var.Eki-594.31)))/(1.0+exp(-0.5143*(x[0]-var.Eki+4.753)));

	kin = aki/(aki+bki);

	var.iki = gki*kin*(x[0]-var.Eki);

	//if(drugchan==5) iki *= blockd;

}

