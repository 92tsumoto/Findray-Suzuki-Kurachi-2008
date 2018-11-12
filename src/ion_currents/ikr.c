#include "syspara.h"

// Rapidly Activating Potassium Current 
// gkr;   // Channel Conductance of Rapidly Activating K Current (mS/uF)
// ekr;   // Reversal Potential of Rapidly Activating K Current (mV)
// xrss;  // Steady-state value of inactivation gate xr
// xrssblk2;  // block parameter
// tauxr; // Time constant of gate xr (ms^-1)
// tauxrblk2;  // block parameter
// r;     // K time-independant inactivation
// x[17] = [K]i (ki)
// x[22] = [K]o (ko)

void comp_ikr (double x[])
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	double gkr;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.Ekr = ((R*T)/F)*log(x[22]/x[17]); // x[17] = [K]i, x[22] = [K]o
	gkr = 0.02614*sqrt(x[22]/5.4); 

	var.xrss = var.Txrss[iV]*d2 + var.Txrss[iV+1]*d1;
	var.tauxr = var.Ttauxr[iV]*d2 + var.Ttauxr[iV+1]*d1;
	var.r = var.Tr[iV]*d2 + var.Tr[iV+1]*d1;
	
	var.ikr = var.ikrf*gkr*x[11]*var.r*(x[0]-var.Ekr);
	//if(drugchan==1) ikr *= blockd;

}
