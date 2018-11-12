#include "syspara.h"

// Slowly Activating Potassium Current 
// gks;   // Channel Conductance of Slowly Activating K Current (mS/uF)
// eks;   // Reversal Potential of Slowly Activating K Current (mV)
// xs1ss;  // Steady-state value of inactivation gate xs1
// tauxs1; // Time constant of gate xs1 (ms^-1)
// xs2ss;  // Steady-state value of inactivation gate xs2
// tauxs2; // Time constant of gate xs2 (ms^-1)
// prnak = 0.01833;  // Na/K Permiability Ratio
// x[16] = [Na]i (nai)
// x[17] = [K]i (ki)
// x[18] = [Ca]i (cai)
// x[21] = [Na]o (nao)
// x[22] = [K]o (ko)
// x[23] = [Ca]o (cao)

void comp_iks (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.xs1ss = var.Txs1ss[iV]*d2 + var.Txs1ss[iV+1]*d1;
	var.xs2ss = var.xs1ss;
	var.tauxs1 = var.Ttauxs1[iV]*d2 + var.Ttauxs1[iV+1]*d1;
	var.tauxs2 = 4.0*var.tauxs1;

	var.Eks = ((R*T)/F)*log((x[22]+var.prnak*x[21])/(x[17]+var.prnak*x[16]));

	var.iks = var.iksf*0.433*(1.0+0.6/(1.0+pow((0.000038/x[18]),1.4)))*x[12]*x[13]*(x[0]-var.Eks);

	//if(drugchan==2) iks *= blockd;

}
