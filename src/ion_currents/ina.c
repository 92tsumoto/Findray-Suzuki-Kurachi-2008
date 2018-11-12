#include "syspara.h"

// Fast Sodium Current (time dependant) */
// gna;    // Max. Conductance of the Na Channel (mS/uF)
// am;     // Na alpha-m rate constant (ms^-1)
// bm;     // Na beta-m rate constant (ms^-1)
// ah;     // Na alpha-h rate constant (ms^-1)
// bh;     // Na beta-h rate constant (ms^-1)
// aj;     // Na alpha-j rate constant (ms^-1)
// bj;     // Na beta-j rate constant (ms^-1)
// x[16]=[Na]i (nai)
// x[21]=[Na]o (nao)

void comp_ina(double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.Ena=((R*T)/F)*log(x[21]/x[16]); // [Na]i = x[16],[Na]o = x[21]

	var.am = var.Tam[iV]*d2 + var.Tam[iV+1]*d1;
	var.bm = var.Tbm[iV]*d2 + var.Tbm[iV+1]*d1;
	var.ah = var.Tah[iV]*d2 + var.Tah[iV+1]*d1;
	var.bh = var.Tbh[iV]*d2 + var.Tbh[iV+1]*d1;
	var.aj = var.Taj[iV]*d2 + var.Taj[iV+1]*d1;
	var.bj = var.Tbj[iV]*d2 + var.Tbj[iV+1]*d1;

	var.ina = var.gna*x[1]*x[1]*x[1]*x[2]*x[3]*(x[0]-var.Ena);
}

