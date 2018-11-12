#include "syspara.h"

// L-type calcium current
// Simplyfied ver of L-type Ca channel from Findray's 2008 model.
// From Findlay et al., 2008, Prog.Biophys. Mol and Biol.
//
//
// ilcatot; // Total current through the L-type Ca channel (uA/uF)
// ibarca;  // Max. Ca current through Ca channel (uA/uF)
// ibarna;  // Max. Na current through Ca channel (uA/uF)
// ibark;   // Max. K current through Ca channel (uA/uF)
// dss;     // Steady-state value of activation gate d 
// taud;    // Time constant of gate d (ms^-1)
// fss;     // Steady-state value of inactivation gate f
// tauf;    // Time constant of gate f (ms^-1)
// fca;     // Ca dependant inactivation gate
// kmca = 0.0006;     // Half-saturation concentration of Ca channel (mM)
// pca = 0.00054;     // Permiability of membrane to Ca (cm/s)
// gacai = 1;         // Activity coefficient of Ca
// gacao = 0.341;     // Activity coefficient of Ca
// pna = 0.000000675; // Permiability of membrane to Na (cm/s)
// ganai = 0.75;     // Activity coefficient of Na
// ganao = 0.75;     // Activity coefficient of Na
// pk = 0.000000193;  // Permiability of membrane to K (cm/s)
// gaki = 0.75;      // Activity coefficient of K
// gako = 0.75;      // Activity coefficient of K
// ratgca = 0.5;      // rate ?
// fcarat = 0.25;     // rate of the slow gate during CDI process
// aCDI = 0.00128;    // rate constant of inactivation
// bCDI = 0.003;      // rate constant of recovery 
// x[18] = [Ca]i (cai)
// x[23] = [Ca]o (cao)

void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	
	double f;   // VDI variable
	double fca; // CDI variable
	double expCa;
     
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

// activation
// see Suzuki et al., 2008, Prog.Biophys. Mol and Biol.
// Eqs.(2) and (3)
	var.dss = var.Tdss[iV]*d2 + var.Tdss[iV+1]*d1;
	var.taud = var.Ttaud[iV]*d2 + var.Ttaud[iV+1]*d1;

// VDI 
// see Suzuki et al., 2008, Prog.Biophys. Mol and Biol.
// Eqs.(5)-(8)
    var.fss1=var.Tfss1[iV]*d2 + var.Tfss1[iV+1]*d1;
    var.fss2=var.Tfss2[iV]*d2 + var.Tfss2[iV+1]*d1;
    var.tauf1=var.Ttauf1[iV]*d2 + var.Ttauf1[iV+1]*d1;
    var.tauf2=var.Ttauf2[iV]*d2 + var.Ttauf2[iV+1]*d1;;

	f = 0.100+0.524*x[5]+0.376*x[6];

	expCa = var.TexpCa[iV]*d2 + var.TexpCa[iV+1]*d1;
	if(fabs(x[0]) > 0.001) {   // NaN will be occured at v=0
		var.ibarca = var.ratgca*var.pca*var.zca*var.zca*
			((x[0]*F*F)/(R*T))*((var.gacai*x[18]*expCa-var.gacao*x[23])/(expCa-1.0));
		var.ibarna = 0;

		if(var.ibarca < 0.0) {
			var.ibark = var.ratgca*540*var.pca*x[0]/(1.0-var.ibarca/(1.227*var.ratgca));
		 }else{
			var.ibark=var.ratgca*540*var.pca*x[0];
		 }
	 }
	 else {
		 var.ibarca = var.ratgca*var.pca*var.zca*F*(
		 	var.gacai*x[18]*(1.0+(var.zca*x[0]*F)/(R*T))-var.gacao*x[23]);
		 var.ibarna = 0;
		 var.ibark = var.ratgca*540*var.pca*x[0]/(1.0-var.ibarca/(1.227*var.ratgca));
	 }
	
	var.ilcarat=x[8]*x[8]*var.fcarat/(x[8]*x[8]*var.fcarat+x[7]*x[7]*(1.0-var.fcarat));

	var.fcatau1=1.0/(-var.aCDI*var.ilca*(1.0-var.ilcarat) + var.bCDI);
	var.fcatau2=1.0/(-var.aCDI*var.ilca*var.ilcarat + var.bCDI);
	var.fcass1=var.bCDI*var.fcatau1;
	var.fcass2=var.bCDI*var.fcatau2;

	fca=x[7]*x[7]*(1.0-var.fcarat)+x[8]*x[8]*var.fcarat;

	var.ilca = x[4]*f*fca*var.ibarca;
	var.ilcana = x[4]*f*fca*var.ibarna;
	var.ilcak = x[4]*f*fca*var.ibark;

}

