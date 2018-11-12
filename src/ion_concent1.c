#include "syspara.h"

// The units of dnai is in mM.  Note that naiont should be multiplied by the
// cell capacitance to get the correct units.  Since cell capacitance = 1 uF/cm^2,
// it doesn't explicitly appear in the equation below.
// This holds true for the calculation of dki and dcai.

// Na concentration dose not use a function.
// K concentration dose not use a function.
void comp_iconcent (double x[])
{

	conc_nsr(x);
	conc_jsr(x);
	conc_itr(x);
	conc_cai(x);
	conc_cleft(x);

}

void conc_nsr (double x[])

// NSR Ca Ion Concentration Changes 
// iupbar = 0.00875;  // Max. current through iup channel (mM/ms)
// iup;               // Ca uptake from myo. to NSR (mM/ms)
// kleak;             // Rate constant of Ca leakage from NSR (ms^-1)
// ileak;             // Ca leakage from NSR to myo. (mM/ms)
// dnsr;              // Change in [Ca] in the NSR (mM)
// kmup = 0.00092;    // Half-saturation concentration of iup (mM)
// nsrbar = 15;       // Max. [Ca] in NSR (mM) 
// x[18] = [Ca]i (cai)
// x[19] = [Ca]_NSR (nsr)
// x[20] = [Ca]_JSR (jsr)

{

	double kleak;

	kleak = var.iupbar/var.nsrbar;
	var.ileak = kleak*x[19];

	var.iup = var.iupbar*x[18]/(x[18]+var.kmup);

}

void conc_jsr (double x[])

// JSR Ca Ion Concentration Changes 
// djsr;     // Change in [Ca] in the JSR (mM)
// tauon = 2;        // Time constant of activation of Ca release from JSR (ms)
// tauoff = 2;       // Time constant of deactivation of Ca release from JSR (ms)
// tcicr;            // t=0 at time of CICR (ms)
// irelcicr;         // Ca release from JSR to myo. due to CICR (mM/ms)
// csqnth = 8.75;    // Threshold for release of Ca from CSQN due to JSR ovreload (mM)
// gmaxrel = 150;    // Max. rate constant of Ca release from JSR due to overload (ms^-1)
// grelbarjsrol;     // Rate constant of Ca release from JSR due to overload (ms^-1)
// greljsrol;        // Rate constant of Ca release from JSR due to CICR (ms^-1)
// tjsrol;           // t=0 at time of JSR overload (ms)
// ireljsrol;        // Ca release from JSR to myo. due to JSR overload (mM/ms)
// swspontan = 0;    // switch of spontaneous release
// csqnbar = 10;     // Max. [Ca] buffered in CSQN (mM)
// kmcsqn = 0.8;     // Equalibrium constant of buffering for CSQN (mM)
// bjsr;             // b Variable for analytical computation of [Ca] in JSR (mM)
// cjsr;             // c Variable for analytical computation of [Ca] in JSR (mM)
// on;               // Time constant of activation of Ca release from JSR (ms)
// off;              // Time constant of deactivation of Ca release from JSR (ms)
// magrel;           // Magnitude of Ca release
// dICa;          // Rate of change of Ca entry
// dICa_new;       // New rate of change of Ca entry
// ICa_total_old;        // Old rate of change of Ca entry
// x[18] = [Ca]i (cai)
// x[19] = [Ca]_NSR (nsr)
// x[20] = [Ca]_JSR (jsr)
// t_cicr = time of CICR (tcicr)
// t_overload = time of JSR overload (tjsrol)

{

	double magrel;
	double on,off;
	double greljsrol;
	

	var.dICa_total_new = (var.ICa_total - var.ICa_total_old)/var.dt;

	if( x[0] > -35.0 && var.dICa_total_new < var.dICa_total && var.boolien == 0){
	//if( x[0] > -35.0 && var.dICa_total_new < var.dICa_total){
		var.boolien = 1;
		var.t_cicr = 0;
		printf("reset t_cicr.\n");
	}

	on = 1.0/(1.0+exp((-1.0*var.t_cicr+4.0)/0.5));
	off = (1.0-1.0/(1.0+exp((-1.0*var.t_cicr+4.0)/0.5)));
	
	magrel = 1.0/(1.0+exp(((var.ilca+var.icab+var.ipca-2.0*var.inaca+var.icat)+5.0)/0.9));

	var.Irel_cicr = var.gmaxrel*on*off*magrel*(x[20]-x[18]);
	
	var.t_cicr += var.dt;

	greljsrol = var.grelbarjsrol*(1.0-exp(-var.t_overload/var.tauon))*exp(-var.t_overload/var.tauoff);
	var.Irel_jsr_overload = greljsrol*(x[20]-x[18]);
	
	var.t_overload += var.dt;

	var.csqn = var.csqnbar*(x[20]/(x[20] + var.kmcsqn));

}
 void conc_itr (double x[])

// Translocation of Ca Ions from NSR to JSR
// itr;      // Translocation current of Ca ions from NSR to JSR (mM/ms)
// tautr = 180;      // Time constant of Ca transfer from NSR to JSR (ms)
// x[19] = [Ca]_NSR (nsr)
// x[20] = [Ca]_JSR (jsr)
{
	var.itr = (x[19]-x[20])/var.tautr; 
}

void conc_cai (double x[])

// Myoplasmic Ca Ion Concentration Changes 
// dcai;    // Change in myoplasmic Ca concentration (mM)
// catotal; // Total myoplasmic Ca concentration (mM)
// bmyo;    // b Variable for analytical computation of [Ca] in myoplasm (mM)
// cmyo;    // c Variable for analytical computation of [Ca] in myoplasm (mM)
// dmyo;    // d Variable for analytical computation of [Ca] in myoplasm (mM)
// gpig;    // Tribute to all the guinea pigs killed for the advancement of knowledge
// cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM)
// trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM)
// kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM)
// kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM)

{

	var.trpn = var.trpnbar*(x[18]/(x[18]+var.kmtrpn));
	var.cmdn = var.cmdnbar*(x[18]/(x[18]+var.kmcmdn));

}


void conc_cleft (double x[])

// Extracellular Ion Concentration Changes 
// dko;    // Change in Cleft K Concentration (mM)
// dnao;   // Change in Cleft Na Concentration (mM)
// dcao;   // Change in Cleft Ca Concentration (mM)
// taudiff = 1000; // Diffusion Constant for Ion Movement from Bulk Medium to Cleft Space

{
	double dummy;

	dummy = 0;

}
