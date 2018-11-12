#include "syspara.h"

// Plateau Potassium Current + Ultra Rapid Potassium current 
// gkp;    // Channel Conductance of Plateau K Current (mS/uF)
// ekp;    // Reversal Potential of Plateau K Current (mV)
// kp;     // K plateau factor      

void comp_ikp (double x[])
{
	double kp;

	MKL_INT iV=0;
	double V1,V2,d1,d2;
     
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.Ekp = var.Ekr;

	kp = var.Tkp[iV]*d2 + var.Tkp[iV+1]*d1;
	
	var.ikp = var.ikxf*0.02*kp*(x[0]-var.Ekp);

}
