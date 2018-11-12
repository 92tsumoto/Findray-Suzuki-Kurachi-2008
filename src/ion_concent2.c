#include "syspara.h"

void comp_iconcent2 (double x[])
//x[20] = var.jsr;
{

// Update Ca_JSR concentration x[20] = [Ca]_JSR
	var.bjsr = var.csqnbar - var.csqn - x[20] + var.kmcsqn;
	var.cjsr = var.kmcsqn*(var.csqn + x[20]);
	x[20] = (sqrt(var.bjsr*var.bjsr+4.0*var.cjsr)-var.bjsr)/2.0;

// Update intracellular Ca concentration x[18] = [Ca]_i

	var.Ca_total = var.trpn + var.cmdn + x[18];
	var.bmyo = var.cmdnbar + var.trpnbar + var.kmtrpn + var.kmcmdn - var.Ca_total;
	var.cmyo = (var.kmcmdn*var.kmtrpn)
	 			-(var.Ca_total*(var.kmtrpn + var.kmcmdn))
	 			+(var.trpnbar*var.kmcmdn)
	 			+(var.cmdnbar*var.kmtrpn);
	var.dmyo = -1.0*var.kmtrpn*var.kmcmdn*var.Ca_total;
	var.gpig = sqrt(var.bmyo*var.bmyo-3.0*var.cmyo);

	x[18] = (2.0*var.gpig/3.0)*cos(acos((9.0*var.bmyo*var.cmyo-2.0*var.bmyo*var.bmyo*var.bmyo-27.0*var.dmyo)/(2.0*pow((var.bmyo*var.bmyo-3.0*var.cmyo),1.5)))/3.0)-(var.bmyo/3.0);

}



