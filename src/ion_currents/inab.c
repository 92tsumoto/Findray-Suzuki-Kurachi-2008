#include "syspara.h"

// Na Background Current 
// gnab;  // Max. conductance of Ca background (mS/uF)
// Enan;  // Nernst potential for Ca (mV)

void comp_inab (double x[])
{

	// var.gnab = 0.004;
	var.Enan = var.Ena;

	var.inab = var.gnab*(x[0] - var.Enan);

}
