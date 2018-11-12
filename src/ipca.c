#include "syspara.h"

// Sarcolemmal Ca Pump 
// ibarpca = 1.15; // Max. Ca current through sarcolemmal Ca pump (uA/uF)
// kmpca = 0.0005; // Half-saturation concentration of sarcolemmal Ca pump (mM)
// x[18] = [Ca]i (cai)

void comp_ipca (double x[])
{

	var.ipca = (var.ibarpca*x[18])/(var.kmpca+x[18]);

}
