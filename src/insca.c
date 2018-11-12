#include "syspara.h"

/* Nonspecific Ca-activated Current */
// ibarnsna;              // Max. Na current through NSCa channel (uA/uF)
// ibarnsk;               // Max. K current through NSCa channel (uA/uF)
// pnsca = 0.000000175;   // Permiability of channel to Na and K (cm/s)
// kmnsca = 0.0012;       // Half-saturation concentration of NSCa channel (mM)
// ganai = 0.75;          // Activity coefficient of Na
// ganao = 0.75;          // Activity coefficient of Na
// gaki = 0.75;           // Activity coefficient of K
// gako = 0.75;           // Activity coefficient of K
// x[16] = [Na]i (nai)
// x[17] = [K]i (ki)
// x[18] = [Ca]i (cai)
// x[21] = [Na]o (nao)
// x[22] = [K]o (ko)
// x[23] = [Ca]o (cao)

void comp_insca (double x[])
{

	var.ibarnsna = var.pnsca*var.zna*var.zna*((x[0]*F*F)/(R*T))*((var.ganai*x[16]*exp((var.zna*x[0]*F)/(R*T))-var.ganao*x[21])/(exp((var.zna*x[0]*F)/(R*T))-1.0));
	var.ibarnsk = var.pnsca*var.zk*var.zk*((x[0]*F*F)/(R*T))*((var.gaki*x[17]*exp((var.zk*x[0]*F)/(R*T))-var.gako*x[22])/(exp((var.zk*x[0]*F)/(R*T))-1.0));

	var.insna = var.ibarnsna/(1.0+pow(var.kmnsca/x[18],3));
	var.insk = var.ibarnsk/(1.0+pow(var.kmnsca/x[18],3));

}
