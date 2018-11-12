#include "syspara.h"

void function(double x[],double f[],double t)
{
   	double f0; 

	comp_ina(x);
	comp_ical(x);
	comp_icat(x);
	comp_ikr(x);
	comp_iks(x);
	comp_iki(x);
	comp_ikp(x);
	comp_ito(x);
	comp_inaca(x);
	comp_inak(x);
	comp_ipca(x);
	comp_icab(x);
	comp_inab(x);
	comp_insca(x);

	var.INa_total = var.ina + var.inab + var.ilcana + 3.0*var.inak + 3.0*var.inaca;
	var.IK_total = var.ikr + var.iks + var.iki + var.ikp + var.ilcak - 2.0*var.inak + var.ito + var.ikna + var.Istim;
	var.ICa_total = var.ilca + var.icab + var.ipca - 2.0*var.inaca + var.icat;

	comp_iconcent(x);

	f[0] = -(var.INa_total+var.IK_total+var.ICa_total);
	f[1] = var.am*(1.0 - x[1]) - var.bm*x[1]; // m 
	f[2] = var.ah*(1.0 - x[2]) - var.bh*x[2]; // h
	f[3] = var.aj*(1.0 - x[3]) - var.bj*x[3]; // j
	f[4] = (var.dss - x[4])/var.taud;
	f[5] = (var.fss1 - x[5])/var.tauf1;
	f[6] = (var.fss2 - x[6])/var.tauf2;
	f[7] = (var.fcass1 - x[7])/var.fcatau1;
	f[8] = (var.fcass2 - x[8])/var.fcatau2;
	f[9] = (var.bss - x[9])/var.taub;         // 
	f[10] = (var.gss - x[10])/var.taug;         //
	f[11] = (var.xrss - x[11])/var.tauxr;       // 
	f[12] = (var.xs1ss - x[12])/var.tauxs1; 
	f[13] = (var.xs2ss - x[13])/var.tauxs2; 
	f[14] = var.azdv*(1.0 - x[14]) - var.bzdv*x[14]; 
	f[15] = var.aydv*(1.0 - x[15]) - var.bydv*x[15]; 
	f[16] = -(var.INa_total*var.acap)/(var.vmyo*var.zna*F);
	f[17] = -(var.IK_total*var.acap)/(var.vmyo*var.zk*F);
	f[18] = -(((var.ICa_total*var.acap)/(var.vmyo*var.zca*F))+((var.iup-var.ileak)*var.vnsr/var.vmyo)-(var.Irel_cicr*var.vjsr/var.vmyo)-(var.Irel_jsr_overload*var.vjsr/var.vmyo)); 
	f[19] = var.iup - var.ileak - var.itr*var.vjsr/var.vnsr;
	f[20] = var.itr - var.Irel_cicr - var.Irel_jsr_overload;
	f[21] = (var.Na_bulk - x[21])/var.tau_diff + var.INa_total*var.acap/(var.vcleft*F); 
	f[22] = (var.K_bulk  - x[22])/var.tau_diff + var.IK_total* var.acap/(var.vcleft*F);
	f[23] = (var.Ca_bulk - x[23])/var.tau_diff + var.ICa_total*var.acap/(var.vcleft*var.zca*F);

}

// Fast Sodium Current (time dependant) */

void comp_ina(double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.Ena=var.RTonF*log(x[21]/x[16]); // [Na]i = x[16],[Na]o = x[21]

	var.am = var.Tam[iV]*d2 + var.Tam[iV+1]*d1;
	var.bm = var.Tbm[iV]*d2 + var.Tbm[iV+1]*d1;
	var.ah = var.Tah[iV]*d2 + var.Tah[iV+1]*d1;
	var.bh = var.Tbh[iV]*d2 + var.Tbh[iV+1]*d1;
	var.aj = var.Taj[iV]*d2 + var.Taj[iV+1]*d1;
	var.bj = var.Tbj[iV]*d2 + var.Tbj[iV+1]*d1;

	var.ina = var.gna*x[1]*x[1]*x[1]*x[2]*x[3]*(x[0]-var.Ena);
}


// L-type calcium current
// Simplyfied ver of L-type Ca channel from Findray's 2008 model.
// From Findlay et al., 2008, Prog.Biophys. Mol and Biol.
//
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

// Current through T-type Ca Channel */

void comp_icat(double x[])
{

MKL_INT iV=0;
double V1,V2,d1,d2;
     
V1 = (x[0]+200)*dvm;
V2 = (int)V1;
d1 = V1-V2;
d2 = 1.0-d1;
iV = (int)V2;

	var.bss = var.Tbss[iV]*d2 + var.Tbss[iV+1]*d1;
	var.taub = var.Ttaub[iV]*d2 + var.Ttaub[iV+1]*d1;
	var.gss = var.Tgss[iV]*d2 + var.Tgss[iV+1]*d1;
	var.taug = var.Ttaug[iV]*d2 + var.Ttaug[iV+1]*d1;
	
	var.Eca = (R*T/(var.zca*F))*log(x[23]/x[18]); // [Ca]i = x[18],[Ca]o = x[23]

	var.icat = var.gcat*x[9]*x[9]*x[10]*(x[0]-var.Eca);

}


// Rapidly Activating Potassium Current 

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

}

// Slowly Activating Potassium Current 

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

}

// Potassium Current (time-independant) Ik1
void comp_iki (double x[])
{
        
        MKL_INT iV=0,iK=0;   
        double V1,V2,d1,d2;
        double K1,K2,dd1,dd2;
        double p1,p2,p3,p4,q1,q2,q3,q4;
        double aki,bki,kin;

        V1 = (x[0]+200)*dvm;
        V2 = (int)V1;
        d1 = V1-V2;
        d2 = 1.0-d1;
        iV = (int)V2;
        p1 = var.T1Ki[iV]*d2 + var.T1Ki[iV+1]*d1;
        p2 = var.T2Ki[iV]*d2 + var.T2Ki[iV+1]*d1;
        p3 = var.T3Ki[iV]*d2 + var.T3Ki[iV+1]*d1;
        p4 = var.T4Ki[iV]*d2 + var.T4Ki[iV+1]*d1;

        //K1 = x[17]*dvm;
        //K2 = (int)K1;
        //dd1 = K1-K2;
        //dd2 = 1.0-dd1;
        //iK = (int)K2;
        //q1 = var.T5Ki[iK]*dd2 + var.T5Ki[iK+1]*dd1;
        //q2 = var.T6Ki[iK]*dd2 + var.T6Ki[iK+1]*dd1;
        //q3 = var.T7Ki[iK]*dd2 + var.T7Ki[iK+1]*dd1;
        //q4 = var.T8Ki[iK]*dd2 + var.T8Ki[iK+1]*dd1;
		q1=exp(-0.2385*var.RTonF*log(x[22]/x[17]));
		q2=exp(-0.08032*var.RTonF*log(x[22]/x[17]));
		q3=exp(-0.06175*var.RTonF*log(x[22]/x[17]));
		q4=exp(0.5143*var.RTonF*log(x[22]/x[17]));

        //aki = 1.02/(1.0+exp(0.2385*(x[0]-var.Eki-59.215)));
        //bki = (0.49124*exp(0.08032*(x[0]-var.Eki+5.476))
        //              +exp(0.06175*(x[0]-var.Eki-594.31)))/(1.0+exp(-0.5143*(x[0]-var.Eki+4.753)));
        
        var.aki = 1.02/(1.0+p1*q1*var.c1_ki);
        var.bki = (0.49124*p2*q2*var.c2_ki+p3*q3*var.c3_ki)/(1.0+p4*q4*var.c4_ki);

        kin = var.aki/(var.aki+var.bki);

		var.gki = 0.75*sqrt(x[22]/5.4);

        var.iki = var.gki*kin*(x[0]-var.Ekr);

}

// Plateau Potassium Current + Ultra Rapid Potassium current 

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

// Ito Transient Outward Current
// Dumaine et al. Circ Res 1999;85:803-809

void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.rvdv = var.Trvdv[iV]*d2 + var.Trvdv[iV+1]*d1;
    var.azdv = var.Tazdv[iV]*d2 + var.Tazdv[iV+1]*d1;
    var.bzdv = var.Tbzdv[iV]*d2 + var.Tbzdv[iV+1]*d1;
    var.aydv = var.Taydv[iV]*d2 + var.Taydv[iV+1]*d1;
    var.bydv = var.Tbydv[iV]*d2 + var.Tbydv[iV+1]*d1;

	var.ito = var.itof*0.5*x[14]*x[14]*x[14]*x[15]*var.rvdv*(x[0]-var.Ekr);

}

// Sodium-Calcium Exchanger V-S

void comp_inaca (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double exp2,exp3; 

	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp2=var.Texp2[iV]*d2 + var.Texp2[iV+1]*d1;
	exp3=var.Texp3[iV]*d2 + var.Texp3[iV+1]*d1;

    var.inaca = var.c1*exp2*((exp3*x[16]*x[16]*x[16]*x[23]-x[21]*x[21]*x[21]*x[18])/
				(1.0+var.c2*exp2*(exp3*x[16]*x[16]*x[16]*x[23]+x[21]*x[21]*x[21]*x[18])));
    
}

// Sodium-Potassium Pump

void comp_inak (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
    double exp0,exp1; 

	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp0=var.Texp0[iV]*d2 + var.Texp0[iV+1]*d1;
	exp1=var.Texp1[iV]*d2 + var.Texp1[iV+1]*d1;

	var.sigma = (exp(x[21]/67.3)-1.0)/7.0;

	var.fnak = 1.0/(exp0+0.0365*var.sigma*exp1);

	var.inak = var.ibarnak*var.fnak*(1.0/(1.0+pow(var.kmnai/x[16],2.0)))*(x[22]/(x[22]+var.kmko));

}

// Sarcolemmal Ca Pump 

void comp_ipca (double x[])
{

	var.ipca = (var.ibarpca*x[18])/(var.kmpca+x[18]);

}

// Ca Background Current 

void comp_icab (double x[])
{

	// var.gcab = 0.003016;
	var.Ecan = ((R*T)/(var.zca*F))*log(x[23]/x[18]); // = Eca;

	var.icab = var.gcab*(x[0]-var.Ecan);

}

// Na Background Current 

void comp_inab (double x[])
{

	// var.gnab = 0.004;
	var.Enan = var.Ena;

	var.inab = var.gnab*(x[0] - var.Enan);

}

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


