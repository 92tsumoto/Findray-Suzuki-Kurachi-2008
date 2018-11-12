/* produced by Tsumoto. K 2008.10.27 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3;
int mode = 1;
int P = 2;
int beats = 30;

typedef double Number;

main(argc,argv)
int argc;
char **argv;
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double time;
	double h;
	double v_old,dvdt,dvdt_new;
	double t_stok;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;

/* Action Potential Duration and Max. Info */
	//double vmax [beats] ; // Max. Voltage (mV)
	//double dvdtmax [beats] ; // Max. dv/dt (mV/ms)
	double *vmax ; // Max. Voltage (mV)
	double *dvdtmax ; // Max. dv/dt (mV/ms)
	//double apd [beats] ; // Action Potential Duration
	double *apd; // Action Potential Duration
	//double toneapd [beats] ; // Time of dv/dt Max.
	double *toneapd; // Time of dv/dt Max.
	//double ttwoapd [beats] ; // Time of 90% Repolarization
	double *ttwoapd; // Time of 90% Repolarization
	//double rmbp [beats] ; // Resting Membrane Potential
	double *rmbp; // Resting Membrane Potential
	//double nair [beats] ; // Intracellular Na At Rest
	//double cair [beats] ; // Intracellular Ca At Rest
	//double kir [beats] ; // Intracellular K At Rest
	double *nair; // Intracellular Na At Rest
	double *cair; // Intracellular Ca At Rest
	double *kir ; // Intracellular K At Rest
	double caimax [beats] ; // Peak Intracellular Ca

	vmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax=(Number *)calloc(beats,sizeof(Number));
	apd=(Number *)calloc(beats,sizeof(Number));
	toneapd=(Number *)calloc(beats,sizeof(Number));
	ttwoapd=(Number *)calloc(beats,sizeof(Number));
	rmbp=(Number *)calloc(beats,sizeof(Number));
	nair=(Number *)calloc(beats,sizeof(Number));
	cair=(Number *)calloc(beats,sizeof(Number));
	kir=(Number *)calloc(beats,sizeof(Number));
	var.Tam=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbm=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tah=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbh=(Number *)calloc(VNMAX,sizeof(Number));
	var.Taj=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbj=(Number *)calloc(VNMAX,sizeof(Number));
	var.Txs1ss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauxs1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Txrss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauxr=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tr=(Number *)calloc(VNMAX,sizeof(Number));
	var.T1Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.T2Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.T3Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.T4Ki=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tkp=(Number *)calloc(VNMAX,sizeof(Number));
	var.Trvdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tazdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbzdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Taydv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbydv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaub=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tgss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaug=(Number *)calloc(VNMAX,sizeof(Number));
	var.TexpCa=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tfss1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tfss2=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauf1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauf2=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp0=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp2=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp3=(Number *)calloc(VNMAX,sizeof(Number));
	if(vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL || ttwoapd==NULL 
		|| rmbp==NULL || nair==NULL || cair==NULL || kir==NULL
		|| var.Tam==NULL || var.Tah==NULL || var.Taj==NULL || var.Tbm==NULL || var.Tbh==NULL || var.Tbj==NULL
		|| var.Txs1ss==NULL || var.Ttauxs1==NULL || var.Txrss==NULL || var.Ttauxr==NULL || var.Tr==NULL
		|| var.Tkp==NULL || var.Trvdv==NULL || var.Tazdv==NULL || var.Tbzdv==NULL || var.Taydv==NULL || var.Tbydv==NULL
		|| var.Tbss==NULL || var.Ttaub==NULL || var.Tgss==NULL || var.Ttaug==NULL
		|| var.TexpCa==NULL || var.Tdss==NULL || var.Ttaud==NULL || var.Tfss1==NULL || var.Tfss2==NULL || var.Ttauf1==NULL || var.Ttauf2==NULL
		|| var.Texp0==NULL || var.Texp1==NULL || var.Texp2==NULL || var.Texp3==NULL ) exit(1);

	var.T5Ki=(Number *)calloc(KIMAX,sizeof(Number));
	var.T6Ki=(Number *)calloc(KIMAX,sizeof(Number));
	var.T7Ki=(Number *)calloc(KIMAX,sizeof(Number));
	var.T8Ki=(Number *)calloc(KIMAX,sizeof(Number));
	if(var.T5Ki==NULL || var.T6Ki==NULL || var.T7Ki==NULL || var.T8Ki==NULL ) exit(1);

//int i; // Stimulation Counter

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("para.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

// parameter inputs
	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}

	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}
		if (var.half){
			h = M_PI / var.m;
		}
		else {
			h = 1.0 / var.m;
		}
		h *= var.tsign[ii];
		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

// invariant constant
		var.RTonF = R*T/F;
// Cell Geometry */
		var.length = 0.01;       // Length of the cell (cm)
		var.a = 0.0011;     // Radius of the cell (cm)
		var.vcell = 1000*M_PI*var.a*var.a*var.length; // Cell Volume:3.801e-5 (uL)
		var.ageo = 2*M_PI*var.a*var.a+2.0*M_PI*var.a*var.length;  // eometric membrane area: 7.671e-5 (cm^2)
		var.acap = var.ageo*2;          // Capacitive membrane area: 1.534e-4 cm^2 (cm^2)
		var.vmyo = var.vcell*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
		var.vmito = var.vcell*0.26;     // Mitochondria volume (uL) = 26% for cell volume
		var.vsr = var.vcell*0.06;       // SR volume (uL)
		var.vnsr = var.vcell*0.0552;    // NSR volume (uL)
		var.vjsr = var.vcell*0.0048;    // JSR volume (uL)
		var.vcleft = var.vcell*0.12/0.88;  // Cleft volume (uL)

// Ion Valences
		var.zna = 1;  // Na valence
		var.zk = 1;   // K valence
		var.zca = 2;  // Ca valence

	// Max conductanve (constant)
		var.gna = 16;
		var.gcat = 0.05;
		var.gkna = 0.12848;
		var.kdkna = 66.0; 

	// classification by cell type
		//var.itof = 1.0; var.ikxf = 1.0; var.iksf = 1.0; var.ikrf = 1.0;    // Base
		//var.itof = 0.2263; var.ikxf = 0.5040; var.iksf = 0.1903; var.ikrf = 1.2800;    // Epi
		//var.itof = 0.1600; var.ikxf = 1.2699; var.iksf = 0.0951; var.ikrf = 0.5382;    // M
		//var.itof = 0.1131; var.ikxf = 1.2699; var.iksf = 0.1903; var.ikrf = 0.6400;    // End
		var.gkr = 0.02614;                      // real value: gkr*ikrf
		var.gkr_max = var.gkr*var.ikrf;
		var.gks = 0.433;                        // real value: gks_max*iksf
		var.gks_max = var.gks*var.iksf; 
		var.gkp = 0.02;                         // real value: gkp_max*ikxf
		var.gkp_max = var.gkp*var.ikxf; 
		var.gkur = 1.1*var.ikxf;
		var.gitodv = 0.5;                       // real value: gitodv*itof
		var.gitodv_max = var.gitodv*var.itof; 

	// L-type calcium current
		var.kmca = 0.0006;     // Half-saturation concentration of Ca channel (mM)
		var.pca = 0.00054;     // Permiability of membrane to Ca (cm/s)
		var.gacai = 1;         // Activity coefficient of Ca
		var.gacao = 0.341;     // Activity coefficient of Ca
		var.pna = 0.000000675; // Permiability of membrane to Na (cm/s)
		var.pk = 0.000000193;  // Permiability of membrane to K (cm/s)
		var.ratgca = 0.5;      // rate ?
		var.fcarat = 0.25;	   // rate of the slow gate during CDI process
		var.aCDI = 0.00128;	   // rate constant of inactivation
		var.bCDI = 0.003;	   // rate constant of recovery	

	//  Slowly Activating Potassium Current: IKs
		var.prnak = 0.01833;     // Na/K Permiability Ratio
	
	// Inward rectifier K current: Ik1
		//var.gki = 0.75*sqrt(var.ko/5.4);
		var.c1_ki = exp(-0.2385*59.215);
		var.c2_ki = exp(0.08032*5.476);
		var.c3_ki = exp(-0.06175*594.31);
		var.c4_ki = exp(-0.5143*4.753);

	//	ATP-Sensitive K Channel
	//	Note: If you wish to use this current in your simulations, there are additional 
	//  changes which must be made to the code as detailed in Cardiovasc Res 1997;35:256-272 

		var.natp = 0.24;           // K dependence of ATP-sensitive K current
	    var.nicholsarea = 0.00005; // Nichol's ares (cm^2)
		var.atpi = 3;              // Intracellular ATP concentraion (mM)
		var.hatp = 2;              // Hill coefficient
		var.katp = 0.250;          // Half-maximal saturation point of ATP-sensitive K current (mM)

	// Sodium-Calcium Exchanger V-S
		var.c1 = 0.00025;   // Scaling factor for inaca (uA/uF)
		var.c2 = 0.0001;    // Half-saturation concentration of NaCa exhanger (mM)
		var.gammas = 0.15;  // Position of energy barrier controlling voltage dependance of inaca

	// Sodium-Potassium Pump
		var.ibarnak = 2.25;   // Max. current through Na-K pump (uA/uF)
		var.kmnai = 10;    // Half-saturation concentration of NaK pump (mM)
		var.kmko = 1.5;    // Half-saturation concentration of NaK pump (mM)

	// Nonspecific Ca-activated Current
		var.pnsca = 0.000000175;   // Permiability of channel to Na and K (cm/s)
		var.kmnsca = 0.0012;       // Half-saturation concentration of NSCa channel (mM)
		var.ganai = 0.75;      // Activity coefficient of Na
		var.ganao = 0.75;      // Activity coefficient of Na
		var.gaki = 0.75;       // Activity coefficient of K
		var.gako = 0.75;       // Activity coefficient of K

	
	// Sarcolemmal Ca Pump
		var.ibarpca = 1.15;  // Max. Ca current through sarcolemmal Ca pump (uA/uF)
		var.kmpca = 0.0005;  // Half-saturation concentration of sarcolemmal Ca pump (mM)

	// Ca Background Current 
		var.gcab = 0.003016; // Max. conductance of Ca background (mS/uF)

	// Na Background Current 
		var.gnab = 0.004;    // Max. conductance of Na background (mS/uF)

	// NSR Ca Ion Concentration Changes 
		var.kmup = 0.00092;   // Half-saturation concentration of iup (mM)
		var.iupbar = 0.00875; // Max. current through iup channel (mM/ms)
		var.nsrbar = 15;      // Max. [Ca] in NSR (mM)

	// JSR Ca Ion Concentration Changes 
		var.tauon = 2.0;        // Time constant of activation of Ca release from JSR (ms)
		var.tauoff = 2.0;       // Time constant of deactivation of Ca release from JSR (ms)
		var.csqnth = 8.75;    // Threshold for release of Ca from CSQN due to JSR ovreload (mM)
		var.gmaxrel = 150;    // Max. rate constant of Ca release from JSR due to overload (ms^-1)
		var.grelbarjsrol = 0; // Rate constant of Ca release from JSR due to overload (ms^-1)
		var.swspontan = 0;    // switch of spontaneous release
		var.csqnbar = 10;     // Max. [Ca] buffered in CSQN (mM)
		var.kmcsqn = 0.8;     // Equalibrium constant of buffering for CSQN (mM)
		var.csqn = 6.97978;

	// Translocation of Ca Ions from NSR to JSR
		var.tautr = 180;      // Time constant of Ca transfer from NSR to JSR (ms)

	// Myoplasmic Ca Ion Concentration Changes 
		var.cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM)
		var.trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM)
		var.kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM)
		var.kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM)
		var.trpn = 0.0143923;
		var.cmdn = 0.00257849;


	// Extracellular Ion Concentration Changes 
		var.tau_diff = 1000; // Diffusion Constant for Ion Movement from Bulk Medium to Cleft Space

	// Initial Bulk MEdium ion concentrations
		var.Na_bulk = 132.0;     // Initial Bulk Medium Na (mM)
		//var.Na_bulk = 140.0;     // Initial Bulk Medium Na (mM)
		var.K_bulk = 4.5;      // Initial Bulk Medium K (mM)
		var.Ca_bulk = 1.8;     // Initial Bulk Medium Ca (mM)

	// Another parameter initial setting
		var.t_cicr = 1000;
		var.t_overload = 1000;
		var.boolien = 1;
		var.dICa_total = 0;

		printf("Istim=%lf\n",var.Istim_base);
		printf("csqn=%lf\n",var.csqn);

// Tablize exp functions.	
	//var.VNMAX=VNMAX;
	//var.dvm=5.0;
	make_ExpTable();

	// Initialization time
		time -= h;
		var.dt = h;
		var.beat = 0;

		ii = 0;
		
		while (1){
			eventloop(fp1,&mode,&P,x);

			for (j = 0; j< (var.m * var.l ); j++){
				t = h*j;
				time += h;

				if ( time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < h ){
					var.boolien = 0;
					apd[var.beat] =0;
					toneapd[var.beat] =0;
					ttwoapd[var.beat] =0;
					rmbp[var.beat] =x[0];
					nair[var.beat] = x[16];
					kir[var.beat] = x[17];
					cair[var.beat] = x[18];
					caimax[var.beat] = x[18];
					vmax[var.beat] = 0;

					printf("%d %lf %lf %lf %lf\n",var.beat,apd[var.beat],toneapd[var.beat],ttwoapd[var.beat],rmbp[var.beat]);
					printf("%lf %lf %f %f %e\n", x[0],x[1],x[2],x[3],x[4]);
					printf("%lf %lf %f %f %e\n", x[5],x[6],x[7],x[8],x[9]);
					printf("%lf %e %lf %lf %lf\n", x[10],x[11],x[12],x[13],x[14]);
					printf("%lf %lf %lf\n", x[15],x[16],x[17]);
					printf("%e %lf %lf\n", x[18],x[19],x[20]);
					printf("%lf %lf %f\n", x[21],x[22],x[23]);
					printf("time=%lf,Istim=%lf\n",time,var.Istim);
					printf("dvdtmax[%d]=%lf\n",var.beat,dvdtmax[var.beat]);
				}
				//if(var.beat < 30){
				if (time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < 0.5){
					var.Istim = var.Istim_base;
				} else {
					var.Istim = 0;
				}
				//}	

				if (fabs(time) > tend &&  tend != 0.0) break;
				var.ca_pre=x[18];
				v_old = x[0];

				eular(NN,h,x,t);

				comp_iconcent2 (x);

				dvdt_new = (x[0]-v_old)/h; // LRd -> dvdtnew
				//printf("dvdt_new=%lf\n",dvdt_new);

				if(var.beat>=0){
					if (x[0] > vmax[var.beat] )
						vmax[var.beat] = x[0];
					if (x[18] > caimax[var.beat] )
						caimax[var.beat] = x[18];
					if (dvdt_new > dvdtmax[var.beat] ){
						dvdtmax[var.beat] = dvdt_new;
						toneapd[var.beat] = time;
					}
					if (dvdt_new < 0 && x[0] >= (vmax[var.beat] -0.9*(vmax[var.beat]-rmbp[var.beat]) ) )
						ttwoapd[var.beat] = time;
				}

				if(var.csqn >= var.csqnth && var.t_overload > 50.0){
					var.grelbarjsrol = 4;
					t_stok = var.t_overload;
					var.t_overload = 0;
					printf("reset t_overload at %lf\n",t_stok);
				}

				if (var.pflag) orbit(&mode,x,dvdt_new);

				if (time>= (beats-1)*var.BCL && time < beats*var.BCL){
				//	fprintf(fp2,"%lf %lf %e %e %lf %lf %lf %lf %e %e %lf %lf\n",time,x[0],var.ina,var.ito,var.ICa_total,var.t_cicr,var.Irel_cicr,x[18],var.ilca,var.Ena,x[19],x[20]);
					fprintf(fp2,"%lf %lf\n",time,x[0]);
				//}
				//if (time>= 0.0*var.BCL && time < 3.0*var.BCL){
				//	fprintf(fp2,"%lf %lf %e %e %e %e %e %lf %lf %lf %lf %lf %lf\n",time,var.ilca,var.icab,var.ipca,var.inaca,var.icat,var.ibarca,var.ibark,x[4],x[5],x[6],x[7],x[8]);
				}
				
				dvdt = dvdt_new;
				var.ICa_total_old = var.ICa_total;
				var.dICa_total = var.dICa_total_new;

			}
			ii++;

			if(ii==var.BCL){
				var.beat++;
				ii = 0;
				if(var.beat > beats){
					for(w=0;w<=beats;w++){
					}
					for(w=0;w<=beats;w++){
						apd[w] = ttwoapd [w] -toneapd [w] ;
						fprintf(fp3,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%e\t%e\t%g\n",w,
							vmax[w],dvdtmax[w],apd[w],toneapd[w],ttwoapd[w],nair[w],kir[w],cair[w],caimax[w],rmbp[w]);
						printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%e\t%e\t%lf\n",w,
							vmax[w],dvdtmax[w],apd[w],toneapd[w],ttwoapd[w],nair[w],kir[w],cair[w],caimax[w],rmbp[w]);
					}
					exit(0);
				}
			}

			draw_p(&mode,P,x,dvdt);
			mouse(&mode,x,dvdt);
			if (fabs(time) > tend &&  tend != 0.0) break;

		}
		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
		free(vmax);free(dvdtmax);free(apd);free(toneapd);free(ttwoapd);
		free(rmbp);free(nair);free(cair);free(kir);
		free(var.Tam);free(var.Tah);free(var.Taj);free(var.Tbm);free(var.Tbh);free(var.Tbj);
		free(var.Txrss);free(var.Ttauxr);free(var.Tr);free(var.Txs1ss);free(var.Ttauxs1);
		free(var.Tkp);free(var.Tbss);free(var.Ttaub);free(var.Tgss);free(var.Ttaug);
		free(var.TexpCa);free(var.Tdss);free(var.Ttaud);free(var.Tfss1);free(var.Ttauf1);free(var.Tfss2);free(var.Ttauf2);
		free(var.Texp0);free(var.Texp1);free(var.Texp2);free(var.Texp3);
	}
}

