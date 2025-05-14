************************* Run RBA model ********************
*	Original Author: Hoang Dinh (for scRBA model)
*	Latest Author: Wheaton Schroeder (for ctherm RBA model)
*	Written for performing bisection method to study C. therm growth on avicel where cellulsome synthesis is required
*
*	Cellulosome synthesis amounts from the Amador-Noguez dataset and from mapping of the cellulosome from Yoav et al 2017
*	and BLASTp from the ATCC 27405 strain enumerated by Blumer-Schuette et al. 2013.
*
* This code uses the bisection method to identify the maximum growth rate. At the highest feasible growth rate found it
* then maximizes the flux through an exchange reaction representing excreted fermentation product. 
* This makes this code useful for evaluating product yeilds at extreme growth
*
*	Latest Version: 05/31/2024
************************************************************

$setGlobal settings "nc_av_r_e_adhE"
$setGlobal identifier "f_nc_av_et_r_e_adhE"

$setGlobal usemodel "RBAmaxETOH"
*$setGlobal usemodel "RBAmaxAC"
*$setGlobal usemodel "RBAmaxETOHac"
*$setGlobal usemodel "RBAmaxH2"
*$setGlobal usemodel "RBAminH2"
*$setGlobal usemodel "RBAmaxGLYC"
*$setGlobal usemodel "RBAminGLYC"
*$setGlobal usemodel "RBAmixGLYC"

$INLINECOM /*  */
$include "./runRBA_settings_%settings%.txt"

* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1E3
$setGlobal nscaleback 1E-3

*set a tolerance for the bisection algorithm
$setGlobal tolerance 1E-7

*set arbitrary large number to use for BOUNDS
$setGlobal bigM 1E3

*maximum number of outer iterations to allow
$setGlobal out_iter_max 30

*maximum number of inner iterations to allow
$setGlobal inner_iter_max 10

options

	LP = cplex					/*Solver selection*/
	limrow = 0					/*number of equations listed, 0 is suppresed*/
	limcol = 0					/*number of variables listed, 0 is suppresed*/
	iterlim = 1000000			/*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8				/*decimal places for display statement*/
	reslim = 1000000			/*wall-clock time limit for solver in seconds*/
	sysout = on					/*solver status file report option*/
	solprint = on				/*solution printing option*/

;

Sets

	i													/*list of all metabolic species*/
$include "%species_path%"
	
	j													/*list of all reactions*/
$include "%rxns_path%"

	carbs(j)											/*list of reactions where the carbon yield is constrained by MFA values*/
$include "%car_frac_j_path%"

	prosyn(j)											/*list of protein synthesis reactions*/
$include "%prosyn_path%"

	prosyn_cell(prosyn)							/*list of protein synthesis reactions for proteins making up the cellulosome*/
$include "%prosyn_cell_rxns%"

	prosyn_met(prosyn)								/*list of metabolic protein synthesis reactions */
$include "%prosyn_met_rxns%"

	prosyn_ribo(prosyn)								/* list of ribosomal protein synthesis reactions */
$include "%prosyn_ribo_rxns%"						

	prowaste(j)											/*list of protein waste reactions*/
$include "%prowaste_path%"

  cellulosome_waste(prowaste)							/*list of protein waste reactions have to allow for some way to balance cellulosome proteins*/
$include "%cellulosome_waste_path%"

	translation(j)										/*list of translation (protein synthesis) reactions same as prosyn may have been needed due to multiple ribosomes in scRBA model may be removable later*/
$include "%trans_path%"

	uptake(j)											/*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"

	media(j)											/*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"

		t													/*bisection driver*/
	/ 1*100 /

;

Parameters

	S(i,j)												/*stoichiometric matrix*/
$include "%sij_path%"
	
	NAA(prosyn)											/*number of amino acids in protein synthesis reactions*/
$include "%prolen_path%"

	c_min_f(carbs)										/*lower bound limits for flux rates taken from MFA tested already in the base model*/
$include "%carb_min_path%"

	c_max_f(carbs)										/*upper bound limits for flux rates taken from MFA tested already in the base model*/
$include "%carb_max_path%"
	
	kapp(j)												/*apparent kinetic parameter*/
$include "%kapp_path%"

        prosyn_met_rates(prosyn_met)										/*rate of protein synthesis as measured by Amador-Noguez group*/
$include "%prosyn_met_exp%"

	prosyn_cell_rates(prosyn_cell)							/*concentration estimates of cellulosome components under avicel growth*/
$include "%prosyn_cell_exp%"

	muLB												/*current lower bound for bisection-based growth rate search*/

	muUB												/*current upper bound for bisection-based growth rate search*/

	mu													/*current lower bound for bisection-based growth rate search*/

	iter												/*iteration tracker for bisection method track outer loop*/

	iter_in												/*iteration counter tracks inner small perturbation loop*/

	iter_temp											/*updater for iteration tracker*/

	iter_soln											/*binary stores if iteration has a solution*/

	mu_temp												/*used to perturb mu when its unknown error*/

	delta_mu											/*used to report tolerance*/

	feas_mu												/*last feasible value of mu used to enforce rate when maximizing a product*/

	c_etoh(j)											/*weight for ethanol objective functions*/

	c_ac(j)												/*weight for acetate objective functions*/

	c_etoh_ac(j)										/*weight for ethanol + acetate objective functions*/

	c_h2(j)												/*weight for hydrogen objective functions*/

	calc_temp											/*temporarily holds the result fo a calculation*/

;

*default should be able to find feasible solution at zero
feas_mu = 0;

*set FBA objective
c_etoh(j) = 0;
c_etoh('RXN-EXCH_etoh_e_FWD-SPONT') = 1;

c_ac(j) = 0;
c_ac('RXN-EXCH_ac_e_FWD-SPONT') = 1;

c_etoh_ac(j) = 0;
c_etoh_ac('RXN-EXCH_etoh_e_FWD-SPONT') = 1;
c_etoh_ac('RXN-EXCH_ac_e_FWD-SPONT') = 1;

c_h2(j) = 0;
c_h2('RXN-EXCH_h2_e_FWD-SPONT') = 1;

*set the reaction to maximize once growth rate is determined

Variables

	z													/*objective variable*/

	v(j)												/*reaction rate variable*/

;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; 
v.up(j) = %bigM% * %nscale%;

*disable all protien waste
v.up(prowaste) = 0;
v.lo(prowaste) = 0;

*allow protein waste for cellulosome proteins since not all are directly linked to reactions or enzymes therefore then cannot be balanced by an enzyme load reaction
v.up(cellulosome_waste) = %bigM% * %nscale%;

* allow growth since won't be set in phenotype file
v.lo('BIOSYN-BIODIL') = 0;
v.up('BIOSYN-BIODIL') = %bigM% * %nscale%;

* Media
* disables all uptakes that are not in media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = %bigM% * %nscale%;

* in phenotype is set:
* substrate
* NGAM
* Biomass composition
* Changed this bit so I can specify phenotype in settings

$include "%phenotype%"

* Set your NGAM in phenotype.txt since NGAM value depends on growth-condition
* Set your biomass composition in phenotype.txt since biomass composition depends on growth-condition

*** EQUATION DEFINITIONS ***

Equations

	etoh_max_obj										/*objective equation for maximizing ethanol yield*/

	ac_max_obj											/*objective equation for maximizing acetate yield*/

	etoh_ac_max_obj										/*objective equation for maximizing acetate + ethanol yield*/

	h2_max_obj											/*objective equation for maximizing hydrogen gas yield*/

	h2_min_obj											/*objective equation for minimizing hydrogen gas yield*/

	glyc_max_obj										/*objective for manipulating glycolysis protein abundances maximization*/

	glyc_min_obj										/*objective for manipulating glycolysis protein abundances maximization*/

	glyc_mix_obj										/*objective for manipulating glycolysis protein abundances minimize upper maximize lower*/
	
	mb_const											/*mass balance equation*/
	
	RiboCapacity										/*ribosome capacity constraint based on number of amino acids*/
	
	NonModelProtAllo									/*limits proteome allocation to metabolic proteins*/
	
	enforce_cell_syn									/* requires cellulosome synthesis*/

	enforce_mu											/*require growth rate to match biomass reaction rate*/
	ribo_prot
$include %enz_cap_declares_path%						/*enzyme capacity constraint declarations applies kapp to the system*/

;

* define objective minimize proteome cost for RBA runs determining maximum growth rate
etoh_max_obj..							z =e= v('RXN-EXCH_etoh_e_FWD-SPONT');
ac_max_obj..							z =e= sum(j, c_ac(j) * v(j));
etoh_ac_max_obj..						z =e= sum(j, c_etoh_ac(j) * v(j));
h2_max_obj..							z =e= sum(j, c_h2(j) * v(j));
h2_min_obj..							z =e= sum(j, -c_h2(j) * v(j));

*glycolytic protein objectives
glyc_max_obj..							z =e= v('PROSYN-Clo1313_1954') + v('PROSYN-Clo1313_0993') + v('PROSYN-Clo1313_0489') + v('PROSYN-Clo1313_1831') + v('PROSYN-Clo1313_2015') + v('PROSYN-Clo1313_0997') + v('PROSYN-Clo1313_1876') + v('PROSYN-Clo1313_0238') + v('PROSYN-Clo1313_1875') + v('PROSYN-Clo1313_2093') + v('PROSYN-Clo1313_2095') + v('PROSYN-Clo1313_2094') + v('PROSYN-Clo1313_0966') + v('PROSYN-Clo1313_2745') + v('PROSYN-Clo1313_2092') + v('PROSYN-Clo1313_2090') + v('PROSYN-Clo1313_0949') + v('PROSYN-Clo1313_0415') + v('PROSYN-Clo1313_1879') + v('PROSYN-Clo1313_0020') + v('PROSYN-Clo1313_0021') + v('PROSYN-Clo1313_0022') + v('PROSYN-Clo1313_0023') + v('PROSYN-Clo1313_0382') + v('PROSYN-Clo1313_0383') + v('PROSYN-Clo1313_0384') + v('PROSYN-Clo1313_0385')  + v('PROSYN-Clo1313_0673') + v('PROSYN-Clo1313_1717') + v('PROSYN-Clo1313_1798') + v('PROSYN-Clo1313_2911');

glyc_min_obj..							z =e= -(v('PROSYN-Clo1313_1954') + v('PROSYN-Clo1313_0993') + v('PROSYN-Clo1313_0489') + v('PROSYN-Clo1313_1831') + v('PROSYN-Clo1313_2015') + v('PROSYN-Clo1313_0997') + v('PROSYN-Clo1313_1876') + v('PROSYN-Clo1313_0238') + v('PROSYN-Clo1313_1875') + v('PROSYN-Clo1313_2093') + v('PROSYN-Clo1313_2095') + v('PROSYN-Clo1313_2094') + v('PROSYN-Clo1313_0966') + v('PROSYN-Clo1313_2745') + v('PROSYN-Clo1313_2092') + v('PROSYN-Clo1313_2090') + v('PROSYN-Clo1313_0949') + v('PROSYN-Clo1313_0415') + v('PROSYN-Clo1313_1879') + v('PROSYN-Clo1313_0020') + v('PROSYN-Clo1313_0021') + v('PROSYN-Clo1313_0022') + v('PROSYN-Clo1313_0023') + v('PROSYN-Clo1313_0382') + v('PROSYN-Clo1313_0383') + v('PROSYN-Clo1313_0384') + v('PROSYN-Clo1313_0385')  + v('PROSYN-Clo1313_0673') + v('PROSYN-Clo1313_1717') + v('PROSYN-Clo1313_1798') + v('PROSYN-Clo1313_2911'));

glyc_mix_obj..							z =e= -(v('PROSYN-Clo1313_1954') + v('PROSYN-Clo1313_0993') + v('PROSYN-Clo1313_0489') + v('PROSYN-Clo1313_1831') + v('PROSYN-Clo1313_2015') + v('PROSYN-Clo1313_0997') + v('PROSYN-Clo1313_1876') + v('PROSYN-Clo1313_0238') + v('PROSYN-Clo1313_1875') + v('PROSYN-Clo1313_2093') + v('PROSYN-Clo1313_2095') + v('PROSYN-Clo1313_2094')) + (v('PROSYN-Clo1313_0966') + v('PROSYN-Clo1313_2745') + v('PROSYN-Clo1313_2092') + v('PROSYN-Clo1313_2090') + v('PROSYN-Clo1313_0949') + v('PROSYN-Clo1313_0415') + v('PROSYN-Clo1313_1879') + v('PROSYN-Clo1313_0020') + v('PROSYN-Clo1313_0021') + v('PROSYN-Clo1313_0022') + v('PROSYN-Clo1313_0023') + v('PROSYN-Clo1313_0382') + v('PROSYN-Clo1313_0383') + v('PROSYN-Clo1313_0384') + v('PROSYN-Clo1313_0385')  + v('PROSYN-Clo1313_0673') + v('PROSYN-Clo1313_1717') + v('PROSYN-Clo1313_1798') + v('PROSYN-Clo1313_2911'));

* require mass balance
mb_const(i)..							sum(j, S(i,j) * v(j)) =e= 0;

* define ribosome capacity
RiboCapacity..							v('RIBOSYN-ribo') * %kribo% =g= mu * sum(prosyn, NAA(prosyn) * v(prosyn));

*require synthesis of the cellulosome
enforce_cell_syn(prosyn_cell)..		v(prosyn_cell) =e= prosyn_cell_rates(prosyn_cell) * mu * %nscale%;

*enforce growth rate to be the value in the mu
enforce_mu..							v('BIOSYN-BIODIL') =e= mu * %nscale%;

* limit proteome allication to metabolic proteins
* this constrain is equivalent to the protein capacity constraint
NonModelProtAllo..						%nonmodeled_proteome_allocation% * v('BIOSYN-PROTCYT') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTDUMMY');

ribo_prot..					.03051 * mu * %nscale% =e= sum(prosyn_ribo, S('BIO-protcyt', prosyn_ribo) * v(prosyn_ribo));

* import the enzyme capacity constraints
$include %enz_cap_eqns_path%

*** BUILD OPTIMIZATION MODEL ***

Model RBAmaxETOH
/

	etoh_max_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

Model RBAmaxAC
/

	ac_max_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

Model RBAmaxETOHac
/

	etoh_ac_max_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

Model RBAmaxH2
/

	h2_max_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

Model RBAminH2
/

	h2_min_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

Model RBAmaxGLYC
/

	glyc_max_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

Model RBAminGLYC
/

	glyc_min_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

Model RBAmixGLYC
/

	glyc_mix_obj
	mb_const
	RiboCapacity
	enforce_cell_syn
	enforce_mu
	NonModelProtAllo
	
$include %enz_cap_declares_path%
	
/;

RBAmaxETOH.optfile = 1;
RBAmaxAC.optfile = 1;
RBAmaxETOHac.optfile = 1;
RBAmaxH2.optfile = 1;
RBAminH2.optfile = 1;
RBAmaxGLYC.optfile = 1;
RBAminGLYC.optfile = 1;
RBAmixGLYC.optfile = 1;

mu = 0.40438140;

SOLVE %usemodel% USING lp MAXIMIZING z;

*report flux rates for maximum product yield

FILE FLUXRESULTS / 'results/%identifier%/flux_%identifier%.txt' /;
PUT FLUXRESULTS;

FLUXRESULTS.lw = 45;
FLUXRESULTS.nr = 2;

PUT "rxn",system.tab,"LB",system.tab,"flux",system.tab,"UB"/;

LOOP(j,

	PUT j.tl,system.tab,(v.lo(j)*%nscaleback%),system.tab,(v.l(j)*%nscaleback%):0:8,system.tab,(v.up(j)*%nscaleback%)/;

);

PUTCLOSE;

*protein synthesis output file
FILE PROTRESULTS / 'results/%identifier%/prot_syn_%identifier%.txt' /;
PUT PROTRESULTS;

PROTRESULTS.lw = 25;
PROTRESULTS.nr = 2;

PUT "protein",system.tab,"synthesis rate"/;

loop(prosyn,

	PUT prosyn.tl,system.tab,(v.l(prosyn)*%nscaleback%):0:8/;

);

PUTCLOSE;

*protein synthesis output file
FILE PROTPERUP / 'results/%identifier%/per_up_%identifier%.txt' /;
PUT PROTPERUP;

PROTPERUP.lw = 25;
PROTPERUP.nr = 2;

PUT "protein",system.tab,"synthesis rate (nmol protein/mmol c6 unit)"/;

loop(prosyn,

	calc_temp = ((v.l(prosyn)*%nscaleback%)/%c6_up_rate%) * 1E6;

	PUT prosyn.tl,system.tab,calc_temp:0:8/;

);

calc_temp = sum(prosyn, v.l(prosyn)*%nscaleback% * 1E6);
PUT "Sum",system.tab,calc_temp:0:8/;

PUTCLOSE;

*carbon products out
FILE CARBOUT / 'results/%identifier%/carb_out_%identifier%.txt' /;
PUT CARBOUT;

CARBOUT.lw = 25;
CARBOUT.nr = 2;

PUT "carbon exchange",system.tab,"output rate (mmol species/mmol c6 unit)"/;

loop(carbs,

	calc_temp = ((v.l(carbs)*%nscaleback%)/%c6_up_rate%);

	PUT carbs.tl,system.tab,calc_temp:0:8/;

);

PUTCLOSE;
