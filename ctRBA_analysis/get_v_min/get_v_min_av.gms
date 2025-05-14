************************* Run RBA parameterization step of minimizing inactive reactions used **************************
*       Origianl Author: Hoang Dinh (for scRBA model)																																	 *
*       Latest Author: Wheaton Schroeder (for ctherm RBA model)																												 *
*				Latest Version:02/15/2024																																										 *
* 			Written as part of the workflow to determine Kapp values, this step is to determine which inactive reactions 	 *
*				are essential							 																																										 *
************************************************************************************************************************

$INLINECOM /*  */
$include "./get_v_min_av_settings.txt"

* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1E3
$setGlobal nscaleback 1E-3

* For second run when the protein capacity is corrected, overwrite this default value of 1
* Value of 1 means the protein capacity is used freely, which leads to cheap proteins that
* catalyzed flux cycle to be produced.
* Flux cycle consists of the forward and reverse directions of a reversible reaction
* catalyzed by a protein that is cheaper than the dummy protein. Why is it cheaper?
* It is because the composition of that protein has higher fractions of cheap
* amino acids compared to the dummy protein which has the composition of the
* experimentally observations.
* This happens to system where protein availability
* is not the limiting factor (e.g., S. cerevisiae where nutrient or rRNA availability is
* the limiting factor).
$setGlobal corrected_protein_capacity_percentage 1

OPTIONS

	LP = soplex 				/*Solver selection*/
	limrow = 0 					/*number of equations listed, 0 is suppresed*/
	limcol = 0 					/*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 	/*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 				/*decimal places for display statement*/
	reslim = 1000000 		/*wall-clock time limit for solver in seconds*/
	sysout = on 				/*solver status file report option*/
	solprint = on 			/*solution printing option*/

;

SETS

	i										/*list of all metabolic species*/
$include "%species_path%"
	
	j										/*list of all reactions*/
$include "%rxns_path%"

	j_inactive(j)				/*list of inactive reactions*/
$include "%inactive_path%"

	prosyn(j)						/*list of protein synthesis reactions*/
$include "%prosyn_path%"

	prowaste(j)					/*list of protein waste reactions*/
$include "%prowaste_path%"

	translation(j)			/*list of translation (protein synthesis) reactions same as prosyn may have been needed due to multiple ribosomes in scRBA model may be removable later*/
$include "%trans_path%"

	uptake(j) 					/*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"

	media(j) 						/*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"

	carbs(j)						/*list of reactions where the carbon yield is constrained by MFA values*/
$include "%car_frac_j_path%"

	prosyn_met(prosyn)					/*list of reactions where protein translation is constrained by A-N measurements values*/
$include "%prot_met_rxns%"

	prosyn_cell(prosyn)					/*list of reactions where protein translation is constrained by A-N measurements values*/
$include "%prot_cell_rxns%"

	prosyn_ribo(prosyn)					/* list of protein translation reactions that make ribosomal proteins */
$include "%prot_ribo_rxns%"

	blocked(j)					/*list of reactions blocked in the stoichiometric model that should also be in the RAM model*/
$include "%blocked_path%"

	inactive(j)
$include "%inactive_path%"

;

PARAMETERS

	S(i,j)								/*stoichiometric matrix*/
$include "%sij_path%"
	
	NAA(prosyn)							/*number of amino acids in protein synthesis reactions*/
$include "%prolen_path%"

	c_min_f(j)							/*lower bound limits for flux rates taken from MFA tested already in the base model*/
$include "%carb_min_path%"

	c_max_f(j)							/*upper bound limits for flux rates taken from MFA tested already in the base model*/
$include "%carb_max_path%"

	prosyn_met_rates(prosyn_met)					/*rate of protein synthesis to impose on protein synthesis reactions estimated from proteomics and growth rate*/
$include "%prot_met_rates%"

	prosyn_cell_rates(prosyn_cell)					/*rate of protein synthesis to impose on protein synthesis reactions estimated from proteomics and growth rate*/
$include "%prot_cell_rates%"

calc_temp								/* holder for temporary calculations */
;

VARIABLES

	z												/*objective variable*/

	v(j)										/*reaction rate variable*/

	S_plus(prosyn)
	S_minus(prosyn)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; 
v.up(j) = 1e3 * %nscale%;

S_plus.lo(prosyn) = 0;
S_minus.lo(prosyn) = 0;
S_plus.up(prosyn) = 0* %nscale%;
S_minus.up(prosyn) = 0 * %nscale%;

* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
*v.fx(j)$prowaste(j) = 0;
v.lo(j)$prowaste(j) = 0;
v.up(j)$prowaste(j) = 1E3 * %nscale%;

* Disable all biomass reactions
* Condition-specific biomass reaction activation has to be done in phenotype.txt file
v.up('BIOSYN-BIODIL') = 0; 
v.lo('BIOSYN-BIODIL') = 0; 

* Media
* disables all uptakes that are not in media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1E3 * %nscale%;

* in phenotype is set:
* substrate
* NGAM
* Biomass composition
* Changed this bit so I can specify phenotype in settings
$include "%phenotype%"

* Set your NGAM in phenotype.txt since NGAM value depends on growth-condition
* Set your biomass composition in phenotype.txt since biomass composition depends on growth-condition

*disable inactive reactions
v.lo(j_inactive) = 0;
v.up(j_inactive) = 0;

*keep blocked reactions as blocked
v.fx(blocked) = 0;

*** EQUATION DEFINITIONS ***

Equations

	Obj													/*Objective equation*/
	Stoic												/*mass balance equation*/
	RiboCapacity								/*ribosome capacity constraint based on number of amino acids*/
	NonModelProtAllo						/*require dummy protein synthesis*/
	shunt_ratio									/*enforces roughly 2:1 flux split between ppdk and malate shunt for pep to pyr conversion*/ 
	ribo_prot

;

* minimize flux sum across the model
*note that this is supposed to be the sum of active reactions but since active + inactive = all reactions by fixing inactive reactions to zero, this is essentially the same thing
Obj..											z =e= sum(j,v(j));

* require mass balance
Stoic(i)..								sum(j, S(i,j)*v(j)) =e= 0;

* define ribosome capacity used to enforce ribosome synthesis
RiboCapacity.. 						v('RIBOSYN-ribo') * %kribo% =g= %mu% * sum(prosyn, NAA(prosyn) * v(prosyn));

* require synthesis of non-modeled protein
NonModelProtAllo..			%nonmodeled_proteome_allocation% * v('BIOSYN-PROTCYT') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTDUMMY');

* require 2:1 flux ratio of PPDK to malate shunt allows estimation on all parts of the shunt
shunt_ratio..						2 * v('RXN-ME2_FWD-maeB') =e= v('RXN-PPDK_FWD-ppdk');

ribo_prot..					.02938 * %mu% * %nscale% =e= sum(prosyn_ribo, S('BIO-protcyt', prosyn_ribo) * v(prosyn_ribo));


*** BUILD OPTIMIZATION MODEL ***

Model GETVMIN
/

	Obj
	Stoic
	RiboCapacity
        NonModelProtAllo
	shunt_ratio
	ribo_prot
/;

GETVMIN.optfile = 1;

*** SOLVE ***
SOLVE GETVMIN USING lp MINIMIZING z;

*write a file on the model status
FILE statFile / 'results/av/get_v_min_modelStat.txt' /;
PUT statFile;

PUT "modelstat: ",GETVMIN.modelStat/;
PUT "solvestat: ",GETVMIN.solveStat/;

PUTCLOSE statFile;

*write a file for flux results just those with flux
file fluxFileMin / 'results/av/get_vmin_flux_min.txt' /;
PUT fluxFileMin;

PUT 'RXN',system.tab,'LB',system.tab,'flux',system.tab,'UB'/;

loop(j,

	/*note: seems to ignore unused reactions in the output*/
	if ( (v.l(j) gt 0),
	
		PUT j.tl:0, system.tab, (v.lo(j) * %nscaleback%), system.tab, (v.l(j) * %nscaleback%):0:8, system.tab,(v.up(j) * %nscaleback%)/;
		
	);
	
);

PUTCLOSE fluxFileMin;

*write a file for flux results all reactions including zeros
file fluxFileAll / 'results/av/get_v_min_flux_all.txt' /;
PUT fluxFileAll;

fluxFileAll.lw = 25;
fluxFileAlL.nr = 2;

PUT 'RXN',system.tab,'LB',system.tab,'flux',system.tab,'UB'/;

loop(j,

	PUT j.tl:0, system.tab, (v.lo(j) * %nscaleback%), system.tab, (v.l(j) * %nscaleback%):0:8, system.tab,(v.up(j) * %nscaleback%)/;
	
);

PUTCLOSE fluxFileAll;

FILE PROTPERUP / 'results/av/per_up_v_min.txt' /;
PUT PROTPERUP;

PROTPERUP.lw = 25;
PROTPERUP.nr = 2;

PUT "protein",system.tab,"synthesis rate (nmol protein/mmol c6 unit)"/;

loop(prosyn,

	calc_temp = ((v.l(prosyn)*%nscaleback%)/ 8.6) * 1E6;

	PUT prosyn.tl,system.tab,calc_temp:0:8/;

);

calc_temp = sum(prosyn, v.l(prosyn)*%nscaleback% *1E6);
PUT "Sum",system.tab,calc_temp /;

PUTCLOSE;
