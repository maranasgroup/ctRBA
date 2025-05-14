#!/usr/bin/python
#! python 3.9
#try to specify that we will use python version 3.9
__author__ = "Wheaton Schroeder"
#latest version: 10/05/2023

#imports, again, these are used previously, not sure of their necessity
import cobra
from fba import FBA
from datetime import datetime
from cobra import Model, Reaction, Metabolite
import re
import os
import warnings

#get the current time and date for tracking solution time
start_time = datetime.now()
print("Starting time: ",start_time)

#get the current directory to use for importing things
curr_dir = os.getcwd()

#import the ctherm model
model = cobra.io.read_sbml_model(curr_dir + "/iCTH669_for_MFA.sbml")

#write a file to put the FBA results
output=open('pFBA_results_ctherm669_for_RBA_MFA_data_clean.txt','w',buffering=1)

#fixed rates dictionary, if we want it
#these make fixed versions of the strains, rather than knockouts so same list of reactions still
fixed_rates = {

    #ensure overflow metabolism by requiring full cellbiose uptake
    "EXCH_glc__D_e":-100,
    "EXCH_cellb_e":0,

}

#load fba
object = FBA(model)

#objective = "BIOMASS"
objective = "MFA_BIOMASS"
#objective = "DM_TEST"
#objective = "PYRt"
dir = "max"

############################ DEFINE CONSTRAINTS FOR APPLYING MFA DATA ############################

################ EXCHANGES ################

######## ETHANOL ########

#ethanol export
#define bounds
object.model.reactions.get_by_id("EXCH_etoh_e").upper_bound = 24.7319
object.model.reactions.get_by_id("EXCH_etoh_e").lower_bound = 20.655

######## ACETATE ########

#acetate export
#define bounds
object.model.reactions.get_by_id("EXCH_ac_e").upper_bound = 34.3793
object.model.reactions.get_by_id("EXCH_ac_e").lower_bound = 29.9675

######## CO2 ########

#CO2 export
#define bounds
object.model.reactions.get_by_id("EXCH_co2_e").upper_bound = 79.6414
object.model.reactions.get_by_id("EXCH_co2_e").lower_bound = 74.0563

######## PYRUVATE ########

#CO2 export
#define bounds
object.model.reactions.get_by_id("EXCH_pyr_e").upper_bound = 6.26
object.model.reactions.get_by_id("EXCH_pyr_e").lower_bound = 5.1795

######## LACTATE ########

#lactate export
#define bounds
object.model.reactions.get_by_id("EXCH_lac__L_e").upper_bound = 8.0832
object.model.reactions.get_by_id("EXCH_lac__L_e").lower_bound = 6.7279

######## FORMATE ########

#formate export
#define bounds
object.model.reactions.get_by_id("EXCH_for_e").upper_bound = 22.1002
object.model.reactions.get_by_id("EXCH_for_e").lower_bound = 19.2576

######## MALATE ########

#malate export
#define bounds
object.model.reactions.get_by_id("EXCH_mal__L_e").upper_bound = 7.1676
object.model.reactions.get_by_id("EXCH_mal__L_e").lower_bound = 0.0000010174

######## HYDROGEN GAS ########

#hydrogen gas export
#define bounds
object.model.reactions.get_by_id("EXCH_h2_e").upper_bound = 129.1105
object.model.reactions.get_by_id("EXCH_h2_e").lower_bound = 118.9361

######## L-ALANINE ########

#L-alanine export
#define bounds
object.model.reactions.get_by_id("EXCH_ala__L_e").upper_bound = 1.2016
object.model.reactions.get_by_id("EXCH_ala__L_e").lower_bound = 0.9946

######## L-VALINE ########

#L-valine export
#define bounds
object.model.reactions.get_by_id("EXCH_val__L_e").upper_bound = 5.5585
object.model.reactions.get_by_id("EXCH_val__L_e").lower_bound = 4.6205

######## L-ASPARTATE ########

#L-aspartate export
#define bounds
object.model.reactions.get_by_id("EXCH_asp__L_e").upper_bound = 6.7393
object.model.reactions.get_by_id("EXCH_asp__L_e").lower_bound = 0

######## L-THREONINE ########

#L-threonine export
#define bounds
object.model.reactions.get_by_id("EXCH_thr__L_e").upper_bound = 0.925
object.model.reactions.get_by_id("EXCH_thr__L_e").lower_bound = 0

######## FUMARATE ########

#fumarate export
#define bounds
object.model.reactions.get_by_id("EXCH_fum_e").upper_bound = 2.6678
object.model.reactions.get_by_id("EXCH_fum_e").lower_bound = 0.000053212

######## L-LEUCINE ########

#L-leucine export
#define bounds
object.model.reactions.get_by_id("EXCH_leu__L_e").upper_bound = 1.193
object.model.reactions.get_by_id("EXCH_leu__L_e").lower_bound = 0

######## L-GLUTAMATE ########

#L-glutamate export
#define bounds
object.model.reactions.get_by_id("EXCH_glu__L_e").upper_bound = 2.3354
object.model.reactions.get_by_id("EXCH_glu__L_e").lower_bound = 0

######## ALPHA-KETOGLUTERATE ########

#alpha-ketogluterate export
#define bounds
object.model.reactions.get_by_id("EXCH_akg_e").upper_bound = 6.2633
object.model.reactions.get_by_id("EXCH_akg_e").lower_bound = 0.00000088218

######## SUCCINATE ########

#succinate export
#alpha-ketogluterate export
#define bounds
object.model.reactions.get_by_id("EXCH_succ_e").upper_bound = 10.5402
object.model.reactions.get_by_id("EXCH_succ_e").lower_bound = 0

######## L-ASPARAGINE ########

#L-asparagine export
#define bounds
object.model.reactions.get_by_id("EXCH_asn__L_e").upper_bound = 0.6049
object.model.reactions.get_by_id("EXCH_asn__L_e").lower_bound = 0.4974



################ GLYCOLYSIS ################

######## PGI ########

#define bounds
object.model.reactions.get_by_id("PGI").upper_bound = 98.6759
object.model.reactions.get_by_id("PGI").lower_bound = 98.6358

######## F6PPT ########

#define bounds
object.model.reactions.get_by_id("F6PPT").upper_bound = 93.1173
object.model.reactions.get_by_id("F6PPT").lower_bound = 92.6655


######## FBA ########

#define bounds
object.model.reactions.get_by_id("FBA").upper_bound = 93.1173
object.model.reactions.get_by_id("FBA").lower_bound = 92.6655

"""
######## TPI ########

#define bounds
object.model.reactions.get_by_id("TPI").upper_bound = 90.7324
object.model.reactions.get_by_id("TPI").lower_bound = 90.3885

"""
"""
######## GAPD ########

#define bounds
object.model.reactions.get_by_id("GAPD").upper_bound = 175.7842
object.model.reactions.get_by_id("GAPD").lower_bound = 174.8892
"""
"""
######## PGK & PGK2 ########

#define constraint, overall conversion
pgk_conv_const = object.model.problem.Constraint(getattr(object.model.reactions,"PGK").flux_expression + getattr(object.model.reactions,"PGK2").flux_expression, lb=174.8892, ub=175.7842)

#add constraint
object.model.add_cons_vars(pgk_conv_const)
"""

######## PGM & ENO ########

#define bounds
object.model.reactions.get_by_id("PGM").upper_bound = 168.5545
object.model.reactions.get_by_id("PGM").lower_bound = 156.888

#define bounds
object.model.reactions.get_by_id("ENO").upper_bound = 168.5545
object.model.reactions.get_by_id("ENO").lower_bound = 156.888

######## MALATE SHUNT AND BYPASS ########

#define constraint, overall conversion
MS_conv_const = object.model.problem.Constraint(getattr(object.model.reactions,"PPDK").flux_expression + getattr(object.model.reactions,"ME2").flux_expression, lb=132.8498, ub=136.7328)

#add constraint
object.model.add_cons_vars(MS_conv_const)


"""
################ TCA CYCLE ################

######## PEPCK ########

#define bounds
#stops growth
object.model.reactions.get_by_id("PEPCK_re").upper_bound = 21.4218
object.model.reactions.get_by_id("PEPCK_re").lower_bound = 18.0848
"""
"""
######## MDH ########

#define bounds
#stops growth
object.model.reactions.get_by_id("MDH").upper_bound = 14.4583
object.model.reactions.get_by_id("MDH").lower_bound = -3.45
"""

######## CS ########

cs_re_flux = -8.83
cs_re_lb = -12.5
cs_re_ub = 12.4
cs_si_flux = 15.9
cs_si_lb = -4.46
cs_si_ub = 20.1

#define bounds
object.model.reactions.get_by_id("CS").upper_bound = cs_re_ub + cs_si_ub
object.model.reactions.get_by_id("CS").lower_bound = cs_re_lb + cs_si_lb


######## ACONT, ICOR, OXSCL ########

#define bounds
object.model.reactions.get_by_id("ACONT").upper_bound = 9.552
object.model.reactions.get_by_id("ACONT").lower_bound = 6.7599

#define bounds
object.model.reactions.get_by_id("ICOR").upper_bound = 9.552
object.model.reactions.get_by_id("ICOR").lower_bound = 6.7599

#define bounds
object.model.reactions.get_by_id("OXSCL").upper_bound = 9.552
object.model.reactions.get_by_id("OXSCL").lower_bound = 6.7599

######## FUM ########

#define bounds
object.model.reactions.get_by_id("FUM").upper_bound = 6.6755
object.model.reactions.get_by_id("FUM").lower_bound = -3.4519

######## FRDx ########

#define bounds
object.model.reactions.get_by_id("FRDx").upper_bound = 0
object.model.reactions.get_by_id("FRDx").lower_bound = -10.5402


################ ONE CARBON METABOLISM ################

######## FTHFLi ########

#define bounds
object.model.reactions.get_by_id("FTHFLi").upper_bound = 0.2144
object.model.reactions.get_by_id("FTHFLi").lower_bound = 0.1961

######## MTHFD ########

#define bounds
object.model.reactions.get_by_id("MTHFD").upper_bound = 2.5305
object.model.reactions.get_by_id("MTHFD").lower_bound = 2.3151

######## MTHFO_nadp ########

#define bounds
object.model.reactions.get_by_id("MTHFO_nadp").upper_bound = 0.8959
object.model.reactions.get_by_id("MTHFO_nadp").lower_bound = 0.8196

################ ELECTRON METABOLISM ################

######## BIF ########

#define bounds
object.model.reactions.get_by_id("BIF").lower_bound = -70.6
object.model.reactions.get_by_id("BIF").upper_bound = -41.4

######## ECH ########

#define bounds
object.model.reactions.get_by_id("ECH").upper_bound = 250.1363
object.model.reactions.get_by_id("ECH").lower_bound = 209.0912

######## FRNDPR2r ########

#define bounds
object.model.reactions.get_by_id("FRNDPR2r").upper_bound = 22.2618
object.model.reactions.get_by_id("FRNDPR2r").lower_bound = -7.448

######## RNF ########

#define bounds
#reduces growth 11.5 -> 11.02
object.model.reactions.get_by_id("RNF").upper_bound = -122.3798
object.model.reactions.get_by_id("RNF").lower_bound = -134.4827

################ AMINO ACID METABOLISM ################

######## GHMT2r ########
"""
#define bounds
object.model.reactions.get_by_id("GHMT2r").upper_bound = -3.1347
object.model.reactions.get_by_id("GHMT2r").lower_bound = -3.4264
"""
######## GLU5K; G5SD; G5SADs; P5CR ########

#define bounds
object.model.reactions.get_by_id("GLU5K").upper_bound = 1.2101
object.model.reactions.get_by_id("GLU5K").lower_bound = 1.1071

#define bounds
object.model.reactions.get_by_id("G5SD").upper_bound = 1.2101
object.model.reactions.get_by_id("G5SD").lower_bound = 1.1071

#define bounds
object.model.reactions.get_by_id("G5SADs").upper_bound = 1.2101
object.model.reactions.get_by_id("G5SADs").lower_bound = 1.1071

#define bounds
object.model.reactions.get_by_id("P5CR").upper_bound = 1.2101
object.model.reactions.get_by_id("P5CR").lower_bound = 1.1071

"""
######## ASPTA ########

#define bounds
object.model.reactions.get_by_id("ASPTA").upper_bound = 15.2898
object.model.reactions.get_by_id("ASPTA").lower_bound = 14.2351

"""

######## ASNS2 ########

#define bounds
#stops growth
object.model.reactions.get_by_id("ASNS2").upper_bound = 2.4638
object.model.reactions.get_by_id("ASNS2").lower_bound = 2.3564

######## GLYTA & AGT ########

#define bounds
#stops growth
object.model.reactions.get_by_id("GLYTA").upper_bound = -6.1196
object.model.reactions.get_by_id("GLYTA").lower_bound = -6.4412

#define bounds
#stops growth
object.model.reactions.get_by_id("AGT").upper_bound = -6.1196
object.model.reactions.get_by_id("AGT").lower_bound = -6.4412

######## CYS ########

#define bounds
object.model.reactions.get_by_id("CYS").upper_bound = 1.3141
object.model.reactions.get_by_id("CYS").lower_bound = 1.2023


######## HSDH, HSK, & THRS ########

#define bounds
object.model.reactions.get_by_id("HSDH").upper_bound = 4.0052
object.model.reactions.get_by_id("HSDH").lower_bound = 1.6423

#define bounds
object.model.reactions.get_by_id("HSK").upper_bound = 4.0052
object.model.reactions.get_by_id("HSK").lower_bound = 1.6423

#define bounds
object.model.reactions.get_by_id("THRS").upper_bound = 4.0052
object.model.reactions.get_by_id("THRS").lower_bound = 1.6423

######## VALTA ########

#define bounds
object.model.reactions.get_by_id("VALTA").upper_bound = -7.0663
object.model.reactions.get_by_id("VALTA").lower_bound = -8.0209

######## METS ########

#define bounds
object.model.reactions.get_by_id("METS").upper_bound = 0.8959
object.model.reactions.get_by_id("METS").lower_bound = 0.8196

######## ACAS_2ahbut, KARI, KARI_23dhmp, DHAD2, ILETA ########

#define bounds
object.model.reactions.get_by_id("ACAS_2ahbut").upper_bound = 4.3122
object.model.reactions.get_by_id("ACAS_2ahbut").lower_bound = 2.7634

#define bounds
object.model.reactions.get_by_id("KARI").upper_bound = 4.3122
object.model.reactions.get_by_id("KARI").lower_bound = 2.7634

#define bounds
object.model.reactions.get_by_id("KARI_23dhmp").upper_bound = 4.3122
object.model.reactions.get_by_id("KARI_23dhmp").lower_bound = 2.7634

#define bounds
object.model.reactions.get_by_id("DHAD2").upper_bound = 4.3122
object.model.reactions.get_by_id("DHAD2").lower_bound = 2.7634

#define bounds
object.model.reactions.get_by_id("ILETA").upper_bound = 4.3122
object.model.reactions.get_by_id("ILETA").lower_bound = 2.7634

######## IPPMIa, IPMD, OMCDC, & LEUTA ########

#define bounds
object.model.reactions.get_by_id("IPPS").upper_bound = 3.8255
object.model.reactions.get_by_id("IPPS").lower_bound = 2.8231

#define bounds
object.model.reactions.get_by_id("IPPMIb").upper_bound = 3.8255
object.model.reactions.get_by_id("IPPMIb").lower_bound = 2.8231

#define bounds
object.model.reactions.get_by_id("IPPMIa").upper_bound = -2.8231
object.model.reactions.get_by_id("IPPMIa").lower_bound = -3.8255

#define bounds
object.model.reactions.get_by_id("IPMD").upper_bound = 3.8255
object.model.reactions.get_by_id("IPMD").lower_bound = 2.8231

#define bounds
object.model.reactions.get_by_id("OMCDC").upper_bound = 3.8255
object.model.reactions.get_by_id("OMCDC").lower_bound = 2.8231

#define bounds
#reducesa growth rates 11.5 -> 10.93
object.model.reactions.get_by_id("LEUTA").upper_bound = 3.8255
object.model.reactions.get_by_id("LEUTA").lower_bound = 2.8231

######## PGLYDH, PSAT, & PSP ########

#define bounds
object.model.reactions.get_by_id("PGLYDH").upper_bound = 18.8877
object.model.reactions.get_by_id("PGLYDH").lower_bound = 8.2423

#define bounds
object.model.reactions.get_by_id("PSAT").upper_bound = -8.2423
object.model.reactions.get_by_id("PSAT").lower_bound = -18.8877

#define bounds
object.model.reactions.get_by_id("PSP").upper_bound = 18.8877
object.model.reactions.get_by_id("PSP").lower_bound = 8.2423

######## ANS, ANPRT, PRAI, IGPS, & SERH ########

#define bounds
object.model.reactions.get_by_id("SERH").upper_bound = 0.3163
object.model.reactions.get_by_id("SERH").lower_bound = 0.2894

#define bounds
object.model.reactions.get_by_id("ANS").upper_bound = 0.3163
object.model.reactions.get_by_id("ANS").lower_bound = 0.2894

#define bounds
object.model.reactions.get_by_id("ANPRT").upper_bound = 0.3163
object.model.reactions.get_by_id("ANPRT").lower_bound = 0.2894

#define bounds
object.model.reactions.get_by_id("PRAI").upper_bound = 0.3163
object.model.reactions.get_by_id("PRAI").lower_bound = 0.2894

#define bounds
object.model.reactions.get_by_id("IGPS").upper_bound = 0.3163
object.model.reactions.get_by_id("IGPS").lower_bound = 0.2894

######## DAPDC ########

#define bounds
object.model.reactions.get_by_id("DAPDC").upper_bound = 2.9092
object.model.reactions.get_by_id("DAPDC").lower_bound = 2.6615

######## LHOR, HISOR, HISTP, HSTPT, IGPDH, IG3PS, PRMICI, PRATPP, & ATPPRT ########

#define bounds
object.model.reactions.get_by_id("LHOR").upper_bound = 0.5101
object.model.reactions.get_by_id("LHOR").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("HISOR").upper_bound = 0.5101
object.model.reactions.get_by_id("HISOR").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("HISTP").upper_bound = 0.5101
object.model.reactions.get_by_id("HISTP").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("HSTPT").upper_bound = 0.5101
object.model.reactions.get_by_id("HSTPT").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("IGPDH").upper_bound = 0.5101
object.model.reactions.get_by_id("IGPDH").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("IG3PS").upper_bound = 0.5101
object.model.reactions.get_by_id("IG3PS").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("PRMICI").upper_bound = 0.5101
object.model.reactions.get_by_id("PRMICI").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("PRAMPC").upper_bound = 0.5101
object.model.reactions.get_by_id("PRAMPC").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("PRATPP").upper_bound = 0.5101
object.model.reactions.get_by_id("PRATPP").lower_bound = 0.4666

#define bounds
object.model.reactions.get_by_id("ATPPRT").upper_bound = 0.5101
object.model.reactions.get_by_id("ATPPRT").lower_bound = 0.4666

######## ORNTAC, ACOTA, AGPR, & ACGK ########

#define bounds
object.model.reactions.get_by_id("ORNTAC").upper_bound = 4.4693
object.model.reactions.get_by_id("ORNTAC").lower_bound = 0

#define bounds
object.model.reactions.get_by_id("ACOTA").upper_bound = 0
object.model.reactions.get_by_id("ACOTA").lower_bound = -4.4693

#define bounds
object.model.reactions.get_by_id("AGPR").upper_bound = 0
object.model.reactions.get_by_id("AGPR").lower_bound = -4.4693

#define bounds
object.model.reactions.get_by_id("ACGK").upper_bound = 4.4693
object.model.reactions.get_by_id("ACGK").lower_bound = 0


######## PHETA1 & PPNDH ########
#define bounds
object.model.reactions.get_by_id("PHETA1").upper_bound = -1.384
object.model.reactions.get_by_id("PHETA1").lower_bound = -1.5128

#define bounds
object.model.reactions.get_by_id("PPNDH").upper_bound = 1.5128
object.model.reactions.get_by_id("PPNDH").lower_bound = 1.384


######## TYRTA & PPND ########
#define bounds
object.model.reactions.get_by_id("TYRTA").upper_bound = -1.38
object.model.reactions.get_by_id("TYRTA").lower_bound = -1.5128

#define bounds
object.model.reactions.get_by_id("PPND").upper_bound = 1.5128
object.model.reactions.get_by_id("PPND").lower_bound = 1.38

"""
######## GLUN ########
#define bounds
#decreases growth rate 11.5 -> 6.58
object.model.reactions.get_by_id("GLUN").upper_bound = -(64.5063+61.8825)
object.model.reactions.get_by_id("GLUN").lower_bound = -(69.1581+66.0177)

#define bounds
object.model.reactions.get_by_id("GLUDy").upper_bound = -61.8825
object.model.reactions.get_by_id("GLUDy").lower_bound = -66.0177
"""

################ PENTOSE PHOSPHATE PATHWAY ################

######## FBA3 ########

#define bounds
object.model.reactions.get_by_id("FBA3").upper_bound = 0.1154
object.model.reactions.get_by_id("FBA3").lower_bound = 0.1055

######## PFK3_ppi ########

#define bounds
object.model.reactions.get_by_id("PFK3_ppi").upper_bound = 0.1154
object.model.reactions.get_by_id("PFK3_ppi").lower_bound = 0.1055

######## RPE ########
#define bounds
object.model.reactions.get_by_id("RPE").upper_bound = -2.8424
object.model.reactions.get_by_id("RPE").lower_bound = -3.1069

######## RPI ########
#define bounds
object.model.reactions.get_by_id("RPI").upper_bound = -2.8351
object.model.reactions.get_by_id("RPI").lower_bound = -3.0989

######## ACALD ########

#define bounds
object.model.reactions.get_by_id("ACALD").upper_bound = 24.7319
object.model.reactions.get_by_id("ACALD").lower_bound = 20.655

######## ALCD2x ########

#define bounds
object.model.reactions.get_by_id("ALCD2x").upper_bound = 24.7319
object.model.reactions.get_by_id("ALCD2x").lower_bound = 20.655

######## LDL_L ########

#define bounds
object.model.reactions.get_by_id("LDH_L").upper_bound = 8.0832
object.model.reactions.get_by_id("LDH_L").lower_bound = 6.7279

######## PFL ########

#define bounds
object.model.reactions.get_by_id("PFL").upper_bound = 22.2495
object.model.reactions.get_by_id("PFL").lower_bound = 20.0477

######## POR ########

#define bounds
object.model.reactions.get_by_id("POR").upper_bound = 67.7686
object.model.reactions.get_by_id("POR").lower_bound = 62.0935

######## PTAr ########

#define bounds
object.model.reactions.get_by_id("PTAr").upper_bound = 34.2002
object.model.reactions.get_by_id("PTAr").lower_bound = 29.945

######## ACKr ########

#define bounds
object.model.reactions.get_by_id("ACKr").upper_bound = 34.2002
object.model.reactions.get_by_id("ACKr").lower_bound = 29.945


######## G3PD1ir ########

#define bounds
object.model.reactions.get_by_id("G3PD1ir").upper_bound = 12.8296
object.model.reactions.get_by_id("G3PD1ir").lower_bound = 0

######## BIOMASS ########

#define bounds
object.model.reactions.get_by_id("MFA_BIOMASS").upper_bound = 8.8263
object.model.reactions.get_by_id("MFA_BIOMASS").lower_bound = 8.0749



#################################### FINISHED MFA CONSTRAINTS #####################################

#note that for some reason the "R_" gets dropped
#run regular FBA, for determining viability
#results = object.run(objective,dir,fixed_rates)

#run pFBA to get flux rates
results = object.run_pFBA(objective,dir,fixed_rates)

#check if an exception occured
if results['exception']:

    #if here, an exception occured, we have no solution
    output.write("time to solve: "+str(results['total_time'])+"\n")
    output.write("exception occured, no solution exception: \n"+results['exception_str']+"\n")

    print("FAILED... exception was thrown:\n",results['exception_str'])

else:
    
    #create a formatted string to report a large number of decimal places
    formatted_string = "{:.8f}".format(results['objective'])

    print("objective value (flux sum): ",formatted_string)
    print("objective value ("+objective+"): "+str(results[objective]['flux']))
    print("model status 1: "+str(results['status1']))
    print("model status 2: "+str(results['status2']))

    #write header
    output.write("reaction\tlb\tlevel\tub\n")

    #also create an output file for flux mapping for escher plot
    escher_flux_f = open("escher_flux_data.txt","w")

    #write a header. I think the first line gets ignored
    escher_flux_f.write("reaction,rate\n")

    #for each reaction in the model
    for rxn in model.reactions:

        #write the FBA results
        output.write(rxn.id+"\t"+str(results[rxn.id]['lb'])+"\t"+str(results[rxn.id]['flux'])+"\t"+str(results[rxn.id]['ub'])+"\n")
        escher_flux_f.write(rxn.id+","+str(results[rxn.id]['flux'])+"\n")