#!/usr/bin/python

################################################
# FluctuProtST, Version 1.2
# Francois Bertaux, Inria Paris-Rocquencourt
# francois.bertaux@inria.fr
# March 2015
################################################


###########################################################################
# Description of the 'hEARM' model of TRAIL-induced apoptosis, proposed in
# Bertaux et al., PLoS Computational Biology, 2014.
#
# It extends the EARM model used in Spencer et al., Nature, 2009,
# by accounting for protein fluctuations.
###########################################################################

####
# MODS ! cPARP stab , non-fitted , no TRAIL back or removed 
####

#### imports
import FluctuProtST as fpst

### define native (i.e. fluctuating) proteins
def defineNativeProteins (model) :
	# direct path to death
	model.addNativeProteinStdFluct ( name="R" , EP=1000. , dilutionHalfLife=27. )
	model.addNativeProteinStdFluct ( name="pC8" , EP=10000. , dilutionHalfLife=27. )
	model.addNativeProtein ( name="pC3" , EP=10000. , HLP=27. , EM=17. , HLM=9. , Ton=0.099 , Toff=3.46 )
	model.addNativeProteinStdFluct ( name="pC6" , EP=10000. , dilutionHalfLife=27. )
	model.addNativeProteinStdFluct ( name="PARP" , EP=1000000. , dilutionHalfLife=27. )
	# direct path to death : inhibitors
	model.addNativeProtein ( name="flip" , EP=2000. , HLP=0.5 , EM=17. , HLM=2. , Ton=0.1 , Toff=2.6 )
	model.addNativeProteinStdFluct ( name="Bar" , EP=1000. , dilutionHalfLife=27. )
	model.addNativeProtein ( name="XIAP" , EP=100000. , HLP=27. , EM=17. , HLM=9. , Ton=0.098 , Toff=3.64 )
	# indirect path to death : mitchondrial pathway
	model.addNativeProtein ( name="Bid" , EP=60000. , HLP=27. , EM=17. , HLM=9. , Ton=0.098 , Toff=3.64 )
	model.addNativeProtein ( name="Bax" , EP=80000. , HLP=27. , EM=17. , HLM=9. , Ton=0.099 , Toff=3.14 )
	model.addNativeProteinStdFluct ( name="Pore" , EP=500000. , dilutionHalfLife=27. )
	model.addNativeProteinStdFluct ( name="smac_m" , EP=100000. , dilutionHalfLife=27. )
	model.addNativeProteinStdFluct ( name="CyC_m" , EP=500000. , dilutionHalfLife=27. )
	model.addNativeProteinStdFluct ( name="pC9" , EP=100000. , dilutionHalfLife=27. )
	model.addNativeProteinStdFluct ( name="Apaf" , EP=100000. , dilutionHalfLife=27. )
	# indirect path to death : inhibitors
	model.addNativeProtein ( name="Mcl1" , EP=20000. , HLP=0.5 , EM=17. , HLM=2. , Ton=0.1 , Toff=2.6  )
	model.addNativeProtein ( name="Bcl2" , EP=30000. , HLP=27. , EM=17. , HLM=9. , Ton=0.098 , Toff=3.73 )


### Signaling reactions : receptor activation and caspase cascade
def defineCaspaseCascadeReactions (model,ref_kb,ref_ku,ref_kc) :
	model.addReversibleReaction (name="TrailBindsReceptor",reactants=["TRAIL","R"],products=["TRAIL_R","TRAIL"],rates=[("kb",4.*ref_kb),("ku",1.e-3*ref_ku)])
	model.addReaction (name="ReceptorActivation",reactants=["TRAIL_R"], products=["R_act"],rate=("kc",0.01*ref_kc))
	model.addReaction (name="ActiveReceptorGivesTrailBack",reactants=["R_act"], products=["R"],rate=("k",10.*ref_ku))
	model.addCatalyticReaction (name="Caspase8Activation",substrate="pC8",catalyst="R_act",product="C8",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCatalyticReaction (name="Caspase3Activation",substrate="pC3",catalyst="C8",product="C3",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCatalyticReaction (name="Caspase6Activation",substrate="pC6",catalyst="C3",product="C6",rates=[("kb",0.),("ku",0.),("kc",0.)])
	model.addCatalyticReaction (name="Caspase8ActivationByC6",substrate="pC8",catalyst="C6",product="C8",rates=[("kb",0.),("ku",0.),("kc",0.)])
	model.addCatalyticReaction (name="PARPCleavage",substrate="PARP",catalyst="C3",product="cPARP",rates=[("kb",10.*ref_kb),("ku",ref_ku),("kc",20.*ref_kc)])

### Signaling reactions : inhibitors of caspase cascade
def defineCaspaseCascadeInhibitionReactions (model,ref_kb,ref_ku,ref_kc) :
	model.addReversibleReaction (name="flipBindsActivatedReceptor",reactants=["flip","R_act"],products=["flip_R_act"],rates=[("kb",10.*ref_kb),("ku",ref_ku)])
	model.addReversibleReaction (name="BarBindsCaspase8",reactants=["Bar","C8"],products=["Bar_C8"],rates=[("kb",10.*ref_kb),("ku",ref_ku)])
	model.addCatalyticReaction (name="XIAPdegradesCaspase3",substrate="C3",catalyst="XIAP",product="C3_ub",rates=[("kb",20.*ref_kb),("ku",ref_ku),("kc",0.1*ref_kc)])

### Signaling reactions : from Bid cleavage to mitochondrial pore formation
def defineUpstreamMompReactions (model,ref_kb,ref_ku,ref_kc,mito_cyto_vol_ratio) :
	model.addCatalyticReaction (name="BidCleavage",substrate="Bid",catalyst="C8",product="tBid",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCatalyticReaction (name="BaxActivation",substrate="Bax",catalyst="tBid",product="Bax_act",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addReversibleReaction (name="BaxTranslocation",reactants=["Bax_act"],products=["Bax_act_m"],rates=[("ktransf",3600.*0.01),("ktransb",3600.*1.)])
	model.addReversibleReaction (name="BaxDimerization",reactants=["Bax_act_m","Bax_act_m"],products=["Bax_act_m_2"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addReversibleReaction (name="BaxTetramerization",reactants=["Bax_act_m_2","Bax_act_m_2"],products=["Bax_act_m_4"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addReversibleReaction (name="BaxTetramerBindsPore",reactants=["Pore","Bax_act_m_4"],products=["Pore_Bax_act_m_4"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addReaction (name="PoreFormation",reactants=["Pore_Bax_act_m_4"],products=["Pore_act"],rate=("kc",ref_kc))


### Signaling reactions : inhibition of MOMP
def defineUpstreamMompInhibitionReactions (model,ref_kb,ref_ku,mito_cyto_vol_ratio) :
	model.addReversibleReaction (name="Mcl1BindstBid",reactants=["Mcl1","tBid"],products=["Mcl1_tBid"],rates=[("kb",10.*ref_kb),("ku",ref_ku)])
	model.addReversibleReaction (name="Bcl2BindsBax",reactants=["Bcl2","Bax_act_m"],products=["Bcl2_Bax_act_m"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addReversibleReaction (name="Bcl2BindsBax2",reactants=["Bcl2","Bax_act_m_2"],products=["Bcl2_Bax_act_m_2"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addReversibleReaction (name="Bcl2BindsBax4",reactants=["Bcl2","Bax_act_m_4"],products=["Bcl2_Bax_act_m_4"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addReversibleReaction (name="XIAPBindsApop",reactants=["XIAP","Apop"],products=["XIAP_Apop"],rates=[("kb",20.*ref_kb),("ku",ref_ku)])

### Signaling reactions : MOMP reactions
def defineMompReactions (model,ref_kb,ref_ku,ref_kc,mito_cyto_vol_ratio) :
	model.addCatalyticReaction (name="CyCRelease1",substrate="CyC_m",catalyst="Pore_act",product="CyC_r",rates=[("kb",20.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku),("kc",10.*ref_kc)])
	model.addCatalyticReaction (name="SmacRelease1",substrate="smac_m",catalyst="Pore_act",product="smac_r",rates=[("kb",20.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku),("kc",10.*ref_kc)])
	model.addReversibleReaction (name="CyCRelease2",reactants=["CyC_r"],products=["CyC"],rates=[("ktransf",3600.*1.),("ktransb",3600.*0.01)])
	model.addReversibleReaction (name="SmacRelease2",reactants=["smac_r"],products=["smac"],rates=[("ktransf",3600.*1.),("ktransb",3600.*0.01)])

### Signaling reactions : Downstream MOMP
def defineDownstreamMompReactions (model,ref_kb,ref_ku,ref_kc) :
	model.addReversibleReaction (name="SmacBindsXIAP",reactants=["smac","XIAP"],products=["smac_XIAP"],rates=[("kb",70.*ref_kb),("ku",ref_ku)])
	model.addCatalyticReaction (name="ApafActivation",substrate="Apaf",catalyst="CyC",product="Apaf_act",rates=[("kb",5.*ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addReversibleReaction (name="ApoptosomeFormation",reactants=["Apaf_act","pC9"],products=["Apop"],rates=[("kb",0.5*ref_kb),("ku",ref_ku)])
	model.addCatalyticReaction (name="ApopCleavesCaspase3",substrate="pC3",catalyst="Apop",product="C3",rates=[("kb",0.05*ref_kb),("ku",ref_ku),("kc",ref_kc)])



### set the degradation of non native protein species
def setDegradationRates (model) :
	model.setAllModifiedProteinDegradation (halfLife=5.)
	model.setModifiedProteinDegradation ( name="TRAIL" , halfLife=9.)
	model.setModifiedProteinDegradation ( name="Pore_act" , halfLife=1.9)
	model.setModifiedProteinDegradation ( name="flip_R_act" , halfLife=0.4)
	model.setModifiedProteinDegradation ( name="Mcl1_tBid" , halfLife=0.4)
	model.setModifiedProteinDegradation ( name="cPARP" , halfLife=27.)


### MAIN

## the model object
model = fpst.FluctuProtSTModel ("hEARM")

## define fluctuation models
defineNativeProteins (model=model)

## references rates parameters
ref_kb = 1.e-7*3600.
ref_ku = 0.001*3600.
ref_kc = 1.*3600.
mito_cyto_vol_ratio = 0.07

## caspase cascade
defineCaspaseCascadeReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc)
defineCaspaseCascadeInhibitionReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc)

## Upstream MOMP
defineUpstreamMompReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc,mito_cyto_vol_ratio=mito_cyto_vol_ratio)
defineUpstreamMompInhibitionReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,mito_cyto_vol_ratio=mito_cyto_vol_ratio)

## MOMP and downstream MOMP
defineMompReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc,mito_cyto_vol_ratio=mito_cyto_vol_ratio)
defineDownstreamMompReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc)

## degradation of modified forms
setDegradationRates (model=model)

## display reactions (not bad to check they are OK !)
for reac in model.signalingReactions :
	print "%s -> %s  __ rate = %e (%s)" % (str(reac.reactants),str(reac.products),reac.rate[1]/3600.,reac.name)

## generate cpp code to simulate that model
fpst.buildCppFromModel ( model=model , targetFolderPath = "hEARM_cpp_code_generated" )




