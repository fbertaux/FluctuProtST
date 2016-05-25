#!/usr/bin/python

################################################
# FluctuProtST, Version 1.2
# Francois Bertaux, Inria Paris-Rocquencourt
# francois.bertaux@inria.fr
# March 2015
################################################




#### imports
import FluctuProtST as fpst
from math import log

# the model object
model = fpst.FluctuProtSTModel ("SensingIP")

# parameters (time unit = hrs)
params = {}

# IP degrad, cell cycle, dilution
params["IP_degrad"] = log(2.) / 23.
params["doubling_time"] = 1.73
params["dilution_rate"] = log(2.) / params["doubling_time"]

# degradations
params["R_degrad"] = params["dilution_rate"] + 6.47 
params["ITF_degrad"] = params["dilution_rate"] + 0.157
params["uGFP_degrad"] = params["dilution_rate"]
params["GFP_degrad"] = params["dilution_rate"] + 0.117
params["IP_R_degrad"] = params["dilution_rate"]
params["IP_R_ITF_degrad"] = params["dilution_rate"]
params["ATF_degrad"] = params["dilution_rate"]

# kinetic reactions
params["folding_rate"] = 3.32
params["uGFP_degrad"] += params["folding_rate"] # weird choice on purpose, think about it
params["kb1"] = 3.63 * 60.
params["ku1"] = 0.54 * 60.
params["kb2"] = 0.012 * 60.
params["ku2"] = 1.88 * 60.
params["kc1"] = 0.716 * 60.
params["sigma"] = 0.128 * 60.

# gene expression parameters
params["mRNA_mean"] = 7.5 # Brown et al, 2013 ?
params["mRNA_HL"] = 2. * 0.17 # Brown et al, 2013 ?
params["mRNA_degrad"] = log(2.) / params["mRNA_HL"]
params["T_on"] = 1.25 
params["T_off"] = 1.5 

# misc
params["uGFP_mean_basal"] = params["GFP_degrad"] / params["folding_rate"]
print "uGFP basal mean level = %s" % params["uGFP_mean_basal"]

# define native (i.e. fluctuating) proteins
model.addNativeProtein ( name="R" , EP=1. , HLP=log(2.)/params["R_degrad"] , EM=params["mRNA_mean"] , HLM=params["mRNA_HL"] , Ton=params["T_on"] , Toff=params["T_off"] )
model.addNativeProtein ( name="ITF" , EP=1. , HLP=log(2.)/params["ITF_degrad"] , EM=params["mRNA_mean"] , HLM=params["mRNA_HL"] , Ton=params["T_on"] , Toff=params["T_off"] )
model.addNativeProtein ( name="uGFP" , EP=params["uGFP_mean_basal"] , HLP=log(2.)/params["uGFP_degrad"] , EM=params["mRNA_mean"] , HLM=params["mRNA_HL"] , Ton=params["T_on"] , Toff=params["T_off"] )


# define signaling reactions
model.addReaction ( name="IP_binding" , reactants=["IP","R"] , products=["IP_R","IP"] , rate = ("kb1",params["kb1"]) )
model.addReaction ( name="IP_unbinding" , reactants=["IP_R"] , products=["R"] , rate = ("ku1",params["ku1"]) )
model.addReaction ( name="TF_binding" , reactants=["IP_R","ITF"] , products=["IP_R_ITF"] , rate = ("kb2",params["kb2"]) )
model.addReaction ( name="TF_unbinding" , reactants=["IP_R_ITF"] , products=["IP_R","ITF"] , rate = ("ku2",params["ku2"]) )
model.addReaction ( name="TF_activation" , reactants=["IP_R_ITF"] , products=["IP_R","ATF"] , rate = ("kc1",params["kc1"]) )
model.addReaction ( name="TF_dep_synth" , reactants=["ATF"] , products=["uGFP","ATF"] , rate = ("sigma",params["sigma"]) )
model.addReaction ( name="folding" , reactants=["uGFP"] , products=["uGFP","GFP"] , rate = ("folding_rate",params["folding_rate"]) ) # removing of uGFP done via degradation rate

# set the degradation of non native protein species
model.setModifiedProteinDegradation ( name="IP" , halfLife=log(2.)/params["IP_degrad"] )
model.setModifiedProteinDegradation ( name="IP_R" , halfLife=log(2.)/params["IP_R_degrad"] )
model.setModifiedProteinDegradation ( name="IP_R_ITF" , halfLife=log(2.)/params["IP_R_ITF_degrad"] )
model.setModifiedProteinDegradation ( name="ATF" , halfLife=log(2.)/params["ATF_degrad"] )
model.setModifiedProteinDegradation ( name="GFP" , halfLife=log(2.)/params["GFP_degrad"] )


# generate cpp code to simulate that model
fpst.buildCppFromModel ( model=model , targetFolderPath = "SensingIP_cpp_code_generated" , redoMain = False )



