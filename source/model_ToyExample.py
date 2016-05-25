#!/usr/bin/python

################################################
# FluctuProtST, Version 1.2
# Francois Bertaux, Inria Paris-Rocquencourt
# francois.bertaux@inria.fr
# March 2015
################################################


#########################################################################################
# Toy example.
#
# It considers two native proteins, DeathReceptor and Caspase.
# They are equipped with standard fluctuations models.
#
# In addition, the following signaling reactions are defined :
# DeathLigand + DeathReceptor <-> DeathLigand_DeathReceptor (0.0001 [Ligand]-1.hrs-1, 0.001 hrs-1)
# DeathLigand_DeathReceptor -> DeathLigand + ActiveDeathReceptor (10. hrs-1)
# ActiveDeathReceptor + Caspase <-> ActiveDeathReceptor_Caspase (0.1 #-1.hrs-1, 0.001 hrs-1)
# ActiveDeathReceptor_Caspase -> ActiveDeathReceptor + CleavedCaspase (10. hrs-1)
#
# Thus, when DeathLigand is non-zero, DeathReceptor is activated and can cleave Caspase.
##########################################################################################


#### imports
import FluctuProtST as fpst

# the model object
model = fpst.FluctuProtSTModel ("ToyModel")

# define native (i.e. fluctuating) proteins
model.addNativeProteinStdFluct ( name="DeathReceptor" , EP=10000. , dilutionHalfLife=27.  )
model.addNativeProteinStdFluct ( name="Caspase" , EP=50000. , dilutionHalfLife=27.  )

# define signaling reactions.
# One can use addReaction for uni-directional reactions, addReversibleReaction for reversible reactions, and addCatalyticReaction for catalytic reactions.
model.addReversibleReaction ( name="DeathLigandBindsDeathReceptor" , reactants=["DeathLigand","DeathReceptor"] , products=["DeathLigand_DeathReceptor"] , rates = [ ("kb",0.0001) , ("ku",0.001) ] ) 
model.addReaction ( name="DeathReceptorActivation" , reactants=["DeathLigand_DeathReceptor"] , products=["ActiveDeathReceptor","DeathLigand"] , rate = ("kc",10.) )
model.addCatalyticReaction ( name="CaspaseCleavage" , substrate="Caspase" , catalyst="ActiveDeathReceptor" , product="CleavedCaspase" , rates = [ ("kb",0.1) , ("ku",0.001) , ("kc",10.) ] )

# set the degradation of non native protein species: dilution by default, rapid degradation for CleavedCaspase
model.setAllModifiedProteinDegradation ( halfLife=27. )
model.setModifiedProteinDegradation ( name="CleavedCaspase" , halfLife=0.5)

# generate cpp code to simulate that model
fpst.buildCppFromModel ( model=model , targetFolderPath = "ToyModel_cpp_code_generated" )
