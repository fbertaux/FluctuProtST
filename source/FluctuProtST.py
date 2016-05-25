#!/usr/bin/python

################################################
# FluctuProtST, Version 1.2
# Francois Bertaux, Inria Paris-Rocquencourt
# francois.bertaux@inria.fr
# March 2015
################################################


#### imports
import os.path
import math
import shutil


####  classes to describe a HybridOdeSge model
class NativeProtein (object) :
	def __init__ ( self , name , kon , koff , ksm , rm , ksp , rp ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		self.name = name
		self.kon = kon
		self.koff = koff
		self.ksm = ksm
		self.rm = rm
		self.ksp = ksp
		self.rp = rp

class ModifiedProtein (object) :
	def __init__ ( self , name , degRate ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		self.name = name
		self.degRate = degRate

class SignalingReaction (object) :
	def __init__ ( self , name , reactants , products , rate ) :
		self.name = name
		self.rate = rate
		self.reactants = reactants
		self.products = products

class FluctuProtSTModel (object) :
	def __init__ ( self , name = "Mymodel" ) :
		self.name = name
		self.nativeProteins = []
		self.modifiedProteins = []
		self.signalingReactions = []
	def addNativeProtein ( self , name , Ton , Toff , EM , HLM , EP , HLP ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		for s in ["*",":"] : 
			if s in name : raise Exception ("Symbols * or , not allowed in names.")
		if name in [ prot.name for prot in self.nativeProteins ] : raise Exception ("Name already taken by a nativeProt.")
		if name in [ prot.name for prot in self.modifiedProteins ] : raise Exception ("Name already taken by a modifiedProt.")
		koff = 1./Ton
		kon = 1./Toff
		EG = kon / (kon+koff)
		rm = math.log(2.) / HLM 
		ksm = EM * rm / EG
		rp = math.log(2.) / HLP
		ksp = EP * rp / EM
		nativeProt = NativeProtein ( name=name , kon=kon , koff=koff , ksm=ksm , rm=rm , ksp=ksp , rp=rp )
		self.nativeProteins.append (nativeProt)
	def addNativeProteinStdFluct ( self , name , EP , dilutionHalfLife ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		if name in [ prot.name for prot in self.nativeProteins ] : raise Exception ("Name already taken by a nativeProt.")
		if name in [ prot.name for prot in self.modifiedProteins ] : raise Exception ("Name already taken by a modifiedProt.")
		self.addNativeProtein ( name=name , Ton=0.1 , Toff=2.58 , EM=17. , HLM=9. , EP=EP , HLP=dilutionHalfLife)		
	def addModifiedProtein ( self , name , degRate=0. ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		if name == "" : raise Exception ("Empty name forbidden.")
		for s in ["*",":"] : 
			if s in name : raise Exception ("Symbols * or , not allowed in names.")
		if name in [ prot.name for prot in self.nativeProteins ] : raise Exception ("Name already taken by a nativeProt.")
		if name in [ prot.name for prot in self.modifiedProteins ] : raise Exception ("Name already taken by a modifiedProt.")
		mprot = ModifiedProtein ( name , degRate )
		self.modifiedProteins.append (mprot)
	def setModifiedProteinDegradation ( self , name , halfLife ) :
		mprot_matches = [ mprot for mprot in self.modifiedProteins if mprot.name == name ]
		if len(mprot_matches) == 0 : raise Exception ("Not existing modifiedProt.")
		mprot_matches[0].degRate = math.log(2.) / halfLife
	def setAllModifiedProteinDegradation ( self , halfLife ) :
		for mprot in self.modifiedProteins : self.setModifiedProteinDegradation ( mprot.name , halfLife )
	def addReaction ( self , name , reactants , products , rate ) :
		if len (reactants) == 0 : raise Exception ("I don't accept 0-th order reaction (for now at least.)")
		for reactant in reactants :
			if reactant not in [ prot.name for prot in self.nativeProteins ] + [ mprot.name for mprot in self.modifiedProteins ] :
				self.addModifiedProtein ( name=reactant )
		for product in products :
			if product not in [ prot.name for prot in self.nativeProteins ] + [ mprot.name for mprot in self.modifiedProteins ] :
				self.addModifiedProtein ( name=product )
		reaction = SignalingReaction ( name=name , reactants=reactants , products=products  , rate=rate)
		self.signalingReactions.append (reaction)
	def addReversibleReaction ( self , name , reactants , products , rates ) :
		self.addReaction ( name=name+"_forward" , reactants=reactants , products=products , rate=rates[0] )
		self.addReaction ( name=name+"_backward" , reactants=products , products=reactants , rate=rates[1] )
	def addCatalyticReaction ( self , name , substrate , catalyst , product , rates ) :
		self.addReaction ( name=name+"_binding" , reactants=[substrate,catalyst] , products=[substrate+"_"+catalyst] , rate=rates[0] )
		self.addReaction ( name=name+"_unbinding" , reactants=[substrate+"_"+catalyst] , products=[substrate,catalyst] , rate=rates[1] )
		self.addReaction ( name=name+"_catalysis" , reactants=[substrate+"_"+catalyst] , products=[catalyst,product] , rate=rates[2] )
	def giveProteinIndexFromName ( self , name ) :
		if name not in [ nprot.name for nprot in self.nativeProteins ] + [ mprot.name for mprot in self.modifiedProteins ] : raise Exception ("Protein does not exist.")
		np_matches = [ nprot for nprot in self.nativeProteins if nprot.name==name ]
		if len(np_matches) > 0 : idxProt = self.nativeProteins.index(np_matches[0])
		mp_matches = [ mprot for mprot in self.modifiedProteins if mprot.name==name ]
		if len(mp_matches) > 0 : idxProt = len(self.nativeProteins)+self.modifiedProteins.index(mp_matches[0])
		return idxProt

#### method to place specific code into template files
def parseTemplateFileAndReplace ( templateFile , replacementList ) :
	with open ( templateFile , 'r' ) as readFile :
		data = readFile.read () 
	for replacement in replacementList :
		if data.find ( replacement[0] ) == -1 : raise Exception ("place-holder not found in template file")
		data = data.replace ( replacement[0] , replacement[1] )
	return data



#### method to build CPP files from a model
def buildCppFromModel ( model , targetFolderPath , redoMain=True ) :

	## create dir if needed
	# if os.path.exists (targetFolderPath) : raise Exception ("Target folder already exists, I don't like that.")
	if not os.path.exists (targetFolderPath) : os.mkdir (targetFolderPath)
	if not os.path.exists (targetFolderPath+"/libs") : os.mkdir (targetFolderPath+"/libs")

	## copy what to copy from template
	templateFolder = "template_cpp_code"
	for f in os.listdir ( templateFolder + "/libs" ) : shutil.copy ( templateFolder + "/libs/" + f , targetFolderPath+"/libs" )
	if redoMain :
		shutil.copy ( templateFolder + "/template_main.cpp" , targetFolderPath + "/main.cpp" )
	shutil.copy ( templateFolder + "/template_CellState.cpp" , targetFolderPath + "/CellState.cpp" )
	shutil.copy ( templateFolder + "/template_ModelParameters.hpp" , targetFolderPath + "/ModelParameters.hpp" )
	shutil.copy ( templateFolder + "/template_HybridSimulator.hpp" , targetFolderPath + "/HybridSimulator.hpp" )
	shutil.copy ( templateFolder + "/template_HybridSimulator.cpp" , targetFolderPath + "/HybridSimulator.cpp" )
	shutil.copy ( templateFolder + "/template_HybridRhs.hpp" , targetFolderPath + "/HybridRhs.hpp" )

	## construct ModelParameters.cpp and write it
	toInsert = "\tmf_NumGenes = " + str(len(model.nativeProteins)) + " ;\n"
	toInsert = toInsert + "\tmf_NumAllProteinSpecies = " + str(len(model.nativeProteins)+len(model.modifiedProteins)) + " ;\n"
	toInsert = toInsert + "\tmf_NumModifiedProteins = " + str(len(model.modifiedProteins)) + " ;\n"
	toInsert = toInsert + "\tmf_NumReacs = " + str(len(model.signalingReactions)) + " ;\n"
	toInsert = toInsert + "\n\tmf_kons = VecDoub ( mf_NumGenes , 0. ) ;\n"
	toInsert = toInsert + "\tmf_koffs = VecDoub ( mf_NumGenes , 0. ) ;\n"
	toInsert = toInsert + "\tmf_ksms = VecDoub ( mf_NumGenes , 0. ) ;\n"
	toInsert = toInsert + "\tmf_rms = VecDoub ( mf_NumGenes , 0. ) ;\n"
	toInsert = toInsert + "\tmf_ksps = VecDoub ( mf_NumGenes , 0. ) ;\n"
	toInsert = toInsert + "\tmf_rps = VecDoub ( mf_NumGenes , 0. ) ;\n"
	toInsert = toInsert + "\tmf_kreacs = VecDoub ( mf_NumReacs , 0. ) ;\n"
	toInsert = toInsert + "\tmf_degrates = VecDoub ( mf_NumModifiedProteins , 0. ) ;\n\n"
	for nativeProt in model.nativeProteins :
		toInsert = toInsert + "\tmf_kons[" + str(model.nativeProteins.index(nativeProt)) + "] = " + str(nativeProt.kon) + " ;\n"
		toInsert = toInsert + "\tmf_koffs[" + str(model.nativeProteins.index(nativeProt)) + "] = " + str(nativeProt.koff) + " ;\n"
		toInsert = toInsert + "\tmf_ksms[" + str(model.nativeProteins.index(nativeProt)) + "] = " + str(nativeProt.ksm) + " ;\n"
		toInsert = toInsert + "\tmf_rms[" + str(model.nativeProteins.index(nativeProt)) + "] = " + str(nativeProt.rm) + " ;\n"
		toInsert = toInsert + "\tmf_ksps[" + str(model.nativeProteins.index(nativeProt)) + "] = " + str(nativeProt.ksp) + " ;\n"
		toInsert = toInsert + "\tmf_rps[" + str(model.nativeProteins.index(nativeProt)) + "] = " + str(nativeProt.rp) + " ;\n\n"
	for reac in model.signalingReactions :
		toInsert = toInsert + "\tmf_kreacs[" + str(model.signalingReactions.index(reac)) + "] = " + str(reac.rate[1]) + " ; " + "//" + reac.name + "\n"
	toInsert = toInsert + "\n"
	for mprot in model.modifiedProteins :
		toInsert = toInsert + "\tmf_degrates[" + str(model.modifiedProteins.index(mprot)) + "] = " + str(mprot.degRate) + " ;\n"

	replacementList = [ ( "placeholder_parameter_init" , toInsert ) ]
	toWrite = parseTemplateFileAndReplace ( templateFile="template_cpp_code/template_ModelParameters.cpp" , replacementList=replacementList )
	writeFile = open ( targetFolderPath + "/ModelParameters.cpp" , 'w' )
	writeFile.write ( toWrite )

	## construct MrnaSimulator.hpp and write it
	toInsert = ""
	commonStr = "(VecDoub &s) {return mf_ModelParameters->mf_"
	for nativeProt in model.nativeProteins :
		idxProt = model.nativeProteins.index(nativeProt)
		toInsert = toInsert + "\tDoub rate" + str(4*idxProt) + commonStr + "koffs[" + str(idxProt) + "]*s[" + str(3*idxProt) + "];}\n"
		toInsert = toInsert + "\tDoub rate" + str(4*idxProt+1) + commonStr + "kons[" + str(idxProt) + "]*s[" + str(3*idxProt+1) + "];}\n"
		toInsert = toInsert + "\tDoub rate" + str(4*idxProt+2) + commonStr + "ksms[" + str(idxProt) + "]*s[" + str(3*idxProt) + "];}\n"
		toInsert = toInsert + "\tDoub rate" + str(4*idxProt+3) + commonStr + "rms[" + str(idxProt) + "]*s[" + str(3*idxProt+2) + "];}\n"
	replacementList = [ ( "placeholder_rate_functions" , toInsert ) ]
	toWrite = parseTemplateFileAndReplace ( templateFile="template_cpp_code/template_MrnaSimulator.hpp" , replacementList=replacementList )
	writeFile = open ( targetFolderPath + "/MrnaSimulator.hpp" , 'w' )
	writeFile.write ( toWrite )

	## construct MrnaSimulator.cpp and write it
	toInsert = ""
	for nativeProt in model.nativeProteins :
		idxProt = model.nativeProteins.index(nativeProt)
		for r in range(0,4) :
			toInsert = toInsert + "\tdispatch[" + str(4*idxProt+r) + "] = &MrnaSimulator::rate" + str(4*idxProt+r) + ";\n"
	replacementList = [ ( "placeholder_dispatch" , toInsert ) ]
	toWrite = parseTemplateFileAndReplace ( templateFile="template_cpp_code/template_MrnaSimulator.cpp" , replacementList=replacementList )
	writeFile = open ( targetFolderPath + "/MrnaSimulator.cpp" , 'w' )
	writeFile.write ( toWrite )

	## construct RHS
	toInsert = ""
	for nativeProt in model.nativeProteins :
		idxProt = model.nativeProteins.index(nativeProt)
		toInsert = toInsert + "\tdydx[" + str(idxProt) + "] = mf_ModelParameters->mf_ksps[" + str(idxProt) +  "] * (*mrnatable)[" + str(idxProt) + "][tindex] "
		toInsert = toInsert + "- mf_ModelParameters->mf_rps[" + str(idxProt) + "] * y[" + str(idxProt) + "] ; //" + nativeProt.name + "\n"
	toInsert = toInsert + "\n"
	for mProt in model.modifiedProteins :
		idxmProt = model.modifiedProteins.index(mProt) + len(model.nativeProteins)
		toInsert = toInsert + "\tdydx[" + str(idxmProt) + "] = - mf_ModelParameters->mf_degrates[" + str(idxmProt-len(model.nativeProteins)) + "] * y[" + str(idxmProt) + "] ; //" + mProt.name + "\n"
	toInsert = toInsert + "\n"
	for reac in model.signalingReactions :
		idxReac = model.signalingReactions.index(reac)
		toInsert = toInsert + "\tmf_computedReactionRates[" + str(idxReac) + "] = mf_ModelParameters->mf_kreacs[" + str(idxReac) + "] *"
		for reactant in reac.reactants :
			# find reactant either in native or modified prots
			nprot_reac_match = [ np for np in model.nativeProteins if np.name == reactant ]
			mprot_reac_match = [ mp for mp in model.modifiedProteins if mp.name == reactant ]
			if len(nprot_reac_match) > 0 : idxReactant = model.nativeProteins.index (nprot_reac_match[0])
			if len(mprot_reac_match) > 0 : idxReactant = model.modifiedProteins.index (mprot_reac_match[0]) + len(model.nativeProteins)
			toInsert = toInsert + " y[" + str(idxReactant) + "] *"
		toInsert = toInsert[:-2]
		toInsert = toInsert + " ; \n"
	toInsert = toInsert + "\n"
	for nprot in model.nativeProteins :
		idxProt = model.nativeProteins.index(nprot)
		toInsert = toInsert + "\tdydx[" + str(idxProt) + "] += "
		for reac in model.signalingReactions :
			idxReac = model.signalingReactions.index(reac)
			reactant_matches = [ reactant for reactant in reac.reactants if reactant == nprot.name ]		
			if len (reactant_matches) > 0 :
				toInsert = toInsert + "- " + str(len(reactant_matches)) + ". * mf_computedReactionRates[" + str(idxReac) + "] "
			product_matches = [ product for product in reac.products if product == nprot.name ]
			if len (product_matches) > 0 :
				toInsert = toInsert + "+ " + str(len(product_matches)) + ". * mf_computedReactionRates[" + str(idxReac) + "] "
		if toInsert[-3:-1] == "+=" : toInsert = toInsert + "0."
		toInsert = toInsert + " ;\n"
	toInsert = toInsert + "\n"
	for mprot in model.modifiedProteins :
		idxProt = model.modifiedProteins.index(mprot) + len(model.nativeProteins)
		toInsert = toInsert + "\tdydx[" + str(idxProt) + "] += "
		for reac in model.signalingReactions :
			idxReac = model.signalingReactions.index(reac)
			reactant_matches = [ reactant for reactant in reac.reactants if reactant == mprot.name ]		
			if len (reactant_matches) > 0 :
				toInsert = toInsert + "- " + str(len(reactant_matches)) + ". * mf_computedReactionRates[" + str(idxReac) + "] "
			product_matches = [ product for product in reac.products if product == mprot.name ]
			if len (product_matches) > 0 :
				toInsert = toInsert + "+ " + str(len(product_matches)) + ". * mf_computedReactionRates[" + str(idxReac) + "] "
		if toInsert[-3:-1] == "+= " : toInsert = toInsert + "0."
		toInsert = toInsert + " ;\n"
	toInsertRhs = toInsert

	## write RHS in HybridRhs.cpp
	replacementList = [ ( "placeholder_hybrid_rhs" , toInsertRhs ) ]
	toWrite = parseTemplateFileAndReplace ( templateFile="template_cpp_code/template_HybridRhs.cpp" , replacementList=replacementList )
	writeFile = open ( targetFolderPath + "/HybridRhs.cpp" , 'w' )
	writeFile.write ( toWrite )

	## construct CellState.hpp and write it
	toInsert = ""
	for nativeProt in model.nativeProteins :
		idxProt = model.nativeProteins.index(nativeProt)
		toInsert = toInsert + "\tinline Doub get_" + nativeProt.name + "_MrnaLevel () { return mf_GeneMrnas[" + str(3*idxProt+2) + "] ; }\n"
		toInsert = toInsert + "\tinline Doub get_" + nativeProt.name + "_Level () { return mf_AllProts[" + str(idxProt) + "] ; }\n"
	toInsert = toInsert + "\n"
	for modProt in model.modifiedProteins :
		idxProt = model.modifiedProteins.index(modProt) + len(model.nativeProteins)
		toInsert = toInsert + "\tinline Doub get_" + modProt.name + "_Level () { return mf_AllProts[" + str(idxProt) + "] ; }\n"
		toInsert = toInsert + "\tinline void set_" + modProt.name + "_Level ( Doub value ) { mf_AllProts[" + str(idxProt) + "] = value ; }\n"
	replacementList = [ ( "placeholder_name_access" , toInsert ) ]
	toWrite = parseTemplateFileAndReplace ( templateFile="template_cpp_code/template_CellState.hpp" , replacementList=replacementList )
	writeFile = open ( targetFolderPath + "/CellState.hpp" , 'w' )
	writeFile.write ( toWrite )



