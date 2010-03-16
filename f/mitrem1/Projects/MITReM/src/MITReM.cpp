//---------------------------------------------------------------------------

#include "MITReM.h"

#include <math.h>
#include <vector>
#include <algorithm>

#include "DatFileReader.h"
#include "MathematicsPhysics.h"

#include "ElectrolyteModel_Ideal.h"
#include "ElectrolyteModel_DH.h"
#include "ElectrolyteModel_MSA.h"
#include "ElectrolyteModel_Exponential.h"
#include "ElecReaction_BV.h"
#include "ElecReaction_BVads.h"
#include "ElecReaction_Linear.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MITReM::MITReM(const std::string &name)
{
	/*
	Required formats:

	*******************************
	.homreactions:
	Versie 1.0
	[nHomReactions] = ...
		<label> = ...	<kf> = ...	<kb> = ...
	[NReagents0] = ...
		<label> = ...	<stoich> = ...
	[NProducts0] = ...
		<label> = ...	<stoich> = ...
	*******************************
	.elecreactions:
	Versie 1.0
	[nElecReactions] = ...
		<label> = ...	<Type> = ...	<nElectrons> = ...	<kOxi> = ...	<kRed> = ...	<aOxi> = ...	<aRed> = ...
	[NRedAgents0] = ...
		<label> = ...	<stoich> = ...	<order> = ...
	[NOxiAgents0] = ...
		<label> = ...	<stoich> = ...	<order> = ...
	*******************************
	.models:
	Versie 1.0
	[Conductivity] = ...
	[electrolyteModel] = ...
	[electrostatics] = ...
	*******************************
	.electrolytesolution:
	Versie 1.0
	[solventDielectricConst] = ...
	[solventDynamicViscosity] = ...
	[solutionDensity] = ...
	[solutionKinematicViscosity] = ...
	[temperature] = ...
	[nIons] = ...
		<label> = ...	<z> = ...	<D> = ...	<d> = ...	<M> = ...	<cInlet> = ...
	*******************************
	*/

	electrolyteSolutionFile = name + ".electrolytesolution";
	modelsFile = name + ".models";
	homReactionsFile = name + ".homreactions";
	elecReactionsFile = name + ".elecreactions";

	std::string choiceElectrolyteModel;
	std::string choiceElectrostaticsModel;

	DatFileReader datFile(electrolyteSolutionFile);
	unsigned nIons = datFile.readMultipleVector_nVectors("[nIons]");
	electrolyteSolution = new ElectrolyteSolution(nIons);

	// read models
	{
		DatFileReader datFile(modelsFile);

		datFile.readScalar("[conductivity]",conductivity);

		datFile.readScalar("[electrolyteModel]",choiceElectrolyteModel);
		datFile.readScalar("[electrostaticsModel]",choiceElectrostaticsModel);

		if (choiceElectrolyteModel == "Ideal")
		{
			electrolyteModel = new ElectrolyteModel_Ideal(electrolyteSolution);
		}
		else if (choiceElectrolyteModel == "DH")
		{
			electrolyteModel = new ElectrolyteModel_DH(electrolyteSolution);
		}
		else if (choiceElectrolyteModel == "MSA")
		{
			electrolyteModel = new ElectrolyteModel_MSA(electrolyteSolution);
		}
		else if (choiceElectrolyteModel == "Exponential")
		{
			electrolyteModel = new ElectrolyteModel_Exponential(electrolyteSolution);
		}
		else errorInvalidModel(modelsFile,choiceElectrolyteModel);

		if (choiceElectrostaticsModel == "Electroneutrality")
		{
			electrostatics = new Electrostatics_Electroneutrality(electrolyteSolution);
		}
		else if (choiceElectrostaticsModel == "Poisson")
		{
			electrostatics = new Electrostatics_Poisson(electrolyteSolution);
		}
		else errorInvalidModel(modelsFile,choiceElectrostaticsModel);
	}
	
	//read electrolyte solution
	{
		double doubleValue;
		std::string stringValue;
		int intValue;

		DatFileReader datFile(electrolyteSolutionFile);

		//nIons = datFile.readMultipleVector_nVectors("[nIons]");
		//electrolyteSolution = new ElectrolyteSolution(nIons);

		datFile.readScalar("[solutionKinematicViscosity]",doubleValue);
		electrolyteSolution->setSolutionKinematicViscosity(doubleValue);
		datFile.readScalar("[solutionTemperature]",doubleValue);
		electrolyteSolution->setSolutionTemperature(doubleValue);
		electrolyteSolution->setSolutionPotential(0.);
		datFile.readScalar("[solventDielectricConstant]",doubleValue);
		electrolyteSolution->setSolventDielectricConstant(doubleValue);
		for (unsigned i=0; i<nIons; i++)
		{
			datFile.readMultipleVector("[nIons]", i, "<label>", stringValue);
			electrolyteSolution->setIonLabel(i,stringValue);
			datFile.readMultipleVector("[nIons]", i, "<z>", intValue);
			electrolyteSolution->setIonChargeNumber(i,intValue);
			datFile.readMultipleVector("[nIons]", i, "<D>", doubleValue);
			electrolyteSolution->setIonDiffusionConstant(i,doubleValue);
			datFile.readMultipleVector("[nIons]", i, "<cInlet>", doubleValue);
			electrolyteSolution->setIonConcentration(i,doubleValue);
			electrolyteSolution->setIonInletConcentration(i,doubleValue);
		}

		if (electrolyteSolution->getIonChargeNumber(0) == 0)
		{
			std::cout << "INPUT FILE ERROR:\nThe first ion can not be a neutral species, "\
			"because this gives rise to linearly dependent equations." << std::endl;
			system("pause");
			exit(1);
		}

		if (choiceElectrolyteModel == "DH")
		{
			datFile.readScalar("[solutionDensity]",doubleValue);
			electrolyteSolution->setSolutionDensity(doubleValue);
			datFile.readScalar("[solventDynamicViscosity]",doubleValue);
			electrolyteSolution->setSolventDynamicViscosity(doubleValue);
			for (unsigned i=0; i<nIons; i++)
			{
				datFile.readMultipleVector("[nIons]", i, "<M>", doubleValue);
				electrolyteSolution->setIonMolarMass(i,doubleValue);
			}
		}
		else if (choiceElectrolyteModel == "MSA")
		{
			datFile.readScalar("[solutionDensity]",doubleValue);
			electrolyteSolution->setSolutionDensity(doubleValue);
			datFile.readScalar("[solventDynamicViscosity]",doubleValue);
			electrolyteSolution->setSolventDynamicViscosity(doubleValue);
			for (unsigned i=0; i<nIons; i++)
			{
				datFile.readMultipleVector("[nIons]", i, "<D>", doubleValue);
				electrolyteSolution->setIonDiffusionConstant(i,doubleValue);
				datFile.readMultipleVector("[nIons]", i, "<d>", doubleValue);
				electrolyteSolution->setIonDiameter(i,doubleValue);
				datFile.readMultipleVector("[nIons]", i, "<M>", doubleValue);
				electrolyteSolution->setIonMolarMass(i,doubleValue);
			}
		}
		else if (choiceElectrolyteModel == "Exponential")
		{
			for (unsigned i=0; i<nIons; i++)
			{
				datFile.readMultipleVector("[nIons]", i, "<DLim>", doubleValue);
				electrolyteSolution->setIonDiffusionConstantLimit(i,doubleValue);
				datFile.readMultipleVector("[nIons]", i, "<DPow>", doubleValue);
				electrolyteSolution->setIonDiffusionConstantPower(i,doubleValue);
			}
		}
	}

	//read homogeneous reactions
	{
		//unsigned nIons = electrolyteSolution->getNIons();

		DatFileReader datFile(homReactionsFile);

		nHomReactions = datFile.readMultipleVector_nVectors("[nHomReactions]");
		homReactions = new HomReaction*[nHomReactions];
		for (unsigned r=0; r<nHomReactions; r++) 
		{
			char buf[10];
			sprintf(buf,"%d",r);
			std::string rString(buf);
			std::string nReagentsr = "[nReagents" + rString + ']';
			std::string nProductsr = "[nProducts" + rString + ']';

			unsigned nReagents = datFile.readMultipleVector_nVectors(nReagentsr);
			if (nReagents == 0) errorZero(homReactionsFile,"nReagents");
			unsigned nProducts = datFile.readMultipleVector_nVectors(nProductsr);
			if (nProducts == 0) errorZero(homReactionsFile,"nProducts");
			homReactions[r] = new HomReaction(electrolyteSolution,electrolyteModel,nReagents,nProducts);

			std::string label;
			datFile.readMultipleVector("[nHomReactions]", r, "<label>", label);
			homReactions[r]->setLabel(label);
			double kf,kb;
			datFile.readMultipleVector("[nHomReactions]", r, "<kf>", kf);
			homReactions[r]->setForwardRateConstant(kf);
			if (kf == 0) errorZero(homReactionsFile,"kf");
			datFile.readMultipleVector("[nHomReactions]", r, "<kb>", kb);
			homReactions[r]->setBackwardRateConstant(kb);
			if (kb == 0) errorZero(homReactionsFile,"kb");
			homReactions[r]->setEquilibriumConstant(kf/kb);

			int stoich;
			for (unsigned i=0; i<nReagents; i++) 
			{
				homReactions[r]->setReagents(i,nIons);
				datFile.readMultipleVector(nReagentsr, i, "<label>", label);
				datFile.readMultipleVector(nReagentsr, i, "<stoich>", stoich);
				if (stoich != -1)
				{
					std::cout << "For the moment, only -1 is allowed for the stoichiometric coefficient of a reagent in a homogeneous reaction..." << std::endl;
					system("pause");
					exit(1);
				}
				homReactions[r]->setStoichReag(i,stoich);
				for (unsigned j=0; j<nIons; j++) 
				{
					if (electrolyteSolution->getIonLabel(j) == label) 
					{
						homReactions[r]->setReagents(i,j);
						break;
					}
				}
				if (homReactions[r]->getReagents(i) == nIons)
					errorConflictingData(homReactionsFile,electrolyteSolutionFile,"ION LABELS");
			}

			for (unsigned i=0; i<nProducts; i++) 
			{
				homReactions[r]->setProducts(i,nIons);
				datFile.readMultipleVector(nProductsr, i, "<label>", label);
				datFile.readMultipleVector(nProductsr, i, "<stoich>", stoich);
				if (stoich != 1)
				{
					std::cout << "For the moment, only +1 is allowed for the stoichiometric coefficient of a product in a homogeneous reaction..." << std::endl;
					system("pause");
					exit(1);
				}
				homReactions[r]->setStoichProd(i,stoich);
				for (unsigned j=0; j<nIons; j++) 
				{
					if (electrolyteSolution->getIonLabel(j) == label) 
					{
						homReactions[r]->setProducts(i,j);
						break;
					}
				}
				if (homReactions[r]->getProducts(i) == nIons)
					errorConflictingData(homReactionsFile,electrolyteSolutionFile,"ION LABELS");
			}

			// check charge conservation for homogenous reactions

			int sumChargeReagents = 0, sumChargeProducts = 0 ;
			
			for(unsigned i=0; i<nReagents; i++ )
			{
				sumChargeReagents += electrolyteSolution->getIonChargeNumber(homReactions[r]->getReagents(i));
			} 	

			for(unsigned j=0; j<nProducts; j++)
			{
				sumChargeProducts += electrolyteSolution->getIonChargeNumber(homReactions[r]->getProducts(j));
			}

			if (sumChargeReagents != sumChargeProducts)
			{
				std::cout << "NO charge conservation)" << std::endl;
				exit(1);
			}

		}
	}

	//read electrochemical reactions
	{
		//unsigned nIons = electrolyteSolution->getNIons();

		double doubleValue;
		std::string stringValue;
		int intValue;
		unsigned unsignedValue;

		DatFileReader datFile(elecReactionsFile);

		nElecReactions = datFile.readMultipleVector_nVectors("[nElecReactions]");
		elecReactions = new ElecReaction*[nElecReactions];
		unsigned rType0 = 0;
		unsigned rType1 = 0;
		for (unsigned r=0; r<nElecReactions; r++)
		{
			char buf[10];
			sprintf(buf,"%d",r);
			std::string rString(buf);
			std::string nAgentsRedr = "[nAgentsRed" + rString + ']';
			std::string nAgentsOxir = "[nAgentsOxi" + rString + ']';
			unsigned nAgentsRed, nAgentsOxi;

			nAgentsRed = datFile.readMultipleVector_nVectors(nAgentsRedr);
			nAgentsOxi = datFile.readMultipleVector_nVectors(nAgentsOxir);
			datFile.readMultipleVector("[nElecReactions]", r, "<type>", stringValue);
			if (stringValue == "Linear")
			{
				elecReactions[r] = new ElecReaction_Linear(electrolyteSolution,nAgentsRed,nAgentsOxi);
				datFile.readMultipleVector("[nElecReactions]", rType0, "<a>", doubleValue);
				elecReactions[r]->setKinParam(0,doubleValue);
				datFile.readMultipleVector("[nElecReactions]", rType0, "<b>", doubleValue);
				elecReactions[r]->setKinParam(1,doubleValue);
				rType0++;
			}
			else if (stringValue == "BV")
			{
				elecReactions[r] = new ElecReaction_BV(electrolyteSolution,nAgentsRed,nAgentsOxi);
				datFile.readMultipleVector("[nElecReactions]", rType1, "<kOxi>", doubleValue);
				elecReactions[r]->setKinParam(0,doubleValue);
				datFile.readMultipleVector("[nElecReactions]", rType1, "<kRed>", doubleValue);
				elecReactions[r]->setKinParam(1,doubleValue);
				datFile.readMultipleVector("[nElecReactions]", rType1, "<aOxi>", doubleValue);
				elecReactions[r]->setKinParam(2,doubleValue);
				datFile.readMultipleVector("[nElecReactions]", rType1, "<aRed>", doubleValue);
				elecReactions[r]->setKinParam(3,doubleValue);
				rType1++;
			}
			else
			{
				errorConflictingData(elecReactionsFile,"MITReM.cpp","ELECTRODE REACTION TYPES");
			}

			datFile.readMultipleVector("[nElecReactions]", r, "<label>", stringValue);
			elecReactions[r]->setLabel(stringValue);
			datFile.readMultipleVector("[nElecReactions]", r, "<nElectrons>", unsignedValue);
			elecReactions[r]->setNElectrons(unsignedValue);

			for (unsigned i=0; i<nAgentsRed; i++)
			{
				elecReactions[r]->setAgentsRed(i,nIons);
				std::string agentRediLabel;
				datFile.readMultipleVector(nAgentsRedr, i, "<label>", agentRediLabel);
				datFile.readMultipleVector(nAgentsRedr, i, "<stoich>", intValue);
				if (intValue >= 0)
				{
					std::cout << "Stoichiometric coefficients of reducing agents must be negative!" << std::endl;
					system("pause");
					exit(1);
				}
				elecReactions[r]->setStoichRed(i,intValue);
				datFile.readMultipleVector(nAgentsRedr, i, "<order>", doubleValue);
				if (doubleValue != 1)
				{
					std::cout << "Only first order electrode reactions are allowed for the moment." << std::endl;
					system("pause");
					exit(1);
				}
				elecReactions[r]->setOrderRed(i,doubleValue);
				for (unsigned j=0; j<nIons; j++)
				{
					if (electrolyteSolution->getIonLabel(j) == agentRediLabel)
					{
						elecReactions[r]->setAgentsRed(i,j);
						break;
					}
				}
				if (elecReactions[r]->getAgentsRed(i) == nIons)
				{
					errorConflictingData(elecReactionsFile,electrolyteSolutionFile,"ION LABELS");
				}
			}

			for (unsigned i=0; i<nAgentsOxi; i++)
			{
				elecReactions[r]->setAgentsOxi(i,nIons);
				std::string agentRediLabel;
				datFile.readMultipleVector(nAgentsOxir, i, "<label>", agentRediLabel);
				datFile.readMultipleVector(nAgentsOxir, i, "<stoich>", intValue);
				if (intValue <= 0)
				{
					std::cout << "Stoichiometric coefficients of oxidizing agents must be positive!" << std::endl;
					system("pause");
					exit(1);
				}
				elecReactions[r]->setStoichOxi(i,intValue);
				datFile.readMultipleVector(nAgentsOxir, i, "<order>", doubleValue);
				if (doubleValue != 1)
				{
					std::cout << "Only first order electrode reactions are allowed for the moment." << std::endl;
					system("pause");
					exit(1);
				}
				elecReactions[r]->setOrderOxi(i,doubleValue);
				for (unsigned j=0; j<nIons; j++)
				{
					if (electrolyteSolution->getIonLabel(j) == agentRediLabel)
					{
						elecReactions[r]->setAgentsOxi(i,j);
						break;
					}
				}
				if (elecReactions[r]->getAgentsOxi(i) == nIons)
				{
					errorConflictingData(elecReactionsFile,electrolyteSolutionFile,"ION LABELS");
				}
			}

			// check charge conservation for electrode reactions

			unsigned sumChargeRed = 0, sumChargeOxi = 0 /*, sumChargeElec = 0*/ ;
			
			for(unsigned i=0; i<nAgentsRed; i++ )
			{
				sumChargeRed += electrolyteSolution->getIonChargeNumber(elecReactions[r]->getAgentsRed(i));
			} 	

			for(unsigned j=0; j<nAgentsOxi; j++)
			{
				sumChargeOxi += electrolyteSolution->getIonChargeNumber(elecReactions[r]->getAgentsOxi(j));
			}

			sumChargeOxi -= elecReactions[r]->getNElectrons();	

			if (sumChargeRed != sumChargeOxi)
			{
				std::cout << "NO charge conservation" << std::endl;
				exit(1);
			}
		}
	}

	//read gas reactions
	{
		//unsigned nIons = electrolyteSolution->getNIons();

		DatFileReader datFile(elecReactionsFile);

		nGasReactions = datFile.readMultipleVector_nVectors("[nGasReactions]");
		gasReactions = new GasReaction*[nGasReactions];
		for (unsigned r=0; r<nGasReactions; r++)
		{
			gasReactions[r] = new GasReaction(electrolyteSolution);

			std::string label;
			datFile.readMultipleVector("[nGasReactions]", r, "<label>", label);
			gasReactions[r]->setLabel(label);

			double k;
			datFile.readMultipleVector("[nGasReactions]", r, "<k>", k);
			gasReactions[r]->setKinParam(0,k);

			double cSat;
			datFile.readMultipleVector("[nGasReactions]", r, "<cSat>", cSat);
			gasReactions[r]->setKinParam(1,cSat);

			gasReactions[r]->setDissolvedGas(nIons);
			std::string dissolvedGasLabel;
			datFile.readMultipleVector("[nGasReactions]", r, "<dissolvedGas>", dissolvedGasLabel);
			for (unsigned j=0; j<nIons; j++)
			{
				if (electrolyteSolution->getIonLabel(j) == dissolvedGasLabel)
				{
					gasReactions[r]->setDissolvedGas(j);
					break;
				}
			}
			if (gasReactions[r]->getDissolvedGas() == nIons)
			{
				errorConflictingData(elecReactionsFile,electrolyteSolutionFile,"ION LABELS");
			}
		}
	}

	//unsigned nIons = electrolyteSolution->getNIons();
	f = new double[nHomReactions];
	x = new double[nHomReactions];
	dfdx = new double*[nHomReactions];
	homReactionStoichMat = new int*[nHomReactions];
	for (unsigned r=0; r<nHomReactions; r++)
	{
		dfdx[r] = new double[nHomReactions];
		homReactionStoichMat[r] = new int[nIons];
	}
	cSave = new double[nIons];

	calcEquilibrium();
	setConductivityCorrectionFactor();
}
//---------------------------------------------------------------------------
MITReM::MITReM(const std::string& s_database, const std::string& s_ec_label)
{
  /*
  Required format:

  <ec label="." version=".">
   <electrolyte model=".">
    <solvent dielectricconstant="." dynamicviscosity="."/>
    <solution kinematicviscosity="." density="." temperature="."/>
    <species label="." z="." D="." cInlet="."/>
    <species label="." z="." D="." cInlet="."/>
   </electrolyte>
   <elecreactions>
    <reaction label="." model="." nElectrons="." E0="." kOxi="." kRed="." aOxi="." aRed=".">
     <agent label="." type="." stoich="." order="."/>
    </reaction>
   </elecreactions>
   <gasreactions>
    <reaction label="." k="." cSat="." dissolvedGas=".">
   </gasreactions>
   <homreactions>
    <reaction label="." kf="." kb=".">
     <reagent label="." stoich="."/>
     <product label="." stoich="."/>
    </reaction>
   </homreactions>
   <models>
    <conductivity value="."/>
    <electrostatics model="."/>
   </models>
  </ec>
  */


  // open database (xdb) and set <ec /> node
  XMLNode xdb = XMLNode::openFileHelper(s_database.c_str(),"ecdb");
  XMLNode xec = xdb.getChildNodeWithAttribute("ec","label",s_ec_label.c_str());
  if (xec.isEmpty()) {
    std::cerr << "database doesn't have <ecdb><ec label=\"" << s_ec_label << "\" /></ecdb>!" << std::endl;
    throw 400;
  }


  // setup electrolyte based on number of ions
  const unsigned nIons = (unsigned) xec.getChildNode("electrolyte").nChildNode("species");
  electrolyteSolution = new ElectrolyteSolution(nIons);


  // setup models
  {
    // <models/> node
    XMLNode m = xec.getChildNode("models");

    conductivity = m.getChildNode("conductivity").getAttribute< double >("value");

    const std::string choiceElectrolyteModel =
      xec.getChildNode("electrolyte").getAttribute< std::string >("model");
    if (choiceElectrolyteModel == "Ideal")
      electrolyteModel = new ElectrolyteModel_Ideal(electrolyteSolution);
    else if (choiceElectrolyteModel == "DH")
      electrolyteModel = new ElectrolyteModel_DH(electrolyteSolution);
    else if (choiceElectrolyteModel == "MSA")
      electrolyteModel = new ElectrolyteModel_MSA(electrolyteSolution);
    else if (choiceElectrolyteModel == "Exponential")
      electrolyteModel = new ElectrolyteModel_Exponential(electrolyteSolution);
    else
      errorInvalidModel(s_ec_label,choiceElectrolyteModel);

    const std::string choiceElectrostaticsModel =
      m.getChildNode("electrostatics").getAttribute< std::string >("model");
    if (choiceElectrostaticsModel == "Electroneutrality")
      electrostatics = new Electrostatics_Electroneutrality(electrolyteSolution);
    else if (choiceElectrostaticsModel == "Poisson")
      electrostatics = new Electrostatics_Poisson(electrolyteSolution);
    else
      errorInvalidModel(s_ec_label,choiceElectrostaticsModel);
  }


  // setup electrolyte solution
  {
    XMLNode s;  // helper node

    // setup solution and solvent
    s = xec.getChildNode("electrolyte").getChildNode("solution");
    electrolyteSolution->setSolutionKinematicViscosity(s.getAttribute< double >("kinematicviscosity"));
    electrolyteSolution->setSolutionTemperature(s.getAttribute< double >("temperature"));
    electrolyteSolution->setSolutionPotential(0.);
    s = xec.getChildNode("electrolyte").getChildNode("solvent");
    electrolyteSolution->setSolventDielectricConstant(s.getAttribute< double >("dielectricconstant"));
  
    // setup present species
    for (unsigned i=0; i<nIons; ++i) {
      s = xec.getChildNode("electrolyte").getChildNode("species",(int) i);
      electrolyteSolution->setIonLabel(i,s.getAttribute< std::string >("label"));
      electrolyteSolution->setIonChargeNumber(i,s.getAttribute< int >("z"));
      electrolyteSolution->setIonDiffusionConstant(i,s.getAttribute< double >("D"));
      electrolyteSolution->setIonInletConcentration(i,s.getAttribute< double >("cInlet"));
      electrolyteSolution->setIonConcentration(i,s.getAttribute< double >("cInlet"));

      electrolyteSolution->setIonTVExpansionCoefficient(i,s.getAttribute< double >("alpha",0.));
      electrolyteSolution->setIonCDensificationCoefficient(i,s.getAttribute< double >("beta",0.));
      electrolyteSolution->setIonMMagneticSusceptibility(i,s.getAttribute< double >("MMChi",0.));
    }
  
    // setup properties specific to model
    const std::string choiceElectrolyteModel = xec.getChildNode("electrolyte").getAttribute< std::string >("model");
    if (choiceElectrolyteModel=="DH" || choiceElectrolyteModel=="MSA") {

      electrolyteSolution->setSolutionDensity(xec.getChildNode("electrolyte").getChildNode("solution").getAttribute< double >("density"));
      electrolyteSolution->setSolventDynamicViscosity(xec.getChildNode("electrolyte").getChildNode("solvent").getAttribute< double >("dynamicviscosity"));
      for (unsigned i=0; i<nIons; ++i) {
        s = xec.getChildNode("electrolyte").getChildNode("species",(int) i);
        electrolyteSolution->setIonMolarMass(i,s.getAttribute< double >("M"));
        if (choiceElectrolyteModel=="MSA") {
          electrolyteSolution->setIonDiffusionConstant(i,s.getAttribute< double >("D"));
          electrolyteSolution->setIonDiameter(i,s.getAttribute< double >("diam"));
        }
      }

    }
    else if (choiceElectrolyteModel=="Exponential") {

      for (unsigned i=0; i<nIons; ++i) {
        s = xec.getChildNode("electrolyte").getChildNode("species",(int) i);
        electrolyteSolution->setIonDiffusionConstantLimit(i,s.getAttribute< double >("DLim"));
        electrolyteSolution->setIonDiffusionConstantPower(i,s.getAttribute< double >("DPow"));
      }

    }

    // display species properties
    std::cout << "ElectrolyteSolution species:" << std::endl;
    for (unsigned i=0; i<nIons; ++i)
      std::cout << "  species "  << i
                << "  label: \"" << getIonLabel(i) << "\""
                << "  z: "       << getIonChargeNumber(i)
                << "  D: "       << getIonDiffusionConstant(i)
                << "  cInlet: "  << getIonInletConcentration(i)
                << "  c: "       << getIonConcentration(i)
                << "  alpha: "   << getIonTVExpansionCoefficient(i)
                << "  beta: "    << getIonCDensificationCoefficient(i)
                << "  MMChi: "   << getIonMMagneticSusceptibility(i)
                << std::endl;

    // check for correctness
    if (electrolyteSolution->getIonChargeNumber(0)==0) {
      std::cout << "The first ion can not be a neutral species, because this gives rise to linearly dependent equations." << std::endl;
      system("pause");
      exit(1);
    }
  }


  // setup homogeneous reactions
  {
    nHomReactions = (unsigned) xec.getChildNode("homreactions").nChildNode("reaction");
    homReactions = new HomReaction*[nHomReactions];

    for (unsigned r=0; r<nHomReactions; ++r) {
      XMLNode xr = xec.getChildNode("homreactions").getChildNode("reaction",(int) r);
      unsigned nReagents = (unsigned) xr.nChildNode("reagent");
      unsigned nProducts = (unsigned) xr.nChildNode("product");
      homReactions[r] = new HomReaction(electrolyteSolution,electrolyteModel,nReagents,nProducts);

      double kf = xr.getAttribute< double >("kf");
      double kb = xr.getAttribute< double >("kb");
      if (kf==0)
        errorZero(s_ec_label,"kf");
      if (kb==0)
        errorZero(s_ec_label,"kb");
      homReactions[r]->setForwardRateConstant(kf);
      homReactions[r]->setBackwardRateConstant(kb);
      homReactions[r]->setEquilibriumConstant(kf/kb);
      homReactions[r]->setLabel(xr.getAttribute< std::string >("label"));

      for (unsigned a=0; a<nReagents; ++a) {
        XMLNode xa = xr.getChildNode("reagent",(int) a);
        homReactions[r]->setReagents(a,nIons);
        std::string stringValue = xa.getAttribute< std::string >("label");
				if (xa.getAttribute< int >("stoich") != -1)
				{
					std::cout << "For the moment, only -1 is allowed for the stoichiometric coefficient of a reagent in a homogeneous reaction..." << std::endl;
					system("pause");
					exit(1);
				}
        homReactions[r]->setStoichReag(a,xa.getAttribute< int >("stoich"));
        for (unsigned j=0; j<nIons; ++j)
          if (electrolyteSolution->getIonLabel(j) == stringValue) {
            homReactions[r]->setReagents(a,j);
            break;
          }
        if (homReactions[r]->getReagents(a) == nIons)
          errorConflictingData(s_ec_label,s_ec_label,"ION LABELS");
      }

      for (unsigned a=0; a<nProducts; ++a) {
        XMLNode xa = xr.getChildNode("product",(int) a);
        homReactions[r]->setProducts(a,nIons);
        std::string stringValue = xa.getAttribute< std::string >("label");
				if (xa.getAttribute< int >("stoich") != 1)
				{
					std::cout << "For the moment, only +1 is allowed for the stoichiometric coefficient of a product in a homogeneous reaction..." << std::endl;
					system("pause");
					exit(1);
				}
        homReactions[r]->setStoichProd(a,xa.getAttribute< int >("stoich"));
        for (unsigned j=0; j<nIons; ++j)
          if (electrolyteSolution->getIonLabel(j) == stringValue) {
            homReactions[r]->setProducts(a,j);
            break;
          }
        if (homReactions[r]->getProducts(a) == nIons)
          errorConflictingData(s_ec_label,s_ec_label,"ION LABELS");
      }

    }
  }


  // setup electrode reactions
  {
    nElecReactions = (unsigned) xec.getChildNode("elecreactions").nChildNode("reaction");
    elecReactions = new ElecReaction*[nElecReactions];

    for (unsigned r=0; r<nElecReactions; ++r) {
      XMLNode xr = xec.getChildNode("elecreactions").getChildNode("reaction",(int) r);

      unsigned nAgentsRed = 0;
      unsigned nAgentsOxi = 0;
      for (int a=0; a<xr.nChildNode("agent"); ++a) {
        std::string atype(xr.getChildNode("agent",a).getAttribute("type"));
        nAgentsRed += (atype=="red"? 1:0);
        nAgentsOxi += (atype=="oxi"? 1:0);
      }

      std::string stringValue = xr.getAttribute< std::string >("model");
      if (stringValue=="Linear")
      {
        elecReactions[r] = new ElecReaction_Linear(electrolyteSolution,nAgentsRed,nAgentsOxi);
        elecReactions[r]->setKinParam(0,xr.getAttribute< double >("a"));
        elecReactions[r]->setKinParam(1,xr.getAttribute< double >("b"));
      }
      else if (stringValue=="BV")
      {
        elecReactions[r] = new ElecReaction_BV(electrolyteSolution,nAgentsRed,nAgentsOxi);
        elecReactions[r]->setKinParam(0,xr.getAttribute< double >("kOxi"));
        elecReactions[r]->setKinParam(1,xr.getAttribute< double >("kRed"));
        elecReactions[r]->setKinParam(2,xr.getAttribute< double >("aOxi"));
        elecReactions[r]->setKinParam(3,xr.getAttribute< double >("aRed"));
      }
			else if (stringValue=="BVads")
      {
        elecReactions[r] = new ElecReaction_BVads(electrolyteSolution,nAgentsRed,nAgentsOxi);
        elecReactions[r]->setKinParam(0,xr.getAttribute< double >("kOxi"));
        elecReactions[r]->setKinParam(1,xr.getAttribute< double >("kRed"));
        elecReactions[r]->setKinParam(2,xr.getAttribute< double >("aOxi"));
        elecReactions[r]->setKinParam(3,xr.getAttribute< double >("aRed"));
        elecReactions[r]->setKinParam(4,xr.getAttribute< double >("KAds"));
        elecReactions[r]->setKinParam(5,xr.getAttribute< double >("aAds"));
      }
      else
      {
        errorConflictingData(s_ec_label,s_ec_label,"ELECTRODE REACTION TYPES");
      }

      elecReactions[r]->setLabel(xr.getAttribute< std::string >("label"));
      elecReactions[r]->setNElectrons(xr.getAttribute< unsigned >("nElectrons"));

      unsigned ired = 0;
      unsigned ioxi = 0;
      for (int a=0; a<xr.nChildNode("agent"); ++a) {
        XMLNode xa = xr.getChildNode("agent",a);
        std::string type  = xa.getAttribute< std::string >("type");
        std::string label = xa.getAttribute< std::string >("label");
        unsigned ulabel = nIons;
        for (unsigned j=0; j<nIons; ++j)
          if (electrolyteSolution->getIonLabel(j) == label) {
            ulabel = j;
            break;
          }
        if (ulabel==nIons)
          errorConflictingData(s_ec_label,s_ec_label,"ION LABELS");
        if (type=="red") {
					if (xa.getAttribute< int >("stoich") >= 0)
					{
						std::cout << "Stoichiometric coefficients of reducing agents must be negative!" << std::endl;
						system("pause");
						exit(1);
					}
          elecReactions[r]->setStoichRed(ired,xa.getAttribute< int >("stoich"));
					if (xa.getAttribute< double >("order") != 1)
					{
						std::cout << "Only first order electrode reactions are allowed for the moment." << std::endl;
						system("pause");
						exit(1);
					}
          elecReactions[r]->setOrderRed(ired,xa.getAttribute< double >("order"));
          elecReactions[r]->setAgentsRed(ired,ulabel);
          ++ired;
        }
        else if (type=="oxi") {
					if (xa.getAttribute< int >("stoich") <= 0)
					{
						std::cout << "Stoichiometric coefficients of oxidizing agents must be positive!" << std::endl;
						system("pause");
						exit(1);
					}
          elecReactions[r]->setStoichOxi(ioxi,xa.getAttribute< int >("stoich"));
					if (xa.getAttribute< double >("order") != 1)
					{
						std::cout << "Only first order electrode reactions are allowed for the moment." << std::endl;
						system("pause");
						exit(1);
					}
          elecReactions[r]->setOrderOxi(ioxi,xa.getAttribute< double >("order"));
          elecReactions[r]->setAgentsOxi(ioxi,ulabel);
          ++ioxi;
        }

      }

    }
  }


  // setup gas reactions
  {
    nGasReactions = (unsigned) xec.getChildNode("gasreactions").nChildNode("reaction");
    gasReactions = new GasReaction*[nGasReactions];

    for (unsigned r=0; r<nGasReactions; ++r) {
      XMLNode xr = xec.getChildNode("gasreactions").getChildNode("reaction",(int) r);
      gasReactions[r] = new GasReaction(electrolyteSolution);

      gasReactions[r]->setLabel(xr.getAttribute< std::string >("label"));
      gasReactions[r]->setKinParam(0,xr.getAttribute< double >("k"));
      gasReactions[r]->setKinParam(1,xr.getAttribute< double >("cSat"));

			std::string dissolvedGasLabel = xr.getAttribute< std::string >("dissolvedGas");
      gasReactions[r]->setDissolvedGas(nIons);
      for (unsigned j=0; j<nIons; j++)
        if (electrolyteSolution->getIonLabel(j)==dissolvedGasLabel) {
          gasReactions[r]->setDissolvedGas(j);
          break;
        }
      if (gasReactions[r]->getDissolvedGas()==nIons)
        errorConflictingData(s_ec_label,s_ec_label,"ION LABELS");
    }
  }


  // allocate homogeneous reactions things...
  f = new double[nHomReactions];
  x = new double[nHomReactions];
  dfdx = new double*[nHomReactions];
  homReactionStoichMat = new int*[nHomReactions];
  for (unsigned r=0; r<nHomReactions; r++) {
    dfdx[r] = new double[nHomReactions];
    homReactionStoichMat[r] = new int[nIons];
  }
  cSave = new double[nIons];


  // set conductivity correction factor
	calcEquilibrium();
  setConductivityCorrectionFactor();
}
//---------------------------------------------------------------------------
MITReM::~MITReM()
{
	// Delete in reverse order of construction!!!
	for (unsigned r=0; r<nElecReactions; r++)
	{
		delete elecReactions[r];
	}
	delete[] elecReactions;
	for (unsigned r=0; r<nHomReactions; r++) 
	{
		delete homReactions[r];
	}
	delete[] homReactions;
	delete electrostatics;
	delete electrolyteModel;
	delete electrolyteSolution;

	delete[] f;
	delete[] x;
	for (unsigned r=0; r<nHomReactions; r++)
	{
		delete[] dfdx[r];
		delete[] homReactionStoichMat[r];
	}
	delete[] dfdx;
	delete[] homReactionStoichMat;
	delete[] cSave;
}
//---------------------------------------------------------------------------


//--- METHODS ---------------------------------------------------------------
void MITReM::init(const double* c, double U, double T, double density)
{
	unsigned nIons = electrolyteSolution->getNIons();

	electrolyteSolution->setSolutionTemperature(T);
	electrolyteSolution->setSolutionDensity(density);
	for (unsigned i=0; i<nIons; i++)
	{
		electrolyteSolution->setIonConcentration(i,c[i]);
	}
	electrolyteSolution->setSolutionPotential(U);
	electrolyteModel->init(false);
}
//---------------------------------------------------------------------------
void MITReM::setConductivityCorrectionFactor()
{
	double conductivityTheoretical = electrolyteModel->calcConductivity();
	if (conductivity == 0)
	{
		std::cout << "\nNo measured conductivity was specified.\nNo corrections will be made to the transport properties."
			<< "\nTheoretical conductivity = " << conductivityTheoretical << " S/m" << std::endl;
	}
	else
	{
		electrolyteModel->init(false);
    double conductivityCorrectionFactor = electrolyteModel->calcConductivityCorrectionFactor(conductivity);
		std::cout << "\nTheoretical conductivity = " << conductivityTheoretical << " S/m" 
			<< "\nExperimental conductivity = " << conductivity << " S/m" 
			<< "\nAll diffusion coefficients will be multiplied with " << conductivityCorrectionFactor << std::endl;
	}
	//system("pause");
}
//---------------------------------------------------------------------------
void MITReM::assembleSystem() const
{
	unsigned nIons = electrolyteSolution->getNIons();

	// Each reaction r progresses by an amount of x[r].
	// This will cause a change in c[i] by an amount of s[r][i]*x[r].
	for (unsigned i=0; i<nIons; i++)
	{
		double dc = 0.;
		for (unsigned r=0; r<nHomReactions; r++)
		{
			dc += homReactionStoichMat[r][i]*x[r]; //quite a lot of zero's will be added...
		}
		electrolyteSolution->setIonConcentration(i, cSave[i]+dc);
	}

	// The functions that you want to minimize are the deviations from equilibrium f[r].
	// The unknowns are the progresses x[s] of the reactions.
	// For the Newton iterations you need df[r]/dx[s].
	electrolyteModel->init(false);
	for (unsigned r=0; r<nHomReactions; r++)
	{
		f[r] = homReactions[r]->calcDeviationFromEquilibrium();
		for (unsigned s=0; s<nHomReactions; s++)
		{
			dfdx[r][s] = homReactions[r]->calcDeviationFromEquilibriumDerivative(homReactionStoichMat[s]);
		}
	}
}
//---------------------------------------------------------------------------
double MITReM::calcStepLength() const
{
	//stepLength is also called the Newton damping factor

	// This function calculates the damping factor such that you will not get negative concentrations.
	
	unsigned nIons = electrolyteSolution->getNIons();

	double stepLength = 1.;
	for (unsigned i=0; i<nIons; i++)
	{
		double dcOld = 0.;
		double dcIncrement = 0.;
		for (unsigned r=0; r<nHomReactions; r++)
		{
			// -f[r] is dx[r]
			dcOld += homReactionStoichMat[r][i]*x[r]; //quite a lot of zero's will be added...
			dcIncrement -= homReactionStoichMat[r][i]*f[r]; //quite a lot of zero's will be added...
		}
		double dcNew = dcOld+dcIncrement;
		if (cSave[i]+dcNew <= 0.)
		{
			double stepLengthTemp = -(cSave[i]+dcOld)/dcIncrement * 0.9; //factor 0.9 to avoid c = 0
			if (stepLengthTemp < stepLength) stepLength = stepLengthTemp;
		}
	}
	return stepLength;
}
//---------------------------------------------------------------------------
void MITReM::calcEquilibrium()
{
	unsigned nIons = electrolyteSolution->getNIons();
	std::vector<double> c(nIons);
	for (unsigned i=0; i<nIons; i++)
	{
		c[i] = electrolyteSolution->getIonConcentration(i);
	}
	calcEquilibrium(c);
}
//---------------------------------------------------------------------------
void MITReM::calcEquilibrium(std::vector<double>& c)
{
	unsigned nIons = electrolyteSolution->getNIons();

	checkElectroneutrality();
	if (nHomReactions > 0)
	{
		std::cout << "Making a list of all the ions that participate in homogeneous reactions...\n";

		std::vector<unsigned> participatingIons;
		for (unsigned r=0; r<nHomReactions; r++)
		{
			unsigned nReagents = homReactions[r]->getNReagents();
			for (unsigned i = 0; i < nReagents; i++)
			{
				unsigned j = homReactions[r]->getReagents(i);
				std::vector<unsigned>::iterator result = find(participatingIons.begin(), participatingIons.end(), j);    
				if(result == participatingIons.end()) 
				{
					participatingIons.push_back(j);
				}                  
			}
			unsigned nProducts = homReactions[r]->getNProducts();
			for (unsigned i = 0; i < nProducts; i++)
			{
				unsigned j = homReactions[r]->getProducts(i);
				std::vector<unsigned>::iterator result = find(participatingIons.begin(), participatingIons.end(), j);    
				if(result == participatingIons.end()) 
				{
					participatingIons.push_back(j);
				}                  
			}
		}
		
		for (unsigned r=0; r<nHomReactions; r++)
		{
			x[r] = 0.;
			for (unsigned i=0; i<nIons; i++)
			{
				homReactionStoichMat[r][i] = homReactions[r]->getStoichOf(i);
			}
		}

		for (unsigned i=0; i<nIons; i++)
		{
			electrolyteSolution->setIonConcentration(i, c[i]);
			cSave[i] = c[i];
		}

		unsigned nZeroConcentrations = 0;
		for (unsigned i=0; i<participatingIons.size(); i++)
		{
			if (cSave[participatingIons[i]] == 0) 
				nZeroConcentrations++;
		}

		std::cout << "Attempting to achieve non-zero concentrations...\n";

		while (nZeroConcentrations > 0)
		{
			for (unsigned r=0; r<nHomReactions; r++)
			{
				electrolyteModel->init(false);
				x[r] += homReactions[r]->progressToEquilibrium();
			}
			unsigned nZeroConcentrationsTemp = 0;
			for (unsigned i=0; i<participatingIons.size(); i++)
			{
				if (electrolyteSolution->getIonConcentration(participatingIons[i]) == 0) 
					nZeroConcentrationsTemp++;
			}
			if (nZeroConcentrationsTemp == nZeroConcentrations)
			{
				std::cout << "INPUT FILE ERROR:\nCould not achieve non-zero concentrations.\n"\
										 "Please change some zero concentrations in the input file "\
										 "so that equilibrium is attainable at the inlet." << std::endl;
				system("pause");
				exit(1);
			}
			else if (nZeroConcentrationsTemp < nZeroConcentrations)
			{
				nZeroConcentrations = nZeroConcentrationsTemp;
			}
			else
			{
				std::cout << "IMPOSSIBLE ERROR:\nAn impossible error occurred while "\
					"making non-zero concentrations. Please contact the programmer." << std::endl;
				system("pause");
				exit(1);
			}
		}

		std::cout << "Computing equilibrium concentrations...\n";
		
		const unsigned iterMax = 1000;
		const double deviationMax = 1e-9;

		// Perform Newton iterations
		unsigned iter = 0;
		double deviation;
		do 
		{
			assembleSystem();
			solveGauss(dfdx,f,nHomReactions);
			double stepLength = calcStepLength();
			deviation = 0.;
			for (unsigned r=0; r<nHomReactions; r++)
			{
				x[r] -= stepLength*f[r];
				double deviationTemp = fabs(homReactions[r]->calcRelativeDeviationFromEquilibrium());
				if (deviationTemp > deviation) deviation = deviationTemp;
			}
			iter++;
		} while ((iter < iterMax) && (deviation > deviationMax));

		// Use the calculated x's to set the equilibrium concentrations and initialize electrolyteModel (with those concentrations)
		for (unsigned i=0; i<nIons; i++)
		{
			double dc = 0.;
			for (unsigned r=0; r<nHomReactions; r++)
			{
				dc += homReactionStoichMat[r][i]*x[r]; //quite a lot of zero's will be added...
			}
			electrolyteSolution->setIonConcentration(i,cSave[i]+dc);
			electrolyteSolution->setIonInletConcentration(i,cSave[i]+dc);
			c[i] = cSave[i]+dc;
		}
		electrolyteModel->init(false);

		// Write result to screen
		std::cout.precision(12);
		std::cout << "\nNumber of iterations needed: " << iter << "\n\nConcentrations:\n";
		for (unsigned i=0; i<nIons; i++)
		{
			std::cout << electrolyteSolution->getIonLabel(i) << '\t' << electrolyteSolution->getIonConcentration(i) << '\n';
		}
		std::cout << "\nDeviations from equilibrium, ln(activity quotient / K):\n";
		for (unsigned r=0; r<nHomReactions; r++)
		{
			std::cout << homReactions[r]->getLabel() << '\t' << homReactions[r]->calcRelativeDeviationFromEquilibrium() << '\n';
		}
		std::cout << std::endl;

		if (deviation > deviationMax)
		{
			std::cout << "\nERROR: Equilibrium at the inlet was not reached." << std::endl;
			system("pause");
			exit(1);
		}

		checkElectroneutrality();
	}
	electrolyteModel->init(false);
	//system("pause");
}
//---------------------------------------------------------------------------
void MITReM::checkElectroneutrality() const
{
	unsigned nIons = electrolyteSolution->getNIons();

	std::cout << "Checking electroneutrality...\n";
	double charge = 0.;
	double chargePos = 0.;
	double chargeNeg = 0.;
	for (unsigned i=0; i<nIons; i++)
	{
		double zici = electrolyteSolution->getIonChargeNumber(i)*electrolyteSolution->getIonConcentration(i);
		charge += zici;
		if (zici > 0.) chargePos += zici;
		else if (zici < 0.) chargeNeg += zici;
	}
	//double f = sqrt(-chargeNeg/chargePos);

	const double maxCharge = 1e-12;
	if (fabs(charge) > maxCharge)
	{
		std::cout << "\nINPUT FILE ERROR: Electroneutrality is violated.\nSum(z_i*c_i) = " << charge << std::endl;
		system("pause");
		exit(1);
		//std::cout << "\nThe adjusted concentrations are\n";
		//for (unsigned i=0; i<nIons; i++)
		//{
		//	int zi = electrolyteSolution->getIonChargeNumber(i);
		//	double ci = electrolyteSolution->getIonConcentration(i);
		//	if (zi > 0.) electrolyteSolution->setIonConcentration(i,ci*f);
		//	else if (zi < 0.) electrolyteSolution->setIonConcentration(i,ci/f);
		//	std::cout << electrolyteSolution->getIonLabel(i) << '\t' << electrolyteSolution->getIonConcentration(i) << std::endl;
		//}
		//std::cout << std::endl;
	}
}
//---------------------------------------------------------------------------
void MITReM::setElectrolyteModel(std::string EM)
{
	delete electrolyteModel;
	if (EM == "MSA") electrolyteModel = new ElectrolyteModel_MSA(electrolyteSolution);
	else if (EM == "Ideal") electrolyteModel = new ElectrolyteModel_Ideal(electrolyteSolution);
	else exit(1);
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void MITReM::errorInvalidModel(const std::string &file, const std::string &model) const
{
	std::cout << "ERROR IN " << file << ".\
			\nTHE MODEL " << model <<" IS NOT ALLOWED." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------
void MITReM::errorConflictingData(const std::string &file1, const std::string &file2, const std::string &conflict) const
{
	std::cout << "CONFLICTING DATA IN " << file1 << " AND " << file2 << ".\
			\nTHE " << conflict <<" ARE INCOMPATIBLE." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------
void MITReM::errorZero(const std::string &file, const std::string &parameter) const
{
	std::cout << "ERROR IN " << file << ".\n" << parameter << " CANNOT BE <= 0." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------

