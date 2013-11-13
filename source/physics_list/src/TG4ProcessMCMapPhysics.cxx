// $Id$

//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007, 2008 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: vmc@pcroot.cern.ch
//-------------------------------------------------

/// \file TG4ProcessMCMapPhysics.cxx
/// \brief Implementation of the TG4ProcessMCMapPhysics class 
///
/// \author I. Hrivnacova; IPN, Orsay

#include "TG4ProcessMCMapPhysics.h"
#include "TG4ProcessMCMap.h"
#include "TG4Globals.h"

#include <TVirtualMCDecayer.h>
#include <TVirtualMC.h>

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>

//_____________________________________________________________________________
TG4ProcessMCMapPhysics::TG4ProcessMCMapPhysics(const G4String& name)
  : TG4VPhysicsConstructor(name)
{
/// Standard constructor

  FillMap();
}

//_____________________________________________________________________________
TG4ProcessMCMapPhysics::TG4ProcessMCMapPhysics(G4int verboseLevel,
                                               const G4String& name)
  : TG4VPhysicsConstructor(name, verboseLevel)
{
/// Standard constructor

  FillMap();
}

//_____________________________________________________________________________
TG4ProcessMCMapPhysics::~TG4ProcessMCMapPhysics() 
{
/// Destructor
}

//
// private methods
//

//_____________________________________________________________________________
void TG4ProcessMCMapPhysics::FillMap()
{
/// Fill the process codes map with known G4 process names

  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();

  mcMap->Add("msc", kPMultipleScattering); 
  mcMap->Add("muMsc", kPMultipleScattering); 
  
  mcMap->Add("CoulombScat", kPCoulombScattering); 
  mcMap->Add("eCoulombScat", kPCoulombScattering); 
  mcMap->Add("nuclearStopping", kPCoulombScattering); 

  mcMap->Add("eIoni", kPEnergyLoss); 
  mcMap->Add("muIoni", kPEnergyLoss); 
  mcMap->Add("hIoni", kPEnergyLoss); 
  mcMap->Add("ionIoni", kPEnergyLoss); 

  mcMap->Add("Decay", kPDecay); 
  mcMap->Add("RadioactiveDecay", kPDecay); 

  mcMap->Add("compt", kPCompton); 
  mcMap->Add("phot", kPPhotoelectric); 
  
  mcMap->Add("eBrem", kPBrem); 
  mcMap->Add("muBrems", kPBrem); 
  mcMap->Add("hBrems", kPBrem); 

  mcMap->Add("annihil", kPAnnihilation); 

  mcMap->Add("conv", kPPair); 
  mcMap->Add("muPairProd", kPPair); 
  mcMap->Add("hPairProd", kPPair); 

  mcMap->Add("nCapture", kPNCapture);
  mcMap->Add("HadronCapture", kPNCapture);
  mcMap->Add("nFission", kPNuclearFission);
  mcMap->Add("HadronFission", kPNuclearFission);
  
  mcMap->Add("PiMinusAbsorptionAtRest", kPNuclearAbsorption); 
  mcMap->Add("PiMinusAbsorptionBertini", kPNuclearAbsorption); 
  mcMap->Add("PionMinusAbsorptionAtRest", kPNuclearAbsorption); 
  mcMap->Add("KaonMinusAbsorption", kPNuclearAbsorption); 
  mcMap->Add("KaonMinusAbsorptionAtRest", kPNuclearAbsorption); 
  mcMap->Add("CHIPSNuclearCaptureAtRest", kPNuclearAbsorption); 
  mcMap->Add("FTFNuclearCaptureAtRest", kPNuclearAbsorption); 
  mcMap->Add("hFritiofCaptureAtRest", kPNuclearAbsorption); 
  mcMap->Add("hBertiniCaptureAtRest", kPNuclearAbsorption); 
  
  mcMap->Add("AntiProtonAnnihilationAtRest", kPPbarAnnihilation); 
  mcMap->Add("AntiNeutronAnnihilationAtRest", kPNbarAnnihilation);

  mcMap->Add("hadElastic", kPHElastic); 
  mcMap->Add("HadronElastic", kPHElastic); 
  mcMap->Add("hElastic", kPHElastic); 
  mcMap->Add("hElasticLHEP", kPHElastic); 
  mcMap->Add("ionElastic", kPHElastic); 
  mcMap->Add("CHIPSElasticScattering", kPHElastic); 
  
  mcMap->Add("inelastic", kPHInhelastic); 
  mcMap->Add("PionMinusInelastic", kPHInhelastic); 
  mcMap->Add("pi-Inelastic", kPHInhelastic); 
  mcMap->Add("PionPlusInelastic", kPHInhelastic); 
  mcMap->Add("pi+Inelastic", kPHInhelastic); 
  mcMap->Add("KaonPlusInelastic", kPHInhelastic); 
  mcMap->Add("kaon+Inelastic", kPHInhelastic); 
  mcMap->Add("KaonMinusInelastic", kPHInhelastic); 
  mcMap->Add("kaon-Inelastic", kPHInhelastic); 
  mcMap->Add("KaonZeroLInelastic", kPHInhelastic); 
  mcMap->Add("kaon0LInelastic", kPHInhelastic); 
  mcMap->Add("KaonZeroSInelastic", kPHInhelastic); 
  mcMap->Add("kaon0SInelastic", kPHInhelastic); 
  mcMap->Add("LambdaInelastic", kPHInhelastic); 
  mcMap->Add("lambdaInelastic", kPHInhelastic); 
  mcMap->Add("ProtonInelastic", kPHInhelastic); 
  mcMap->Add("protonInelastic", kPHInhelastic); 
  mcMap->Add("AntiProtonInelastic", kPHInhelastic); 
  mcMap->Add("anti_protonInelastic", kPHInhelastic); 
  mcMap->Add("NeutronInelastic", kPHInhelastic); 
  mcMap->Add("neutronInelastic", kPHInhelastic); 
  mcMap->Add("AntiNeutronInelastic", kPHInhelastic); 
  mcMap->Add("anti_neutronInelastic", kPHInhelastic); 
  mcMap->Add("LambdaInelastic", kPHInhelastic); 
  mcMap->Add("AntiLambdaInelastic", kPHInhelastic); 
  mcMap->Add("anti_lambdaInelastic", kPHInhelastic); 
  mcMap->Add("SigmaMinusInelastic", kPHInhelastic); 
  mcMap->Add("sigma-Inelastic", kPHInhelastic); 
  mcMap->Add("AntiSigmaMinusInelastic", kPHInhelastic); 
  mcMap->Add("anti_sigma-Inelastic", kPHInhelastic); 
  mcMap->Add("SigmaPlusInelastic", kPHInhelastic); 
  mcMap->Add("sigma+Inelastic", kPHInhelastic); 
  mcMap->Add("AntiSigmaPlusInelastic", kPHInhelastic); 
  mcMap->Add("anti_sigma+Inelastic", kPHInhelastic); 
  mcMap->Add("sigma0Inelastic", kPHInhelastic); 
  mcMap->Add("XiMinusInelastic", kPHInhelastic); 
  mcMap->Add("xi-Inelastic", kPHInhelastic); 
  mcMap->Add("AntiXiMinusInelastic", kPHInhelastic); 
  mcMap->Add("anti_xi-Inelastic", kPHInhelastic); 
  mcMap->Add("XiZeroInelastic", kPHInhelastic); 
  mcMap->Add("xi0Inelastic", kPHInhelastic); 
  mcMap->Add("AntiXiZeroInelastic", kPHInhelastic); 
  mcMap->Add("anti_xi0Inelastic", kPHInhelastic); 
  mcMap->Add("OmegaMinusInelastic", kPHInhelastic); 
  mcMap->Add("omega-Inelastic", kPHInhelastic); 
  mcMap->Add("AntiOmegaMinusInelastic", kPHInhelastic); 
  mcMap->Add("anti_omega-Inelastic", kPHInhelastic); 
  mcMap->Add("DeuteronInelastic", kPHInhelastic); 
  mcMap->Add("dInelastic", kPHInhelastic); 
  mcMap->Add("AntiDeuteronInelastic", kPHInhelastic); 
  mcMap->Add("anti_deuteronInelastic", kPHInhelastic); 
  mcMap->Add("TritonInelastic", kPHInhelastic); 
  mcMap->Add("tInelastic", kPHInhelastic); 
  mcMap->Add("AntiTritonInelastic", kPHInhelastic); 
  mcMap->Add("AntiTritonInelasticProcess", kPHInhelastic); 
  mcMap->Add("anti_tritonInelastic", kPHInhelastic); 
  mcMap->Add("AlphaInelastic", kPHInhelastic); 
  mcMap->Add("alphaInelastic", kPHInhelastic); 
  mcMap->Add("AntiAlphaInelasticProcess", kPHInhelastic); 
  mcMap->Add("anti_alphaInelastic", kPHInhelastic); 
  mcMap->Add("He3Inelastic", kPHInhelastic); 
  mcMap->Add("AntiHe3InelasticProcess", kPHInhelastic); 
  mcMap->Add("anti_He3Inelastic", kPHInhelastic); 
  mcMap->Add("ionInelastic", kPHInhelastic); 
  mcMap->Add("CHIPS_LowEnergyIonIonInelastic", kPHInhelastic); 
  mcMap->Add("CHIPS_IonIonElasticScattering", kPHInhelastic); 
  mcMap->Add("hInelastic", kPHInhelastic); 
  mcMap->Add("CHIPS_Inelastic", kPHInhelastic); 
  mcMap->Add("PhotonInelastic", kPPhotonInhelastic); 

  mcMap->Add("nKiller", kPHadronic); 

  mcMap->Add("muNucl", kPMuonNuclear); 
  mcMap->Add("muMinusCaptureAtRest", kPMuonNuclear); 
  mcMap->Add("PositronNuclear", kPPositronNuclear); 
  mcMap->Add("ElectroNuclear", kPElectronNuclear); 
  mcMap->Add("photoNuclear", kPPhotoNuclear);
  
  mcMap->Add("Cerenkov", kPCerenkov);
  mcMap->Add("Scintillation", kPScintillation);
  mcMap->Add("OpAbsorption", kPLightAbsorption);
  mcMap->Add("OpRayleigh", kPRayleigh);
  mcMap->Add("Rayl", kPRayleigh);
  mcMap->Add("OpBoundary", kPLightScattering);
  mcMap->Add("OpMieHG", kPLightScattering);
  mcMap->Add("OpWLS", kPNull);
             /// \todo Add kPOpticalWavelengthShifting

  mcMap->Add("SynRad", kPSynchrotron);
  mcMap->Add("CHIPS_SynchrotronRadiation", kPSynchrotron);

  mcMap->Add("Transportation", kPTransportation);
  
  mcMap->Add("StepLimiter", kStepMax);
  mcMap->Add("specialCutForChargedHadron", kPStop); 
  mcMap->Add("specialCutForElectron", kPStop); 
  mcMap->Add("specialCutForEplus", kPStop); 
  mcMap->Add("specialCutForGamma", kPStop); 
  mcMap->Add("specialCutForMuon", kPStop); 
  mcMap->Add("specialCutForNeutralHadron", kPStop); 
  mcMap->Add("specialCutForForOther", kPStop); 
  mcMap->Add("specialCutForElectron", kPStop); 
  mcMap->Add("specialCutForElectron", kPStop); 
  mcMap->Add("MinEkineCuts", kPStop); 
  mcMap->Add("MaxTimeCuts", kPStop); 
  mcMap->Add("stackPopper", kPUserDefined);   
}  
//
// protected methods
//

//_____________________________________________________________________________
void TG4ProcessMCMapPhysics::ConstructParticle()
{
/// Instantiate particles - nothing to be done here
}

//_____________________________________________________________________________
void TG4ProcessMCMapPhysics::ConstructProcess()
{
/// Loop over all particles and their processes and check if
/// the process is present in the map

  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  G4bool success = true;

  aParticleIterator->reset();
  while ((*aParticleIterator)())
  {
    G4ProcessVector* processVector 
      = aParticleIterator->value()->GetProcessManager()->GetProcessList();
    for (G4int i=0; i<processVector->length(); i++) {
    
      G4String processName = (*processVector)[i]->GetProcessName();

      if ( mcMap->GetMCProcess(processName) == kPNoProcess &&
           processName != "photoNuclear" ) {
        G4String text = "Unknown process code for ";
        text += processName;
        TG4Globals::Warning(
          "TG4ProcessMCMapPhysics", "ConstructProcess", text);
        success = false;  
      }
    }              
  }
    
  if ( VerboseLevel() > 0 ) {
    if ( success )    
      G4cout << "### Processes mapped to VMC codes ok." << G4endl;
    else 
      G4cout << "### Processes mapped to VMC codes: some proceesses unknown." 
             << G4endl;
  }    
}
