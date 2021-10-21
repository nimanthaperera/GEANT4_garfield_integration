#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
	fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
	G4LogicalVolume * volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();

	if (volume != fScoringVolume)
		return;

	G4double edep = step->GetTotalEnergyDeposit();
	fEventAction->AddEdep(edep);

	G4int parentID = step->GetTrack()->GetParentID() ;
	G4int particlePDG = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

	if(parentID != 0 && particlePDG == 11 && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {//secondary electrons exiting FR4 layer
		if(step->GetStepLength()>0.){
			G4double xposition = step->GetPreStepPoint()->GetPosition().x();
			G4double yposition = step->GetPreStepPoint()->GetPosition().y();
			G4double zposition = step->GetPreStepPoint()->GetPosition().z();
			G4double xmomentum = step->GetPreStepPoint()->GetMomentumDirection().x();
			G4double ymomentum = step->GetPreStepPoint()->GetMomentumDirection().y();
			G4double zmomentum = step->GetPreStepPoint()->GetMomentumDirection().z();
			G4double kenergy = step->GetPreStepPoint()->GetKineticEnergy();
			//G4double tenergy = step->GetPreStepPoint()->GetTotalEnergy();
			//G4double diff = tenergy/MeV-kenergy/MeV;
			//G4cout<<"Total "<< tenergy/eV << "  KE "<<kenergy/eV<< "diff "  << diff/MeV << G4endl;
			G4double time = step->GetPreStepPoint()->GetGlobalTime();
			//Local time - start when particle is created
			//global time - reset to 0 when new event starts
			
			G4int evt =  G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
			
			G4AnalysisManager *man = G4AnalysisManager::Instance();

			man->FillNtupleIColumn(0, 0, evt);
			man->FillNtupleDColumn(0, 1, xposition/cm);
			man->FillNtupleDColumn(0, 2, yposition/cm);
			man->FillNtupleDColumn(0, 3, zposition/cm);
			man->FillNtupleDColumn(0, 4, xmomentum);
			man->FillNtupleDColumn(0, 5, ymomentum);
			man->FillNtupleDColumn(0, 6, zmomentum);
			man->FillNtupleDColumn(0, 7, time/ns); //
			man->FillNtupleDColumn(0, 8, kenergy/eV);
			man->AddNtupleRow(0);

		}
	}

	if(particlePDG == 22 && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && step->GetStepLength()>0.){//all photons exiting FR4 layer
	
	}

}

