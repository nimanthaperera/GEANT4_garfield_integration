#include "run.hh"

MyRunAction::MyRunAction()
{}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
	G4AnalysisManager *man = G4AnalysisManager::Instance();
    
	G4int runID = run->GetRunID();
    
	std::stringstream strRunID;
	strRunID << runID;
    
	man->OpenFile("output"+strRunID.str()+".root");
	
	man->CreateNtuple("Electrons", "Electrons");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleDColumn("fX");
	man->CreateNtupleDColumn("fY");
	man->CreateNtupleDColumn("fZ");
	man->CreateNtupleDColumn("fpx");
	man->CreateNtupleDColumn("fpy");
	man->CreateNtupleDColumn("fpz");
	man->CreateNtupleDColumn("fT");
	man->CreateNtupleDColumn("fKE");
	man->FinishNtuple(0); //0 - beacause this is ntuple number 0
	
/*	man->CreateNtuple("Hits", "Hits");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleDColumn("fX");
	man->CreateNtupleDColumn("fY");
	man->CreateNtupleDColumn("fZ");
	man->FinishNtuple(1); //1 - beacause this is ntuple number 1
*/
	man->CreateNtuple("Scoring", "Scoring");
	man->CreateNtupleDColumn("fEdep");
	man->FinishNtuple(1);

}

void MyRunAction::EndOfRunAction(const G4Run*)
{
	G4AnalysisManager *man = G4AnalysisManager::Instance();
    
	man->Write();
	man->CloseFile();
}
