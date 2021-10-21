#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"

#include "detector.hh"

#include "G4GenericMessenger.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	MyDetectorConstruction();
	~MyDetectorConstruction();

	G4LogicalVolume *GetScoringVolume() const {return fScoringVolume;}

	virtual G4VPhysicalVolume *Construct();
private:
	G4LogicalVolume *logicDetector;
	//virtual void ConstructSDandField();

	G4GenericMessenger *fMessenger;
	G4double fFr4Thickness;
	G4Box *solidWorld, *solidTop, *solidFR4, *solidBot, *solidDetector;
	G4LogicalVolume *logicWorld, *logicTop, *logicFR4, *logicBot;
	G4VPhysicalVolume *physWorld, *physTop, *physFR4, *physBot, *physDetector;

	G4Material *SiO2, *Epoxy, *FR4, *worldMat, *Ar, *CO2, *ArCO2, *topMat, *botMat;
	//G4Element *Ar; //*C;	
	void DefineMaterials();

	G4LogicalVolume *fScoringVolume;
};







#endif

