#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
	fFr4Thickness = 2.e-3*m;
	
	fMessenger = new G4GenericMessenger(this, "/detector/","DetectorConstruction"); //object,folder,helptext
	fMessenger->DeclareProperty("fFr4Thickness",fFr4Thickness,"Thickness of FR4");

	DefineMaterials();
}

MyDetectorConstruction::~MyDetectorConstruction()
{}

void MyDetectorConstruction::DefineMaterials()
{	G4NistManager *nist = G4NistManager::Instance();
	worldMat = nist->FindOrBuildMaterial("G4_Galactic");


	//FR4 material
	SiO2 = new G4Material("SiO2",2.201*g/cm3,2);
	SiO2->AddElement(nist->FindOrBuildElement("Si"),1);
	SiO2->AddElement(nist->FindOrBuildElement("O"),2);

	Epoxy = new G4Material("Epoxy",1.2*g/cm3,2);
	Epoxy->AddElement(nist->FindOrBuildElement("H"),2);
	Epoxy->AddElement(nist->FindOrBuildElement("C"),2);

	//C = nist->FindOrBuildElement("C");		

	FR4 = new G4Material("FR4",1.86*g/cm3,2);
	FR4->AddMaterial(SiO2, 52.8*perCent);
	FR4->AddMaterial(Epoxy, 47.2*perCent);

	//end of FR4 definiton

	//Gas Mixture
	Ar = nist->FindOrBuildMaterial("G4_Ar");
	Ar->GetIonisation()->SetMeanExcitationEnergy(15.7596*eV);
	
	CO2 = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
	
	G4double mixtureDensity = (Ar->GetDensity() * 70/100.0 + CO2->GetDensity() * 30/100.0) ;
	ArCO2 = new G4Material("ArCO2",mixtureDensity,2,kStateGas,293.15*kelvin,1.*atmosphere) ;
	ArCO2->AddMaterial(Ar, 0.7) ;
	ArCO2->AddMaterial(CO2, 0.3) ;

	//end of gas mixture
	//
	topMat = worldMat;
	botMat = ArCO2;
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
	//world volume
	//
	G4double topThickness = 1.0*mm;
	G4double botThickness = 1.0*mm;
	G4double detThickness = 0.001*m;
	
	G4double xWorld = 1.0*m;
	G4double yWorld = 1.0*m;
	G4double zWorld = topThickness + fFr4Thickness + botThickness + detThickness;

	solidWorld = new G4Box("solidWorld",xWorld/2,yWorld/2,zWorld*2.5);//half_length

	logicWorld =  new G4LogicalVolume(solidWorld, worldMat, "logicWorld");

	physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "phyWorld", 0, false, 0, true);//rotation,positon,loic volume, physic volume name, mother volume, implement boolean, copy number(if used several times), chechoverlaps
	
	//Top
	solidTop = new G4Box("solidTop",xWorld/2,yWorld/2,topThickness/2);

	logicTop =  new G4LogicalVolume(solidTop, topMat, "logicTop");

	physTop = new G4PVPlacement(0, G4ThreeVector(0., 0., -topThickness/2), logicTop, "phyTop", logicWorld, false, 0, true);	

	//FR4
	solidFR4 = new G4Box("solidFR4",xWorld/2,yWorld/2,fFr4Thickness/2);

	logicFR4 =  new G4LogicalVolume(solidFR4, FR4, "logicFR4");

	physFR4 = new G4PVPlacement(0, G4ThreeVector(0., 0., -topThickness - fFr4Thickness/2), logicFR4, "phyFR4", logicWorld, false, 0, true);	


	//Bottom
	solidBot = new G4Box("solidBot",xWorld/2,yWorld/2,botThickness/2);

	logicBot =  new G4LogicalVolume(solidBot, botMat, "logicBot");

	physBot = new G4PVPlacement(0, G4ThreeVector(0., 0., -topThickness - fFr4Thickness - botThickness/2), logicBot, "phyTop", logicWorld, false, 0, true);	

/*	//Adding sensitive  detecor
	//
	G4int nRows = 100;
	G4int nCols = 100;

	solidDetector = new G4Box("solidDetector",xWorld/nRows, yWorld/nCols, detThickness);

	logicDetector = new G4LogicalVolume(solidDetector, botMat, "logicDetector");

	for (G4int i = 0; i<nRows; i++)
	{
		for (G4int j = 0; j<nCols; j++)
		{
			physDetector = new G4PVPlacement(0, G4ThreeVector(-0.5*m+(i+0.5)*m/nRows, -0.5*m+(j+0.5)*m/nCols, -topThickness - fFr4Thickness - botThickness - detThickness/2), logicDetector, "physDetector", logicWorld, false, j+i*nCols, true);
		}	
	}
*/
	fScoringVolume = logicBot;

	return physWorld; //should be the highest mother volume

}

/*void MyDetectorConstruction::ConstructSDandField()
{
	MySensitiveDetector* sensDet = new MySensitiveDetector("SensitiveDetector");

	logicDetector->SetSensitiveDetector(sensDet);

}*/
