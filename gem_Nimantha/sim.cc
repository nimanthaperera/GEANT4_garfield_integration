#include <iostream>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "construction.hh"
#include "physics.hh"
#include "action.hh"

int main(int argc, char** argv)
{
	#ifdef G4MULTITHREADED
		G4MTRunManager *runManager = new G4MTRunManager();
	#else
		G4RunManager *runManager = new G4RunManager();
	#endif

	runManager->SetUserInitialization(new MyDetectorConstruction());//2
	runManager->SetUserInitialization(new MyPhysicsList());//3
	runManager->SetUserInitialization(new MyActionInitialization());

	//runManager->Initialize(); //Removed and added to run.mac for multithreading - should be added in other mac files as well
	//
	
	G4UIExecutive *ui = 0;

	if (argc == 1) //Executeif only the program name is given
	{
		ui = new G4UIExecutive(argc, argv);
	}

	G4VisManager *visManager = new G4VisExecutive();
	visManager->Initialize();

	G4UImanager *UImanager = G4UImanager::GetUIpointer();
	if(ui)
	{
		
		//UImanager->ApplyCommand("/vis/open OGL");
		//UImanager->ApplyCommand("/vis/drawVolume");
		//UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1");
		//UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
		//UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
		//UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");	

		UImanager->ApplyCommand("/control/execute vis.mac"); //the vis.mac contains the above commands
	
		ui->SessionStart();
	}
	else
	{
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager-> ApplyCommand(command+fileName);
	}
	delete ui;

	delete runManager;
	delete visManager;
	return 0;
}
