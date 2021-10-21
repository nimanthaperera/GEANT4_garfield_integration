#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <limits>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>

#include <TFile.h>

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"

#include "TTree.h"

using namespace Garfield;

int main(int argc, char * argv[]) {
/*
  unsigned int nEvents = 1;
  if (argc > 1) {
    std::stringstream os;
    os <<argv[1];
    os >> nEvents;
  }
*/
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  const bool debug = true;

  // Load the field map.
  ComponentAnsys123* fm = new ComponentAnsys123();
  const std::string efile = "ELIST.lis";
  const std::string nfile = "NLIST.lis";
  const std::string mfile = "MPLIST.lis";
  const std::string sfile = "PRNSOL.lis";
  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  //Set the periodicities
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  //Print some information about the cell dimensions.
  fm->PrintRange();

  // Dimensions of the GEM - (garfield uses cm opposed to mm in ansys)
  const double pitch = 0.014;
  const double kapton = 50.e-4; //0.005
  const double metal = 5.e-4; //0.0005
  const double outdia = 70.e-4; //0.007
  const double middia = 50.e-4; //0.005
  const double drift = 0.3;
  const double transferx = 0.1;
  const double transfery = 0.2;
  const double induct = 0.1;
  const double rim = 0.008;
  const double xWorld = 100.;
  const double yWorld = 100.;	
  //////Draws the electric field inside the gem
  const bool plotField = false;
  if (plotField) {
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(fm);
    //Set the viewing plane (normal vector) and ranges.
    fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
    fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
    fieldView->SetVoltageRange(-160., 160.);
    TCanvas* cF = new TCanvas();
    fieldView->SetCanvas(cF);
    fieldView->PlotContour();
  }

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 70., "co2", 30.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->EnableDebugging();
  gas->Initialise();  
  gas->DisableDebugging();
  // Set the Penning transfer efficiency.
  const double rPenning = 0.57;
  //Mean distance from the point of excitation.
  const double lambdaPenning = 0.;
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  gas->LoadIonMobility("IonMobility_Ar+_Ar.txt");
  
  // Associate the gas with the corresponding field map material. 
  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }

  fm->PrintMaterials();

  // Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-xWorld/2, -yWorld/2, -kapton/2-metal-transfery-metal-kapton-metal-induct-0.1,    //sensor->SetArea(-5 * pitch, -5 * pitch, -0.03,
                   xWorld/2,  yWorld/2,  kapton/2+metal+transferx+metal+kapton+metal+drift+0.1);   //5 * pitch,  5 * pitch,  0.03);

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

/*  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  //Integrate in constant (2 micro meter) distance intervals.
  drift->SetDistanceSteps(2.e-4);

  const bool plotDrift = false;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) {
    driftView->SetArea(-2 * pitch, -2 * pitch, -0.02,
                        2 * pitch,  2 * pitch,  0.02);
    // Plot every 10 collisions (in microscopic tracking).
    aval->SetCollisionSteps(10); 
    aval->EnablePlotting(driftView);
    drift->EnablePlotting(driftView);
  }
*/
  // Histograms
  TFile* f= new TFile("ht_hole_diameter3.root","update");

  //int nBinsGain = 101;
  //double gmin =   -0.5;
  //double gmax = 100.5;
  //TH1F* hElectrons = new TH1F("hElectrons", "Number of events producing given number of electrons",
  //                            nBinsGain, gmin, gmax);
  //TH1F* hIons = new TH1F("hIons", "Number of events producing given number of ions",
  //                       nBinsGain, gmin, gmax);

  //int nBinsChrg = 100;
  //TH1F* hChrgE = new TH1F("hChrgE", "Electrons on plastic",
  //                        nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
  //TH1F* hChrgI = new TH1F("hChrgI", "Ions on plastic", 
  //                        nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);

  TH1F* hgain2 = new TH1F("hgain2","Gain Distribution (Number of events with the given gain) for a voltage difference of 500V (Drift 50V)",1000001,-0.5,1000000.5);

  TH1F* heffgain2 = new TH1F("heffgain2","Effective Gain Distribution (Number of events with the given gain) for a voltage difference of 500V (Drift 50V)",1000001,-0.5,1000000.5);
  

  double sumIonsTotal = 0.;
  double sumIonsDrift = 0.;
  double sumIonsPlastic = 0.;

  double sumElectronsTotal = 0.;
  double sumElectronsPlasticTop = 0.;
  double sumElectronsUpperMetalTop = 0.;
  double sumElectronsLowerMetalTop = 0.;
  double sumElectronsTransferx = 0.;
  double sumElectronsPlasticMiddle = 0.;
  double sumElectronsUpperMetalMiddle = 0.;
  double sumElectronsLowerMetalMiddle = 0.;
  double sumElectronsTransfery = 0.;
  double sumElectronsPlasticBottom = 0.;
  double sumElectronsUpperMetalBottom = 0.;
  double sumElectronsLowerMetalBottom = 0.;
  double sumElectronsInduct = 0.;
  double sumElectronsOther = 0.;
  float sumElectronsreadout = 0.;

  double nl_dbl_min = std::numeric_limits< double >::min();
  


  //counts the number that is out of the range
  //int neo=1000000;
  
  //Reading from geant4 tree	

  TFile* file = new TFile("output0.root");
  TTree* t = (TTree*)file->Get("Electrons");

  int fEvent, nentries;
  double fX, fY, fZ, fpx, fpy, fpz, fT, fKE;

  t->SetBranchAddress("fEvent",&fEvent);
  t->SetBranchAddress("fX",&fX);
  t->SetBranchAddress("fY",&fY);
  t->SetBranchAddress("fZ",&fZ);
  t->SetBranchAddress("fpx",&fpx);
  t->SetBranchAddress("fpy",&fpy);
  t->SetBranchAddress("fpz",&fpz);
  t->SetBranchAddress("fT",&fT);
  t->SetBranchAddress("fKE",&fKE);

  nentries = t->GetEntries();

  if(argc!=0){
  	std::stringstream os;
  	os <<argv[1];
  	os >> nentries;
  }
  std::cout<<"Number of entries: "<<nentries<<std::endl;
 

  //const int nEvents = 100;
  //for (int i = nEvents; i--;) { 
  for(int i = 0; i < nentries; i++) { 
    //if (debug || i % 10 == 0) std::cout << i << "/" << nEvents << "\n";
    // Setting the initial position
    t->GetEntry(i);
    //std::cout<<fX<<" "<<fY<<" "<<fKE<<" "<<fpx<<" "<<fpy<<" "<<fpz<<std::endl;	
    const double red = nl_dbl_min*pow(10,292);
    double temp0 = drift - red;
    std::cout<<red<<std::endl;
    double z0 = kapton/2 + metal + transferx + metal + kapton + metal + temp0;// * RndmUniform();  //double z0 = 0.025; 
    std::cout<<std::endl;
    
    //Calculate an electron avalanche
    //fX=-1.76117;fY=2.05465;fT=0.0100069;fKE=8.20311e+10;fpx=-2.20302e-05;fpy=1.94251e-05;fpz=-1;
    //std::cout<<fX<<" "<<fY<<" "<<fKE<<" "<<fpx<<" "<<fpy<<" "<<z0<<std::endl;
    aval->AvalancheElectron(fX, fY, z0, fT, fKE, fpx, fpy, fpz);
    int ne = 0, ni = 0;
    //Get the number of electrons and ions produced in the avalanche
    aval->GetAvalancheSize(ne, ni);
    //hElectrons->Fill(ne);

    hgain2 -> Fill(ne);
/*    if (ne>neo){ 
        std::cout<<ne<<std::endl;
        break;
    }
*/    
    //for (int j=ne;j>0;j--) hgain -> Fill(i);
    
    //hIons->Fill(ni);

    //Get the number of electron tracks.
    const int np = aval->GetNumberOfElectronEndpoints();
	    //Get the starting and end point of the first electron.
	    double xe1, ye1, ze1, te1, e1;
	    double xe2, ye2, ze2, te2, e2;
	    double xi1, yi1, zi1, ti1;
	    double xi2, yi2, zi2, ti2;
	    int status;
	    //Loop over the endpoints, i.e. the electron drift lines.
	    std::cout<<"Number of electrons "<<ne<<" Number of endpoints "<<np<<std::endl;
	    for (int j = np; j--;) {
	      //Get the start and the endpoint of the electrons.
	      aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
					   xe2, ye2, ze2, te2, e2, status);

	      
	      std::cout<<ze2<<std::endl;


	      sumElectronsTotal += 1.;
	      if (ze2 >= kapton / 2. + metal + transferx + metal + kapton && ze2 <= kapton / 2. + metal + transferx + metal + kapton + metal) {
		sumElectronsUpperMetalTop += 1.;
		std::cout<<"UpperMetaltop"<<std::endl;
	      } else if (ze2 > kapton / 2. + metal + transferx + metal && ze2 < kapton / 2. + metal + transferx +metal + kapton) {
		sumElectronsPlasticTop += 1.;
		std::cout<<"kaptonTop"<<std::endl;
	      } else if (ze2 >= kapton / 2. + metal + transferx && ze2 <= kapton / 2. + metal + transferx + metal) {
		sumElectronsLowerMetalTop += 1.;
		std::cout<<"lowermetalTop"<<std::endl;
	      } else if (ze2 > kapton / 2. + metal && ze2 < kapton / 2. + metal + transferx) {
		sumElectronsTransferx += 1.;
		std::cout<<"transferx"<<std::endl;
	      } else if (ze2 >= kapton / 2. && ze2 <= kapton  / 2. + metal) {
		sumElectronsUpperMetalMiddle += 1.;
		std::cout<<"uppermetalMiddle"<<std::endl;
	      } else if (ze2 > -kapton / 2. && ze2 < kapton / 2.) {
		sumElectronsPlasticMiddle += 1.;
		std::cout<<"kaptonMiddle"<<std::endl;
	      } else if (ze2 <= -kapton / 2. && ze2 >= -kapton / 2. - metal) {
		sumElectronsLowerMetalBottom += 1.;
		std::cout<<"lowermetalMiddle"<<std::endl;
	      } else if (ze2 < -kapton / 2. - metal && ze2 > -kapton / 2. - metal - transfery) {
		sumElectronsTransfery += 1.;
		std::cout<<"transfery"<<std::endl;
	      } else if (ze2 <= -kapton / 2. - metal - transfery && ze2 >= -kapton / 2. - metal - transfery - metal) {
		sumElectronsUpperMetalMiddle += 1.;
		std::cout<<"uppermetalBottom"<<std::endl;
	      } else if (ze2 < -kapton / 2. - metal - transfery - metal && ze2 > -kapton / 2. - metal - transfery - metal - kapton) {
		sumElectronsPlasticMiddle += 1.;
		std::cout<<"kaptonBottom"<<std::endl;
	      } else if (ze2 <= -kapton / 2. - metal - transfery - metal - kapton && ze2 >= -kapton / 2. - metal - transfery - metal - kapton - metal) {
		sumElectronsLowerMetalBottom += 1.;
		std::cout<<"lowermetalBottom"<<std::endl;
	      } else if (fabs(ze2+kapton/2.+metal+transfery+metal+kapton+metal+induct)<1.e-3) {
		sumElectronsreadout += 1.;
		std::cout<<"readout"<<std::endl;
	      } else if (ze2 < -kapton / 2. - metal - transfery - metal - kapton - metal) {
		sumElectronsInduct += 1.;
		std::cout<<"induct"<<std::endl;
	      } else {
		sumElectronsOther += 1.;
	      }


	/*      //Claculate an ion drift line from the creation point of each electron
	      drift->DriftIon(xe1, ye1, ze1, te1);
	      drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
				       xi2, yi2, zi2, ti2, status);
	      if (zi1 < 0.01) {
		sumIonsTotal += 1.;
		if (zi2 > 0.01) sumIonsDrift += 1.;
	      }
	      if (zi2 > -kapton / 2. && zi2 < kapton / 2.) {
		hChrgI->Fill(zi2 * 1.e4);
		sumIonsPlastic += 1.;
	      }*/


	    }
	    
	    std::cout<<std::endl<<"uppermetaltop: "<<sumElectronsUpperMetalTop<<std::endl;
	    std::cout<<std::endl<<"kaptontop: "<<sumElectronsPlasticTop<<std::endl;
	    std::cout<<std::endl<<"lowermetaltop: "<<sumElectronsLowerMetalTop<<std::endl;
	    std::cout<<std::endl<<"transferx: "<<sumElectronsTransferx<<std::endl;
	    std::cout<<std::endl<<"uppermetalmiddle: "<<sumElectronsUpperMetalMiddle<<std::endl;
	    std::cout<<std::endl<<"kaptonmiddle: "<<sumElectronsPlasticMiddle<<std::endl;
	    std::cout<<std::endl<<"lowermetalmiddle: "<<sumElectronsLowerMetalMiddle<<std::endl;
	    std::cout<<std::endl<<"transfery: "<<sumElectronsTransfery<<std::endl;
	    std::cout<<std::endl<<"uppermetalbottom: "<<sumElectronsUpperMetalBottom<<std::endl;
	    std::cout<<std::endl<<"kaptonbottom: "<<sumElectronsPlasticBottom<<std::endl;
	    std::cout<<std::endl<<"lowermetalbottom: "<<sumElectronsLowerMetalBottom<<std::endl;
	    std::cout<<std::endl<<"induct: "<<sumElectronsInduct<<std::endl;
	    std::cout<<std::endl<<"readout: "<<sumElectronsreadout<<std::endl;
	    std::cout<<std::endl<<"other: "<<sumElectronsOther<<std::endl;
	    heffgain2->Fill(sumElectronsreadout);
	    sumElectronsreadout = 0.;



	  }
        //  file->Close();
	  const double negainMean = hgain2->GetMean();
	  std::cout<<"mean: "<<negainMean<<std::endl;

	  const double negainError = hgain2->GetRMS()/sqrt(nentries);
	  std::cout<<"error: "<<negainError<<std::endl;

	  const double effgainMean = heffgain2->GetMean();
	  std::cout<<"mean (eff. gain): "<<effgainMean<<std::endl;

	  const double effgainError = heffgain2->GetRMS()/sqrt(nentries);
	  std::cout<<"error (eff. gain): "<<effgainError<<std::endl;

	/*  double fFeedback = 0.;
	  if (sumIonsTotal > 0.) fFeedback = sumIonsDrift / sumIonsTotal;
	  std::cout << "Fraction of ions drifting back: " << fFeedback << "\n"; 
	*/
	/*
	  const double neMean = hElectrons->GetMean();
	  std::cout << "Mean number of electrons: " << neMean << "\n";
	  const double niMean = hIons->GetMean();
	  std::cout << "Mean number of ions: " << niMean << "\n";
	*/
	/*  std::cout << "Mean number of electrons on plastic: "
		    << sumElectronsPlastic / nEvents << "\n";
	  std::cout << "Mean number of ions on plastic: "
		    << sumIonsPlastic / nEvents << "\n";
	 
	  std::cout << "Electron endpoints:\n";
	  const double fUpperMetal = sumElectronsUpperMetal / sumElectronsTotal;
	  const double fPlastic = sumElectronsPlastic / sumElectronsTotal;
	  const double fLowerMetal = sumElectronsLowerMetal / sumElectronsTotal;
	  const double fTransfer = sumElectronsTransfer / sumElectronsTotal;
	  const double fOther = sumElectronsOther / sumElectronsTotal;
	  std::cout << "    upper metal: " << fUpperMetal * 100. << "%\n";
	  std::cout << "    plastic:     " << fPlastic * 100. << "%\n";
	  std::cout << "    lower metal: " << fLowerMetal * 100. << "%\n";
	  std::cout << "    transfer:    " << fTransfer * 100. << "%\n";
	  std::cout << "    other:       " << fOther * 100. << "%\n";
	*/

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 /* TCanvas* cD = new TCanvas();
	  const bool plotGeo = false;
	  if (plotGeo && plotDrift) {
	    // Build the geometry in Root.
	    TGeoManager* geoman = new TGeoManager("world", "geometry");
	    TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
	    TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
	    TGeoMaterial* matKapton = new TGeoMaterial("Kapton", 12, 6, 1.42);
	    TGeoMedium* medKapton = new TGeoMedium("Kapton", 2, matKapton);
	    TGeoMaterial* matCopper = new TGeoMaterial("Copper", 63, 29, 8.94);
	    TGeoMedium* medCopper = new TGeoMedium("Copper", 3, matCopper);
	    TGeoVolume* volTop = geoman->MakeBox("TOP", 
						 medVacuum, pitch, pitch, 0.02);
	    volTop->SetVisibility(0);
	    TGeoBBox* shpKapton = new TGeoBBox("K", pitch / 2., 
						    pitch / 2., 
						    kapton / 2.);
	    TGeoPcon* shpHole = new TGeoPcon("H", 0., 360., 3);
	    shpHole->DefineSection(0, -kapton / 2., 0., outdia / 2.);
	    shpHole->DefineSection(1,           0., 0., middia / 2.);
	    shpHole->DefineSection(2,  kapton / 2., 0., outdia / 2.);

	    TGeoCompositeShape* shpGem = new TGeoCompositeShape("G", "K - H");
	    TGeoVolume* volKapton = new TGeoVolume("Kapton", shpGem, medKapton);
	    volKapton->SetLineColor(kGreen);
	    volKapton->SetTransparency(50);

	    TGeoBBox* shpMetal = new TGeoBBox("M", pitch / 2., 
						   pitch / 2., 
						   metal / 2.);
	    TGeoTube* shpTube = new TGeoTube("T", 0., outdia / 2., metal / 2.);
	    TGeoCompositeShape* shpElectrode = new TGeoCompositeShape("E", "M - T");
	    TGeoVolume* volElectrode = new TGeoVolume("Electrode", 
						      shpElectrode, medCopper);
	    volElectrode->SetLineColor(kBlue);
	    volElectrode->SetTransparency(50);

	    TGeoVolumeAssembly* volGem = new TGeoVolumeAssembly("Gem");
	    const double shift =  0.5 * (metal + kapton);
	    volGem->AddNode(volKapton, 1);
	    volGem->AddNode(volElectrode, 2, new TGeoTranslation(0., 0.,  shift));
	    volGem->AddNode(volElectrode, 3, new TGeoTranslation(0., 0., -shift));

	    volTop->AddNode(volGem, 1);
	    volTop->AddNode(volGem, 2, new TGeoTranslation(-pitch, 0., 0.));
	    volTop->AddNode(volGem, 3, new TGeoTranslation(+pitch, 0., 0.));
	    volTop->AddNode(volGem, 4, 
		       new TGeoTranslation(-pitch / 2., sqrt(3) * pitch / 2., 0.));
	    volTop->AddNode(volGem, 5, 
		       new TGeoTranslation(+pitch / 2., sqrt(3) * pitch / 2., 0.));
	    volTop->AddNode(volGem, 6,
		       new TGeoTranslation(-pitch / 2., -sqrt(3) * pitch / 2., 0.));
	    volTop->AddNode(volGem, 7,
		       new TGeoTranslation(+pitch / 2., -sqrt(3) * pitch / 2., 0.));
	    geoman->SetVerboseLevel(0);
	    geoman->SetTopVolume(volTop);
	    geoman->CloseGeometry();
	    geoman->CheckOverlaps(0.1e-4);
	    geoman->SetNmeshPoints(100000);
	    cD->cd();
	    geoman->GetTopVolume()->Draw("ogl");
	  }*/
	/////////////////////////////////////////////////////////////////////////////////////////////
	 /* if (plotDrift) {
	    driftView->SetCanvas(cD);
	    driftView->Plot(true, true);
	    std::cout<<"Plotting is complete";
	  }*/

	  const bool plotHistogram = false;
	  if (plotHistogram) {
	    TCanvas* cH = new TCanvas("cH", "Histograms", 800, 700);
	    cH->Divide(2, 3);
	    cH->cd(1);
	    //hElectrons->Draw();
	    cH->cd(2);
	    //hIons->Draw();
	    cH->cd(3);
	    //hChrgE->Draw();
	    cH->cd(4);
	    //hChrgI->Draw();
	    cH->cd(5);
	    hgain2->Draw();

	    hgain2->Write();
	    heffgain2->Write();
	    f->Close();
	    
	  }

	  app.Run(kTRUE);

	}
