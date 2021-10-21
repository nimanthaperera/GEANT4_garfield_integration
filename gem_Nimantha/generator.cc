#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
	fParticleGun = new G4ParticleGun(1);  // 1 Number of particles per event


	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName="gamma";
	G4ParticleDefinition *particle = particleTable->FindParticle("gamma");

	G4ThreeVector pos(0., 0., 0.);
	G4ThreeVector mom(0., 0., -1.);

	fParticleGun->SetParticlePosition(pos);
 	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleMomentum(500.*keV);
	fParticleGun->SetParticleDefinition(particle);

}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
	fParticleGun->GeneratePrimaryVertex(anEvent);

}
