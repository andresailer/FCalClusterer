#include "FCalAna.hh"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <marlin/ProcessorEventSeeder.h>
#include <marlin/Global.h>

#include <TTree.h>
#include <TFile.h>

#include <iostream>
#include <iomanip>

using namespace lcio ;
using namespace marlin ;

FCalAna aFCalAna ;

FCalAna::FCalAna() : Processor("FCalAna"),
		     m_nameOutputFile(""),
		     m_mcParticleName(""),
		     m_recoParticleName(""),
		     m_nRun(0),
		     m_nEvt(0),
		     m_file(NULL),
		     m_tree(NULL),
		     m_thetaMC(0.0),
		     m_thetaReco(0.0),
		     m_eMC(0.0),
		     m_eReco(0.0),
		     m_phiMC(0.0),
		     m_phiReco(0.0)

 {

  // modify processor description
  _description = "FCalAna reads MCParticles and Reconstructed Particles and fills them in to a root tree" ;

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" ,
			   "Name of MCParticle Collection"  ,
			   m_mcParticleName,
			   std::string("MCParticle") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "RecoParticleCollection",
			   "Name of ReconstructedParticle Collection"  ,
			   m_recoParticleName,
			   std::string("LumiCalRecoParticles") ) ;

  registerProcessorParameter ("RootFile",
			      "Root OutputFile",
			      m_nameOutputFile,
			      std::string("FCalAna.root") ) ;

}

void FCalAna::init() {

  Global::EVENTSEEDER->registerProcessor(this);

  // usually a good idea to
  printParameters() ;

  m_nEvt = 0;

  m_file = TFile::Open( m_nameOutputFile.c_str(), "RECREATE" );
  m_tree = new TTree("fcalAna","fcalAna");

  m_tree->Branch( "thetaMC"   , &m_thetaMC  , "thetaMC/D");
  m_tree->Branch( "thetaReco" , &m_thetaReco, "thetaReco/D");
  m_tree->Branch( "eMC"       , &m_eMC	    , "eMC/D");
  m_tree->Branch( "eReco"     , &m_eReco    , "eReco/D");
  m_tree->Branch( "phiMC"     , &m_phiMC    , "phiMC/D");
  m_tree->Branch( "phiReco"   , &m_phiReco  , "phiReco/D");









}//init


template< class T > double getTheta( T* p ){

  const double* mom = p->getMomentum();

  const double r = sqrt( mom[0]*mom[0] + mom[1]*mom[1] );
  return std::atan( r / mom[2] );

}

template< class T > double getPhi( T* p ){

  const double* mom = p->getMomentum();

  return std::atan2( mom[1], mom[0] );

}


void FCalAna::processRunHeader( LCRunHeader* ) {
  //  streamlog_out (DEBUG) << "Runnumber "<< _nRun << std::endl;
  //   if(_nRun % 4 == 0) {
}

void FCalAna::processEvent( LCEvent * evt ) {
  LCCollection *mcParticles=NULL, *recoParticles=NULL;

  m_nEvt++;

  try {

    mcParticles = evt->getCollection( m_mcParticleName ) ;

  } catch (Exception &e) {
    streamlog_out( WARNING ) << "Not all needed collections present " << std::endl;
    return;
  }
  m_thetaMC = 0.0;
  m_thetaReco = 0.0;
  m_eMC = 0.0;
  m_eReco = -1;
  m_phiMC = 0.0;
  m_phiReco = 0.0;

  for(int iMC = 0; iMC < mcParticles->getNumberOfElements() ; ++iMC ) {
    MCParticle* mcp = static_cast<MCParticle*>(mcParticles->getElementAt(iMC));
    m_thetaMC = getTheta( mcp );
    m_phiMC = getPhi( mcp );
    m_eMC = mcp->getEnergy();
    break;
  }

  try {

    recoParticles= evt->getCollection( m_recoParticleName ) ;

  } catch (Exception &e) {
    streamlog_out( WARNING ) << "No reconstructed particles" << std::endl;
    m_tree->Fill();
    return;
  }


  for (int iReco = 0; iReco < recoParticles->getNumberOfElements() ;++iReco) {
    ReconstructedParticle* p = static_cast<ReconstructedParticle*>(recoParticles->getElementAt(iReco));
    m_thetaReco = getTheta( p );
    m_phiReco = getPhi( p );
    m_eReco = p->getEnergy();
    break;
  }



  m_tree->Fill();

  return;
}//processEvent



void FCalAna::check( LCEvent * ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FCalAna::end(){

  streamlog_out ( MESSAGE ) << __PRETTY_FUNCTION__ << " " << name()
			    << " processed " << m_nEvt << " events."
			    << std::endl ;

  m_tree->Write();
  m_file->Write();
  m_file->Close();

  delete m_file;
}
