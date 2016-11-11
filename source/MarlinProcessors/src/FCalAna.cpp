#include "FCalAna.hh"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <marlin/Global.h>
#include <marlin/ProcessorEventSeeder.h>

#include <TFile.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>

using namespace lcio;
using namespace marlin;

template <class t> inline void RotateToLumiCal(const t* vector, t* rotated, double angle) {
  const double CosAngle = cos(double(angle) / 1000.0);
  const double SinAngle = sin(double(angle) / 1000.0);
  rotated[0]            = CosAngle * vector[0] - SinAngle * vector[2];
  rotated[1]            = vector[1];
  rotated[2]            = SinAngle * vector[0] + CosAngle * vector[2];
}

template <class t> inline void RotateToLumiCal(const t& vector, t& rotated, double angle) {
  const double CosAngle = cos(double(angle) / 1000.0);
  const double SinAngle = sin(double(angle) / 1000.0);
  rotated[0]            = CosAngle * vector[0] - SinAngle * vector[2];
  rotated[1]            = vector[1];
  rotated[2]            = SinAngle * vector[0] + CosAngle * vector[2];
}

FCalAna aFCalAna;

FCalAna::FCalAna()
    : Processor("FCalAna"),
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
  _description = "FCalAna reads MCParticles and Reconstructed Particles and fills them in to a root tree";

  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollection", "Name of MCParticle Collection", m_mcParticleName,
                          std::string("MCParticle"));

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "RecoParticleCollection", "Name of ReconstructedParticle Collection",
                          m_recoParticleName, std::string("LumiCalRecoParticles"));

  registerProcessorParameter("RootFile", "Root OutputFile", m_nameOutputFile, std::string("FCalAna.root"));
}

void FCalAna::init() {
  Global::EVENTSEEDER->registerProcessor(this);

  // usually a good idea to
  printParameters();

  m_nEvt = 0;

  m_file = TFile::Open(m_nameOutputFile.c_str(), "RECREATE");
  m_tree = new TTree("fcalAna", "fcalAna");

  m_tree->Branch("thetaMC", &m_thetaMC, "thetaMC/D");
  m_tree->Branch("thetaReco", &m_thetaReco, "thetaReco/D");
  m_tree->Branch("thetaN", &m_thetaN, "thetaN/D");
  m_tree->Branch("thetaE", &m_thetaE, "thetaE/D");
  m_tree->Branch("thetaLog", &m_thetaLog);
  m_tree->Branch("eMC", &m_eMC, "eMC/D");
  m_tree->Branch("eReco", &m_eReco, "eReco/D");
  m_tree->Branch("phiMC", &m_phiMC, "phiMC/D");
  m_tree->Branch("phiReco", &m_phiReco, "phiReco/D");

}  //init

template <class T> double getTheta(T& p);
template <class T> double getTheta(T* p);

template <> double getTheta(const double* momGlob) {
  double mom[3] = {0.0, 0.0, 0.0};
  RotateToLumiCal(momGlob, mom, 10.0);
  const double r = sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
  return std::atan(r / mom[2]);
}

template <class T> double getTheta(T* p) {
  const double* momGlob = p->getMomentum();
  return getTheta(momGlob);
}

template <> double getTheta(std::vector<double> const& momGlob) { return getTheta(&momGlob[0]); }

template <class T> double getPhi(T* p) {
  const double* momGlob = p->getMomentum();
  double        mom[3]  = {0.0, 0.0, 0.0};
  RotateToLumiCal(momGlob, mom, 10.0);
  return std::atan2(mom[1], mom[0]);
}

void FCalAna::processRunHeader(LCRunHeader*) {
  //  streamlog_out (DEBUG) << "Runnumber "<< _nRun << std::endl;
  //   if(_nRun % 4 == 0) {
}

void FCalAna::processEvent(LCEvent* evt) {
  LCCollection *mcParticles = NULL, *recoParticles = NULL;

  m_nEvt++;

  try {
    mcParticles = evt->getCollection(m_mcParticleName);

  } catch (Exception& e) {
    streamlog_out(WARNING) << "Not all needed collections present " << std::endl;
    return;
  }
  m_thetaMC   = 0.0;
  m_thetaReco = 0.0;
  m_thetaN    = 0.0;
  m_thetaE    = 0.0;
  m_thetaLog.clear();
  m_eMC     = 0.0;
  m_eReco   = -1;
  m_phiMC   = 0.0;
  m_phiReco = 0.0;

  for (int iMC = 0; iMC < mcParticles->getNumberOfElements(); ++iMC) {
    MCParticle* mcp = static_cast<MCParticle*>(mcParticles->getElementAt(iMC));
    m_thetaMC       = getTheta(mcp);
    m_phiMC         = getPhi(mcp);
    m_eMC           = mcp->getEnergy();
    break;
  }

  try {
    recoParticles = evt->getCollection(m_recoParticleName);

  } catch (Exception& e) {
    streamlog_out(WARNING) << "No reconstructed particles" << std::endl;
    m_tree->Fill();
    return;
  }

  for (int iReco = 0; iReco < recoParticles->getNumberOfElements(); ++iReco) {
    ReconstructedParticle* p = static_cast<ReconstructedParticle*>(recoParticles->getElementAt(iReco));
    m_thetaReco              = getTheta(p);
    m_phiReco                = getPhi(p);
    m_eReco                  = p->getEnergy();

    auto& clusters = p->getClusters();
    for (auto* clu : clusters) {
      auto const& position  = getClusterPosition(clu);
      auto const& positionE = getEWeigthedPosition(clu);
      for (double lc = 0.1; lc < 10.0; lc += 0.1) {
        //auto const& positionL = getLogEWeigthedPosition(clu, lc);
        //m_thetaLog.push_back(getTheta(positionL));
        m_thetaLog.push_back(getThetaAverage(clu, lc));
      }
      m_thetaN = getTheta(position);
      m_thetaE = getTheta(positionE);
    }

    break;
  }

  m_tree->Fill();

  return;
}  //processEvent

template <class T> std::vector<double> FCalAna::getLogEWeigthedPosition(T* cluster, double logConstant) {
  std::vector<double> positionL(3, 0.0);

  double totalLoggedEnergy = 0.0;

  const double clusterEnergy = cluster->getEnergy();

  for (auto* hit : cluster->getCalorimeterHits()) {
    double logWeight = std::max(0.0, logConstant + std::log(double(hit->getEnergy() / clusterEnergy)));
    positionL[0] += hit->getPosition()[0] * logWeight;
    positionL[1] += hit->getPosition()[1] * logWeight;
    positionL[2] += hit->getPosition()[2] * logWeight;
    totalLoggedEnergy += logWeight;
  }  //for all hits

  positionL[0] /= totalLoggedEnergy;
  positionL[1] /= totalLoggedEnergy;
  positionL[2] /= totalLoggedEnergy;

  return positionL;
}

template <class T> double FCalAna::getThetaAverage(T* cluster, double logConstant) {
  double       totalLoggedEnergy(0.0);
  double       theta(0.0);
  const double clusterEnergy = cluster->getEnergy();
  for (auto* hit : cluster->getCalorimeterHits()) {
    double logWeight = std::max(0.0, logConstant + std::log(double(hit->getEnergy() / clusterEnergy)));
    if (not(logWeight > 0.0))
      continue;
    float posLoc[3] = {0.0, 0.0, 0.0};
    RotateToLumiCal(hit->getPosition(), posLoc, 10.0);
    // const double r2 = hit->getPosition()[0]*hit->getPosition()[0]+hit->getPosition()[1]*hit->getPosition()[1];
    // const auto thetaTemp = atan( sqrt(r2) / fabs(hit->getPosition()[2]));

    const double r         = sqrt(posLoc[0] * posLoc[0] + posLoc[1] * posLoc[1]);
    const auto   thetaTemp = atan(r / fabs(posLoc[2]));
    // logWeight *= 100.0/r;
    theta += thetaTemp * logWeight;
    totalLoggedEnergy += logWeight;
  }  //for all hits

  theta /= totalLoggedEnergy;
  return theta;
}

template <class T> std::vector<double> FCalAna::getClusterPosition(T* cluster) {
  std::vector<double> position(3, 0.0);

  int numberOfHits = 0;

  for (auto* hit : cluster->getCalorimeterHits()) {
    position[0] += hit->getPosition()[0];
    position[1] += hit->getPosition()[1];
    position[2] += hit->getPosition()[2];
    numberOfHits++;
  }

  position[0] /= double(numberOfHits);
  position[1] /= double(numberOfHits);
  position[2] /= double(numberOfHits);

  return position;
}

template <class T> std::vector<double> FCalAna::getEWeigthedPosition(T* cluster) {
  std::vector<double> positionW(3, 0.0);
  double              energyUsed = 0.0;
  for (auto* hit : cluster->getCalorimeterHits()) {
    double energy = hit->getEnergy();
    //cut off to remove very low energy hits
    if (energy > 0.3) {
      positionW[0] += hit->getPosition()[0] * energy;
      positionW[1] += hit->getPosition()[1] * energy;
      positionW[2] += hit->getPosition()[2] * energy;
      energyUsed += energy;
    }
  }  //for all hits

  positionW[0] /= energyUsed;
  positionW[1] /= energyUsed;
  positionW[2] /= energyUsed;

  return positionW;
}

void FCalAna::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void FCalAna::end() {
  streamlog_out(MESSAGE) << __PRETTY_FUNCTION__ << " " << name() << " processed " << m_nEvt << " events." << std::endl;

  m_tree->Write();
  m_file->Write();
  m_file->Close();

  delete m_file;
}
