/**
 *  Occupancies calculates the occupancy in the LumiCal or BeamCal given a list
 *  of background files produced by the ReadBeamCal processor
 * Arguments: [BeamCal|LumiCal] <Threshold[GeV]> <CompactFile> <Background1.root> [<Background2.root ...>]
 * The larger the number of background files the better the averaging of course
 */

#include "OccupancyUtilities.hh"

#include <BeamCalGeoDD.hh>
#include <BCPadEnergies.hh>
#include <BeamCal.hh>
#include <RootUtils.hh>

#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TPaletteAxis.h>

#include <DD4hep/Detector.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>

void calculateOccupancy(const Parameters& par);

int main (int argc, char **args) {
  if (argc < 6) {
    std::cout << "Not enough arguments "
              << "Occupancies [BeamCal|LumiCal] <Threshold[GeV]> MinZ MaxZ BXPerTrain compactFile  Background1.root "
                 "[Background2.root ...]"
              << std::endl;
    return 1;
  }

  Parameters par;

  par.detectorName = std::string(args[1]);
  par.threshold = std::atof(args[2]);
  par.minZRange = std::atof(args[3]);
  par.maxZRange = std::atof(args[4]);
  par.bxPerTrain = std::atoi(args[5]);
  par.compactFile = std::string(args[6]);
  par.files.emplace_back(std::vector<std::string>());
  for (int i = 7; i < argc; ++i) {
    par.files[0].emplace_back(args[i]);
  }


  RootUtils::SetStyle();

  try{
    calculateOccupancy(par);
  } catch(std::exception &e) {
    std::cout << "Exception " << e.what()  << std::endl;
    return 1;
  }

}

void calculateOccupancy(Parameters const& par) {
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  theDetector.fromCompact(par.compactFile);

  BeamCalGeoDD bcg(theDetector, par.detectorName, par.detectorName + "Collection");

  const int maxFullTrainBXs(int(par.files[0].size() / par.bxPerTrain) * par.bxPerTrain);
  std::cout << "******************************************************************************" << std::endl;
  std::cout << par.detectorName << "  " << par.files[0][0] << std::endl;
  std::cout << "Found " << par.files[0].size() << " files" << std::endl;
  std::cout << "Limiting to " << maxFullTrainBXs << " BXs" << std::endl;

  int nBX(0), maxPad(0), maxEPad(0);
  int nPads = bcg.getPadsPerBeamCal();
  VD countLeft(nPads, 0), countRight(nPads, 0);
  VVD trainLeft(nPads), trainRight(nPads);
  double totalEnergy(0.0), maxOccupancy(0.0), averageOccypancy(0.0), maxEnergy(0.0);

  for (auto const& backgroundFile : par.files[0]) {
    double totalEnergyBX = 0.0;
    double occupancyBX = 0;
    readBackgroundFile(backgroundFile, countLeft, countRight, trainLeft, trainRight, par.bxPerTrain, par.threshold, maxEPad,
                       nBX, totalEnergyBX, occupancyBX, maxEnergy);
    totalEnergy += totalEnergyBX;
    averageOccypancy += occupancyBX;

    if (nBX == maxFullTrainBXs)
      break;
  }

  std::cout << "Have " << nBX << " bunch crossings"  << std::endl;
  std::cout << "Have " << countLeft.size() << " pads"  << std::endl;

  //Calculate average occupancy given threshold
  totalEnergy /= (double(nBX)*2.0); //counting both sides
  for (size_t nPad = 0; nPad < countLeft.size();++nPad) {
    if (countLeft[nPad] > maxOccupancy || countRight[nPad] > maxOccupancy) {
      maxOccupancy = std::max(countLeft[nPad], std::max(countRight[nPad], maxOccupancy));
      maxPad = nPad;
    }
    countLeft[nPad] /= double(nBX);
    countRight[nPad] /= double(nBX);
  }

  BCPadEnergies left(bcg), right(bcg);

  averageOccypancy /= double(nPads);
  averageOccypancy /= double(nBX);

  left.setEnergies(countLeft);
  right.setEnergies(countRight);

  drawOccupancy(par, bcg, "Left", left);
  drawOccupancy(par, bcg, "Right", right);

  std::cout << "******************************************************************************" << std::endl;
  std::cout << "\t\t\t" << par.detectorName << std::endl;
  std::cout << "Total Energy deposit per BX in _one_ Detector: " << totalEnergy << " GeV" << std::endl;
  std::cout << "Total Energy deposit in pad " << maxEPad << ": " << maxEnergy << " GeV" << std::endl;

  std::cout << "Maximum occupancy in pad " << maxPad << " with " << maxOccupancy << " entries"
            << " in " << nBX << " bunch crossings" << std::endl;
  std::cout << "Maximum occupancy in pad " << maxPad << " with " << 100.0 * maxOccupancy / double(nBX) << "% of BX"
            << std::endl;
  std::cout << "Average occupancy " << 100.0 * averageOccypancy << "% of pads per BX" << std::endl;

  trainLeft.getOccupancy();
  trainRight.getOccupancy();
  std::cout << "******************************************************************************" << std::endl;
}
