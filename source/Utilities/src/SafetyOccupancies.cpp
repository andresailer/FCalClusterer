/**
 *  MergedOccupancies calculates the occupancy in the LumiCal or BeamCal given a list
 *  of two background folders produced by the ReadBeamCal processor
 * Arguments: [BeamCal|LumiCal] <Threshold[GeV]> MinZ MaxZ nBX <CompactFile> folder1 folder2
 * The larger the number of background files the better the averaging of course
 * files are capped at integer number of bunch trains
 */

#include "OccupancyUtilities.hh"

#include <BCPadEnergies.hh>
#include <BeamCal.hh>
#include <BeamCalGeoDD.hh>
#include <RootUtils.hh>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TTree.h>

#include <DD4hep/Detector.h>

#include <dirent.h>

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void calculateOccupancy(const Parameters& par);

void getAveragePadOccupancyPerTrain(Parameters const& par, VVD& trainLeft, VVD& trainRight, int sourceIndex,
                                    int const nPads) {
  auto const& backgroundFiles = par.files[sourceIndex];
  const int maxFullTrainBXs(int(backgroundFiles.size() / par.bxPerTrain) * par.bxPerTrain);

  VD countLeft(nPads, 0), countRight(nPads, 0);
  int nBX(0), maxPad(0), maxEPad(0);
  double totalEnergy(0.0), maxOccupancy(0.0), averageOccypancy(0.0), maxEnergy(0.0);

  std::cout << "******************************************************************************" << std::endl;
  std::cout << par.detectorName << "  " << par.folders[sourceIndex] << std::endl;

  std::cout << "Found " << backgroundFiles.size() << " files" << std::endl;
  std::cout << "Limiting to " << maxFullTrainBXs << " BXs" << std::endl;

  for (auto const& backgroundFile : backgroundFiles) {
    double totalEnergyBX = 0.0;
    double occupancyBX = 0;
    readBackgroundFile(backgroundFile, countLeft, countRight, trainLeft, trainRight, par.bxPerTrain, par.threshold, maxEPad,
                       nBX, totalEnergyBX, occupancyBX, maxEnergy);

    totalEnergy += totalEnergyBX;
    averageOccypancy += occupancyBX;

    if (nBX == maxFullTrainBXs)
      break;

  }  //for each file
  //Calculate average occupancy given threshold
  totalEnergy /= (double(nBX));  //counting both sides
  maxEnergy /= (double(nBX));    //counting both sides
  for (size_t nPad = 0; nPad < countLeft.size(); ++nPad) {
    if (countLeft[nPad] > maxOccupancy || countRight[nPad] > maxOccupancy) {
      maxOccupancy = std::max(countLeft[nPad], std::max(countRight[nPad], maxOccupancy));
      maxPad = nPad;
    }
    countLeft[nPad] /= double(nBX);
    countRight[nPad] /= double(nBX);
  }

  averageOccypancy /= double(nPads);
  averageOccypancy /= double(nBX);

  std::cout << "Total Energy deposit per BX in both Detectors: " << totalEnergy << " GeV" << std::endl;
  std::cout << "Energy deposit per BX in pad " << maxEPad << ": " << maxEnergy << " GeV" << std::endl;
  std::cout << "Maximum occupancy in pad " << maxPad << " with " << maxOccupancy << " entries"
            << " in " << nBX << " bunch crossings" << std::endl;
  std::cout << "Maximum occupancy in pad " << maxPad << " with " << 100.0 * maxOccupancy / double(nBX) << "% of BX"
            << std::endl;
  std::cout << "Average occupancy " << 100.0 * averageOccypancy << "% of pads per BX" << std::endl;
  trainLeft.getOccupancy();
  trainRight.getOccupancy();
  std::cout << "******************************************************************************" << std::endl;
}

int main(int argc, char** args) {
  if (argc < 6) {
    std::cout << "Not enough arguments "
              << "Occupancies [BeamCal|LumiCal] <Threshold[GeV]> MinZ MaxZ BXPerTrain compactFile Folder1 Folder2"
              << std::endl;
    return 1;
  }
  bool debug = false;
  Parameters par;
  int argCount = 1;
  par.detectorName = std::string(args[argCount++]);
  par.threshold = std::atof(args[argCount++]);
  par.minZRange = std::atof(args[argCount++]);
  par.maxZRange = std::atof(args[argCount++]);
  par.bxPerTrain = std::atoi(args[argCount++]);
  par.compactFile = std::string(args[argCount++]);
  par.folders.emplace_back(std::string(args[argCount++]));
  par.folders.emplace_back(std::string(args[argCount++]));
  par.files.emplace_back(getFilesFromFolder(par.folders[0], ".root", 950, debug));
  par.files.emplace_back(getFilesFromFolder(par.folders[1], ".root", 1200, debug));

  RootUtils::SetStyle();

  try {
    calculateOccupancy(par);
  } catch (std::exception& e) {
    std::cout << "Exception " << e.what() << std::endl;
    return 1;
  }
}

void calculateOccupancy(Parameters const& par) {
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  theDetector.fromCompact(par.compactFile);
  BeamCalGeoDD bcg(theDetector, par.detectorName, par.detectorName + "Collection");
  int nPads = bcg.getPadsPerBeamCal();

  std::vector<VVD> trainsLeft(par.files.size(), VVD(nPads)), trainsRight(par.files.size(), VVD(nPads));

  for (size_t sourceIndex = 0; sourceIndex < par.files.size(); ++sourceIndex) {
    getAveragePadOccupancyPerTrain(par, trainsLeft[sourceIndex], trainsRight[sourceIndex], sourceIndex, nPads);
  }  //for each source

  double averageSafetiedOccupancy(0.0);
  double trainCount = 0;
  for (auto const& trains : {trainsLeft, trainsRight}) {
    trainCount += 1;
    for (int j = 0; j < nPads; ++j) {
      double combinedOccupancy(1.0);
      for (size_t sourceIndex = 0; sourceIndex < par.sf.size(); ++sourceIndex) {
        double thisOccupancy = trains[sourceIndex].getPadOccupancy(j);
        double safetiedOccupancy(safetyFactored(thisOccupancy, par.sf[sourceIndex]));
        combinedOccupancy *= (1.0 - safetiedOccupancy);
      }
      averageSafetiedOccupancy += 1.0 - combinedOccupancy;
    }
  }
  std::cout << "The average occupancy per train (" << par.bxPerTrain << " BXs)"
            << " with safety factors is " << averageSafetiedOccupancy / double(nPads) * 100 / double(trainCount) << "%"
            << " for the " << par.detectorName << std::endl;
  std::cout << "******************************************************************************" << std::endl;
}
