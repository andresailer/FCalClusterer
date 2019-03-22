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
void readBackgroundFiles(StringVec const& backgroundFiles, VD& countLeft, VD& countRight, VVD& trainLeft, VVD& trainRight,
                         int bxPerTrain, VVD& bxLeft, VVD& bxRight, double threshold, int maxEPad, int& nBX,
                         double& totalEnergyBX, double& occupancyBX, double& maxEnergy);

int main(int argc, char** args) {
  if (argc < 6) {
    std::cout << "Not enough arguments "
              << "Occupancies [BeamCal|LumiCal] <Threshold[GeV]> MinZ MaxZ BXPerTrain compactFile Folder1 Folder2"
              << std::endl;
    return 1;
  }

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
  par.files.emplace_back(getFilesFromFolder(par.folders[0], ".root", 950));
  par.files.emplace_back(getFilesFromFolder(par.folders[1], ".root", 1200));

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

  const int maxFullTrainBXs(std::min(int(par.files[0].size() / par.bxPerTrain) * par.bxPerTrain,
                                     int(par.files[1].size() / par.bxPerTrain) * par.bxPerTrain));
  std::cout << "Limiting to " << maxFullTrainBXs << " BXs" << std::endl;
  int nBX(0), maxPad(0), maxEPad(0);
  int nPads = bcg.getPadsPerBeamCal();
  VD countLeft(nPads, 0), countRight(nPads, 0);
  VVD trainLeft(nPads), trainRight(nPads);
  VVD bxLeft(nPads), bxRight(nPads);
  double totalEnergy(0.0), maxOccupancy(0.0), averageOccypancy(0.0), maxEnergy(0.0);

  std::cout << "Have " << par.files[0].size() << " Files1" << std::endl;
  std::cout << "Have " << par.files[1].size() << " Files2" << std::endl;
  size_t index = 0;
  for (auto const& backgroundFile : par.files[0]) {
    if (index >= par.files[0].size() or index >= par.files[1].size()) {
      std::cout << "Exhausted one of the file vectors." << std::endl;
      break;
    }
    double totalEnergyBX = 0.0;
    double occupancyBX = 0;
    readBackgroundFiles({backgroundFile, par.files[1][index]}, countLeft, countRight, trainLeft, trainRight, par.bxPerTrain,
                        bxLeft, bxRight, par.threshold, maxEPad, nBX, totalEnergyBX, occupancyBX, maxEnergy);

    totalEnergy += totalEnergyBX;
    averageOccypancy += occupancyBX;

    index++;
    if (nBX == maxFullTrainBXs)
      break;
  }

  std::cout << "Have " << nBX << " bunch crossings" << std::endl;
  std::cout << "Have " << countLeft.size() << " pads" << std::endl;

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

  BCPadEnergies left(bcg), right(bcg);

  averageOccypancy /= double(nPads);
  averageOccypancy /= double(nBX);

  left.setEnergies(countLeft);
  right.setEnergies(countRight);

  drawOccupancy(par, bcg, "Left", left);
  drawOccupancy(par, bcg, "Right", right);

  std::cout << "******************************************************************************" << std::endl;
  std::cout << "\t\t\t" << par.detectorName << std::endl;
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

  if (false) {
    int layer(0), ring(0), sector(0);
    std::cout << "| " << std::setw(8) << "pad"
              << " | " << std::setw(5) << "layer"
              << " | " << std::setw(4) << "ring"
              << " | " << std::setw(6) << "sector"
              << " | " << std::setw(8) << "Phi [deg]"
              << " | " << std::setw(8) << "R [mm]"
              << " | " << std::setw(13) << "Occupancy [%]"
              << " |" << std::endl;
    for (int pad = 0; pad < nPads; ++pad) {
      bcg.getLayerRingPad(pad, layer, ring, sector);
      std::cout << "| " << std::setw(8) << pad << " | " << std::setw(5) << layer << " | " << std::setw(4) << ring << " | "
                << std::setw(6) << sector << " | " << std::setw(9) << bcg.getPadMiddlePhi(ring, sector) << " | "
                << std::setw(8) << bcg.getPadMiddleR(ring, sector) << " | " << std::setw(13) << countLeft[pad] * 100
                << " | ";
      // for (int bx = 0; bx < nBX; ++bx) {
      //   if(bxLeft[bx][pad] > 0) {
      //     std::cout << std::setw(13) << bxLeft[bx][pad];
      //   }
      //   std::cout << " | ";
      // }
      std::cout << std::endl;
    }
  }
}

void readBackgroundFiles(StringVec const& backgroundFiles, VD& countLeft, VD& countRight, VVD& trainLeft, VVD& trainRight,
                         int bxPerTrain, VVD& bxLeft, VVD& bxRight, double threshold, int maxEPad, int& nBX,
                         double& totalEnergyBX, double& occupancyBX, double& maxEnergy) {
  VD leftTotal(countLeft.size(), 0.0), rightTotal(countRight.size(), 0.0);

  for (auto const& backgroundFile : backgroundFiles) {
    TTree* tree;
    VD* depLeft = NULL;
    VD* depRight = NULL;

    TFile* file = TFile::Open(backgroundFile.c_str());
    if (not file) {
      std::cerr << "File not found: " << backgroundFile << std::endl;
      throw std::runtime_error("Failed to find file");
    }

    file->GetObject("bcTree", tree);
    if (not tree) {
      std::cerr << "Tree not found in file " << backgroundFile << std::endl;
      file->Close();
      delete file;
      throw std::runtime_error("Failed to find tree in file");
    }

    if (tree->GetEntries() > 1) {
      throw std::runtime_error("Too many entries in the tree");
    }

    tree->SetBranchAddress("vec_left", &depLeft);
    tree->SetBranchAddress("vec_right", &depRight);
    for (int i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      for (size_t nPad = 0; nPad < countLeft.size(); ++nPad) {
        leftTotal[nPad] += (*depLeft)[nPad];
        rightTotal[nPad] += (*depRight)[nPad];
      }
    }  //for all entries
    file->Close();
  }  //for both files

  ++nBX;
  for (size_t nPad = 0; nPad < leftTotal.size(); ++nPad) {
    totalEnergyBX += leftTotal[nPad];
    totalEnergyBX += rightTotal[nPad];

    if (int(nPad) == maxEPad) {
      maxEnergy += leftTotal[nPad];
      maxEnergy += rightTotal[nPad];
    }

    if (leftTotal[nPad] > threshold) {
      countLeft[nPad] += 1;
      trainLeft[(nBX - 1) / bxPerTrain][nPad] = 1;
      bxLeft[nBX - 1][nPad] = leftTotal[nPad];
      occupancyBX += 1;
    }
    if (rightTotal[nPad] > threshold) {
      countRight[nPad] += 1;
      trainRight[(nBX - 1) / bxPerTrain][nPad] = 1;
      bxRight[nBX - 1][nPad] = rightTotal[nPad];
      occupancyBX += 1;
    }
  }  //for all pads

  occupancyBX /= 2.0;
}
