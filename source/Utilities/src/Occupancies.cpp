/**
 *  Occupancies calculates the occupancy in the LumiCal or BeamCal given a list
 *  of background files produced by the ReadBeamCal processor
 * Arguments: [BeamCal|LumiCal] <Threshold[GeV]> <CompactFile> <Background1.root> [<Background2.root ...>]
 * The larger the number of background files the better the averaging of course
 */


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

using VD=std::vector<double>;

class VVD {
private:
  std::map<int, VD> m_vvd{};
  int m_vecSize;

public:
  VD& operator[](int selection) {
    if (m_vvd.find(selection) == m_vvd.end()) {
      m_vvd[selection] = VD(m_vecSize);
    }
    return m_vvd[selection];
  }
  explicit VVD(int vecSize) : m_vecSize(vecSize) {}

  void getOccupancy() const {
    int occupied = 0;
    int trainPads = 0;
    for (auto const& vd : m_vvd) {
      for (auto const& pad : vd.second) {
        occupied += pad;
        trainPads += 1;
      }
    }
    std::cout << "Occupied Per Train: " << double(occupied) / double(m_vvd.size()) << std::endl;
    std::cout << "Train Occupancy is: " << double(occupied) / double(trainPads) << std::endl;
    std::cout << "NTrains: " << m_vvd.size() << std::endl;
  }
};

struct Parameters {
  double minZRange = 1e-4;
  double maxZRange = 1.0;
  double threshold = 0.0;
  int bxPerTrain = 0.0;
  std::string detectorName = "BeamCal";
  std::string compactFile = "";
  std::vector<std::string> backgroundFiles{};
};

void calculateOccupancy(const Parameters& par);
void readBackgroundFile(std::string const& backgroundFile, VD& countLeft, VD& countRight, VVD& trainLeft, VVD& trainRight,
                        int bxPerTrain, double threshold, int maxEPad, int& nBX, double& totalEnergyBX, double& occupancyBX,
                        double& maxEnergy);
void drawOccupancy(const Parameters& par, BeamCalGeo const& bcg, std::string const& name, BCPadEnergies const& bcp);

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
  for (int i = 7; i < argc; ++i) {
    par.backgroundFiles.emplace_back(args[i]);
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

  const int maxFullTrainBXs(int(par.backgroundFiles.size() / par.bxPerTrain) * par.bxPerTrain);
  std::cout << "Limiting to " << maxFullTrainBXs << " BXs" << std::endl;
  int nBX(0), maxPad(0), maxEPad(0);
  int nPads = bcg.getPadsPerBeamCal();
  VD countLeft(nPads, 0), countRight(nPads, 0);
  VVD trainLeft(nPads), trainRight(nPads);
  double totalEnergy(0.0), maxOccupancy(0.0), averageOccypancy(0.0), maxEnergy(0.0);

  std::cout << "Have " << par.backgroundFiles.size() << " Files" << std::endl;
  for (auto const& backgroundFile : par.backgroundFiles) {
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

void drawOccupancy(const Parameters& par, BeamCalGeo const& bcg, std::string const& name, BCPadEnergies const& bcp) {
  TCanvas canv("cB", "cB", 800, 700);
  canv.SetRightMargin(0.21);
  canv.SetLeftMargin(0.15);
  TH2D occupancyHisto("hOcc","Occupancy;Layer;Radial Pad;Pad Occupancy [1/BX]", bcg.getBCLayers(), 0.5, bcg.getBCLayers()+0.5,
                      bcg.getBCRings(), -0.5, bcg.getBCRings());

  for (int layer = 0; layer < bcg.getBCLayers(); ++layer) {
    for (int ring = 0; ring < bcg.getBCRings() ; ++ring) {
      double averageOccupancy(0);
      for (int pad = 0; pad <  bcg.getPadsInRing(ring); ++pad) {
        averageOccupancy += bcp.getEnergy(layer, ring, pad);
      }
      averageOccupancy /= double(bcg.getPadsInRing(ring));
      occupancyHisto.Fill(layer+1, ring, averageOccupancy);
    }//for each ring
  }//for each layer
  occupancyHisto.Draw("colz");
  canv.SetLogz();
  canv.Update();
  RootUtils::MovePaletteHorizontally(&occupancyHisto, 0.0);
  occupancyHisto.GetZaxis()->SetLabelOffset(-0.01);
  occupancyHisto.GetZaxis()->SetRangeUser(par.minZRange, par.maxZRange);
  gPad->Update();
  canv.Modified();
  canv.Update();
  canv.SaveAs(Form("%s_%sOccupancies.eps", par.detectorName.c_str(), name.c_str()));
  canv.SaveAs(Form("%s_%sOccupancies.C", par.detectorName.c_str(), name.c_str()));
}

void readBackgroundFile(std::string const& backgroundFile, VD& countLeft, VD& countRight, VVD& trainLeft, VVD& trainRight,
                        int bxPerTrain, double threshold, int maxEPad, int& nBX, double& totalEnergyBX, double& occupancyBX,
                        double& maxEnergy) {
  TTree* tree;
  VD *depLeft=NULL;
  VD *depRight=NULL;

  TFile* file = TFile::Open(backgroundFile.c_str());
  if (not file) {
    std::cerr << "File not found: "<< backgroundFile  << std::endl;
    throw std::runtime_error("Failed to find file");
  }

  file->GetObject("bcTree", tree);
  if (not tree) {
    std::cerr << "Tree not found in file " << backgroundFile  << std::endl;
    file->Close();
    delete file;
    throw std::runtime_error("Failed to find tree in file");
  }

  tree->SetBranchAddress("vec_left" , &depLeft);
  tree->SetBranchAddress("vec_right", &depRight);
  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    nBX++;
    for (size_t nPad = 0; nPad < countLeft.size();++nPad) {
      totalEnergyBX += (*depLeft)[nPad];
      totalEnergyBX += (*depRight)[nPad];

      if (int(nPad) == maxEPad) {
        maxEnergy += (*depLeft)[nPad];
        maxEnergy += (*depRight)[nPad];
      }

      if((*depLeft)[nPad] > threshold) {
        countLeft[nPad] += 1;
        trainLeft[(nBX - 1) / bxPerTrain][nPad] = 1;
        occupancyBX += 1;
      }
      if((*depRight)[nPad] > threshold) {
        countRight[nPad] += 1;
        trainRight[(nBX - 1) / bxPerTrain][nPad] = 1;
        occupancyBX += 1;
      }
    }
  }
  occupancyBX /= 2.0;
  file->Close();
}
