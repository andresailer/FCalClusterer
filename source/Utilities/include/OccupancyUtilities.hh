/** OccupancyUtilities.hh --- 
 *
 * Copyright (C) 2019 Andre Sailer
 *
 * Author: Andre Sailer <andre.philippe.sailer@cern.ch>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 */


#ifndef OCCUPANCYUTILITIES_HH
#define OCCUPANCYUTILITIES_HH 1

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

#include <dirent.h>

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

using VD=std::vector<double>;
using StringVec = std::vector<std::string>;


using VD=std::vector<double>;

class VVD {
private:
  std::map<int, VD> m_vvd{};
  int m_vecSize;
public:
  VD& operator[](int selection){
    if(m_vvd.find(selection) == m_vvd.end()) {
      m_vvd[selection] = VD(m_vecSize);
    }
    return m_vvd[selection];
  }
  explicit VVD(int vecSize): m_vecSize(vecSize) {}

  void getOccupancy() const {
    int occupied=0;
    int trainPads = 0;
    for (auto const& vd: m_vvd) {
      for (auto const& pad: vd.second) {
        occupied += pad;
        trainPads += 1;
      }
    }
    std::cout << "Occupied Per Train: " << double(occupied)/double(m_vvd.size())
              << " / " << m_vecSize
              << std::endl;
    std::cout << "Train Occupancy is: " << double(occupied)/double(trainPads) * 100 << "%"  << std::endl;
    std::cout << "NTrains: " << m_vvd.size()  << std::endl;
  }

  double getPadOccupancy(size_t pad) const {
    double occupied = 0;
    double trainPads = 0;
    for (auto const& vd: m_vvd) {
      occupied += vd.second[pad];
      trainPads += 1;
    }
    return occupied/trainPads;
  }


};


///Apply safety factor to occupancy
double safetyFactored(double const occupancy, double const safetyFactor) {
  return 1.0 - std::pow(1.0 - occupancy, safetyFactor);
}

struct Parameters 
{     
  double minZRange=1e-4;
  double maxZRange=1.0;
  double threshold=0.0;
  int bxPerTrain=0.0;
  std::string detectorName="BeamCal";
  std::string compactFile="";
  std::vector<std::string> folders={};
  std::vector<std::vector<std::string>> files{};
  std::vector<double> sf = {5, 2};
};

std::vector<std::string> getFilesFromFolder(std::string const& folder
                                            ,std::string const& selector
                                            ,int maxFiles=624
                                            ,bool debug=true
                                            ) {
  StringVec files;
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(folder.c_str())) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir(dir)) != NULL and (maxFiles==-1 or files.size() < size_t(maxFiles))) {
      std::string filename(ent->d_name);
      if (filename.find(selector)!=std::string::npos ) {
        if (debug) { std::cout << "Found file: " << filename  << std::endl; }
        files.push_back(folder+filename);
      }
    }//for all entries
    closedir (dir);
  } else {
    /* could not open directory */
    perror ("");
    std::cerr << "Could not open folder " << folder << std::endl;
    throw EXIT_FAILURE;
  }
  return files;
}

void drawOccupancy(const Parameters& par, BeamCalGeo const& bcg,
                   std::string const& name, BCPadEnergies const& bcp) {

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


void readBackgroundFile(std::string const& backgroundFile, VD& countLeft, VD& countRight,
                        VVD& trainLeft, VVD& trainRight, int bxPerTrain,
                        double threshold, int maxEPad, int& nBX, double& totalEnergyBX, double& occupancyBX, double& maxEnergy){
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
        trainLeft[(nBX-1) / bxPerTrain][nPad] = 1;
	occupancyBX += 1;
      }
      if((*depRight)[nPad] > threshold) {
        countRight[nPad] += 1;
        trainRight[(nBX-1) / bxPerTrain][nPad] = 1;
	occupancyBX += 1;
      }
    }
  }
  occupancyBX /= 2.0;
  file->Close();
}


#endif // OCCUPANCYUTILITIES_HH
