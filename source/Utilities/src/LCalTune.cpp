
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>

#include <iostream>
#include <vector>
#include <map>
#include <string>

TH1D* makeEnergyResolution( double energy,
			    std::string filename,
			    std::string histoName,
			    std::string histoTitle ) {

  TFile *file = TFile::Open( filename.c_str() );
  if(not file) return NULL;
  TTree *tree;
  file->GetObject( "fcalAna", tree );
  if( not tree) return NULL;

  double thetaMC, thetaReco, eMC, eReco, phiMC, phiReco;
  tree->SetBranchAddress( "thetaMC"   , &thetaMC  );
  tree->SetBranchAddress( "thetaReco" , &thetaReco);
  tree->SetBranchAddress( "eMC"       , &eMC      );
  tree->SetBranchAddress( "eReco"     , &eReco    );
  tree->SetBranchAddress( "phiMC"     , &phiMC    );
  tree->SetBranchAddress( "phiReco"   , &phiReco  );

  double minE= 0.0;
  double maxE = 550;
  int nBins = 2*maxE;

  TH1D* eHisto = new TH1D( histoName.c_str(), histoTitle.c_str(), nBins, minE, maxE );

  std::cout << "Running over" << tree->GetEntries() << " events."   << std::endl;

  for (int i = 0; i < tree->GetEntries() ;++i) {
    tree->GetEntry(i);

    eHisto->Fill( eReco );

  }

  return eHisto;  

}



int main (int argc, char **args) {
  

  std::map< int, TH1D*> histograms;
  std::string axes = "Reconstructed Energy; E_{Reco}[GeV]; N";

  // histograms[10] = makeEnergyResolution(   10.0, "FCalAna_10.root", "e10", axes);
  // histograms[20] = makeEnergyResolution(   20.0, "FCalAna_20.root", "e20", axes);
  // histograms[30] = makeEnergyResolution(   30.0, "FCalAna_30.root", "e30", axes);
  // histograms[100] = makeEnergyResolution( 100.0, "FCalAna_100.root", "e100", axes);
  // histograms[150] = makeEnergyResolution( 150.0, "FCalAna_150.root", "e150", axes);
  // histograms[200] = makeEnergyResolution( 200.0, "FCalAna_200.root", "e200", axes);
  // histograms[250] = makeEnergyResolution( 250.0, "FCalAna_250.root", "e250", axes);

  // histograms[10] = makeEnergyResolution(   10.0, "FCalAna_55_10.root", "e10", axes);
  // histograms[20] = makeEnergyResolution(   20.0, "FCalAna_55_20.root", "e20", axes);
  // histograms[30] = makeEnergyResolution(   30.0, "FCalAna_55_30.root", "e30", axes);
  // histograms[100] = makeEnergyResolution( 100.0, "FCalAna_55_100.root", "e100", axes);
  // histograms[150] = makeEnergyResolution( 150.0, "FCalAna_55_150.root", "e150", axes);
  // histograms[200] = makeEnergyResolution( 200.0, "FCalAna_55_200.root", "e200", axes);
  // histograms[250] = makeEnergyResolution( 250.0, "FCalAna_55_250.root", "e250", axes);


  histograms[10] = makeEnergyResolution(   10.0, "FCalAna_55_2_10.root", "e10", axes);
  histograms[20] = makeEnergyResolution(   20.0, "FCalAna_55_2_20.root", "e20", axes);
  histograms[30] = makeEnergyResolution(   30.0, "FCalAna_55_2_30.root", "e30", axes);
  histograms[100] = makeEnergyResolution( 100.0, "FCalAna_55_2_100.root", "e100", axes);
  histograms[150] = makeEnergyResolution( 150.0, "FCalAna_55_2_150.root", "e150", axes);
  histograms[200] = makeEnergyResolution( 200.0, "FCalAna_55_2_200.root", "e200", axes);
  histograms[250] = makeEnergyResolution( 250.0, "FCalAna_55_2_250.root", "e250", axes);



  // histograms[30] = makeEnergyResolution( 100.0, "FCAlAna_100.root", "e100", axes);


  std::vector<double> energies;
  std::vector<double> means;


  TCanvas c("c","c", 800, 700);

  bool first = true;
  for (std::map< int, TH1D*>::iterator it = histograms.begin(); it != histograms.end(); ++it) {
    if ( not it->second ) continue;

    energies.push_back( it->first );
    //means.push_back( it->second->GetMean() );
    means.push_back( double(it->second->GetMaximumBin())*0.5 );
    std::cout << energies.back() << "   " << means.back() << "   " << it->second->Integral()  << std::endl;

    

    if(first) {
      it->second->Draw();
      first = false;
      std::cout << "drawing first"  << std::endl;
    } else {
      it->second->Draw("same");
      std::cout << "Drawing same"  << std::endl;
    }

  }
  

  c.SaveAs( "Energies.eps");

  TGraph resolutions(energies.size(), &energies[0], &means[0] );
  
  resolutions.Draw("APL");
  resolutions.SetMarkerStyle(kOpenStar);
  gPad->Update();
  resolutions.Fit("pol1");

  c.SaveAs("Resolutions.eps");

}
