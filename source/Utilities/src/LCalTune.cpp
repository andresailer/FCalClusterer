/* To run reconstruction over simulated files:

 parallel --files --joblog Rec.log -j3 --results Rec  RunReco {1} {2} \
 ::: 10 20 30 50 150 200 250 300 400 500 750 1000 1250 1500 \
 ::: `ls --color=never /eos/experiment/clicdp/grid/ilc/user/s/sailer/LumiCalSim/LC3/100/` &

 */

#include "RootUtils.hh"

#include "TError.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <sys/stat.h>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>

class IgnoreError {
  public:
  IgnoreError(): originalErrorLevel(gErrorIgnoreLevel)
  {
    gErrorIgnoreLevel=kError;
  }
  ~IgnoreError() {
    gErrorIgnoreLevel=originalErrorLevel;
  }
private:
  int originalErrorLevel;
};



struct Results {

  double energy;
  double mradAngle;
  std::map< std::string, TH1D*> hs;
  std::map< std::string, TH2D*> hs2D;

  Results( double e, double a ): energy(e), mradAngle(a), hs(), hs2D() {}

};

typedef std::vector<double> VecD;
typedef std::map<int, TEfficiency*> Efficiencies;
TEfficiency* angleEff( int energy ){
  int nbins(200);
  double minAngle(0.03), maxAngle(0.120);
  return new TEfficiency(Form("eff%d",energy),Form("%d GeV;#theta [rad];Efficiency",energy),
			 nbins,
			 minAngle,
			 maxAngle );
}


const double legendEntryHeight(0.0625);

void drawEfficiencies( RootUtils::PDFFile&, Efficiencies& efficiencies);
void drawEnergyResolutionFunction( RootUtils::PDFFile&, VecD& selectedEnergyEnergy, VecD& selectedSigmaEoverE, TF1&);


void saveCanvas(RootUtils::PDFFile& pdf, TCanvas& canv, int energy, TString title, int nHistos, double legStartX=0.65) {
  IgnoreError IE;
  const double startY = 0.9;
  auto* leg = RootUtils::BuildLegend(canv, legStartX, startY, legStartX+0.3, startY-(nHistos+1)*legendEntryHeight);
  RootUtils::AddLegendHeader( leg, {Form("%d GeV", energy )} );
  RootUtils::AutoSetYRange( canv, 1.1);

  pdf.AddCanvas( canv, title );

}

const int ANGLE_CONVERSION=100000;

Results* makeResults( double energy,
		      double mradAngle,
		      TString filename,
		      TEfficiency* eff) {

  auto results = new Results( energy, mradAngle );

  //Energy ResolutionHistogram
  {
    TString histoName( Form("e%da%d", int(energy), int(mradAngle*ANGLE_CONVERSION)) );
    TString histoTitle( Form( "%3.1f mrad; E_{Reco}[GeV]; N ", mradAngle*1000));
    double minE(0.0), maxE(energy*1.2);
    int nBins = 300;
    // if( energy > 330 ) {
    //   minE = 0.85 * energy;
    //   maxE = 1.15 * energy;
    // }
    results->hs["eRes"]=new TH1D( histoName, histoTitle, nBins, minE, maxE );
  }

  {
    TString histoName( Form("relE%da%d", int(energy), int(mradAngle*ANGLE_CONVERSION)) );
    TString histoTitle( Form( "%3.1f mrad; (E_{Reco} - E_{MC})/E_{MC}; N ", mradAngle*1000));
    double maxE = 0.35/sqrt(energy/20.0);
    double minE = -maxE;
    int nBins = 100;
    results->hs["eResRel"]=new TH1D( histoName, histoTitle, nBins, minE, maxE );
  }

  //Theta ResolutionHistogram
  {
    TString histoName( Form("t%da%d", int(energy), int(mradAngle*ANGLE_CONVERSION)) );
    TString histoTitle( Form( "%3.1f mrad; #theta_{Reco}[rad]; N ", mradAngle*1000));
    double minT= -0.002;
    double maxT = -minT;
    int nBins = 200;
    results->hs["tRes"]=new TH1D( histoName, histoTitle, nBins, minT, maxT );
  }


  //Phi ResolutionHistogram
  {
    TString histoName( Form("p%da%d", int(energy), int(mradAngle*ANGLE_CONVERSION)) );
    TString histoTitle( Form( "%3.1f mrad; #phi_{Reco}[rad]; N ", mradAngle*1000));
    double minP= -0.5;
    double maxP = -minP;
    int nBins = 200;
    results->hs["pRes"]=new TH1D( histoName, histoTitle, nBins, minP, maxP );
  }

  TFile *file = TFile::Open( filename );
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

  TH1D* eHisto = results->hs["eRes"];
  TH1D* rHisto = results->hs["eResRel"];
  TH1D* tHisto = results->hs["tRes"];
  TH1D* pHisto = results->hs["pRes"];

  for (int i = 0; i < tree->GetEntries() ;++i) {
    if( i > 2000) break;
    tree->GetEntry(i);

    eHisto->Fill( eReco );
    rHisto->Fill( (eReco - eMC) / eMC );
    tHisto->Fill( thetaReco-thetaMC );
    pHisto->Fill( phiReco-phiMC );

    // wihthin ten percent
    eff->Fill( fabs((eReco - eMC) / eMC) < 0.1 , mradAngle);
  }

  file->Close();
  return results;
}


double fitGauss( double* x, double* p ) {
  return p[2]/(p[1]*sqrt(2*M_PI)) * std::exp( - (x[0]-p[0]) * (x[0]-p[0]) / ( 2.0 * ( p[1]*p[1] ) ) );
}

double fitResolution( double* x, double* p ) {
  return p[0] / sqrt(x[0]) + p[1];
}


int main (int, char**) {

  RootUtils::SetStyle();
  RootUtils::PDFFile pdf("LCalTune.pdf");
  const bool allHistograms = not true;
  std::map< int, std::map< int, Results*> > histograms;
  Efficiencies efficiencies;
  //  std::vector<int> energies = {10, 20, 30, 50, 100, 150, 500};
  //#include "LCalTune_SmallSet.h"
  //#include "LCalTune_FullSet.h"
  #include "LCalTune_SelectedSet.h"


  //#include "LCalTune_HESet.h"

  // energies = {10, 100, 1000};

  //  energies = {1500};

  //TF1 thetaGauss("thetaGauss", fitGauss, -0.03, 0.03, 3);
  //TF1 phiGauss("phiGauss", fitGauss, -0.3, 0.3, 3);
  TF1 energyGauss("energyGauss", fitGauss, -0.3, 0.3, 3);
  TF1 fitRes("fitRes", fitResolution, 0, 1550, 2);
  fitRes.SetParameter(0, 0.22);
  fitRes.SetParameter(1, 0.0);

  for (auto energy: energies){

    efficiencies.insert( std::make_pair( energy, angleEff(energy) )
			);


    std::cout <<  "Running over " << energy << " GeV files "   << std::endl;
    for (auto mradAngleStr: angles) {

      int iAngle = int(atof(mradAngleStr.c_str())*ANGLE_CONVERSION );
      //if( iAngle > 9000) break;
      TString filename = Form("FCA_%d_%s.root", energy, mradAngleStr.c_str());

      // std::cout << "Opening "
      //		<< std::setw(50) << filename
      //		<< std::setw(10) << energy
      //		<< std::setw(10) << iAngle
      //		<< std::endl;

      histograms[energy][iAngle] = makeResults( double(energy),  atof(mradAngleStr.c_str()), filename, efficiencies[energy]);

    }

  }

  // histograms[30] = makeEnergyResolution( 100.0, "FCAlAna_100.root", "e100");

  VecD selectedEnergyMean, selectedEnergyDiff, selectedEnergyVariance, selectedEnergyEnergy;
  VecD selectedSigmaEoverE{};
  // Draw Energy Related Histograms
  {
    for (auto const& pIntMapIntResults: histograms ) {
      const int energy = pIntMapIntResults.first;
      const double legStartX = ( energy < 310 ) ? 0.20 : 0.20;
      RootUtils::Colors colors;

      std::vector<double> recEnergies;
      std::vector<double> recMeans;
      std::vector<double> recMeanErr;
      std::vector<double> recVariances;
      std::vector<double> dAngles;
      std::vector<double> zeros( angles.size(), 0.0 );

      TCanvas canvE("cE","cE");

      bool first = true; TString opt("");
      int histCounter = 0;
      for (auto const& pIntResults: pIntMapIntResults.second ) {
	const int iAngle = pIntResults.first;
	Results* results = pIntResults.second;
	if( not results ) continue;

	TH1D* eHisto = results->hs["eResRel"];
	energyGauss.SetParameter(0,eHisto->GetMean());
	energyGauss.SetParameter(1,eHisto->GetRMS());
	energyGauss.SetParameter(2,eHisto->GetBinContent(eHisto->GetMaximumBin()));

	//const double fitRange= 0.50;
	//auto fitresultptr = eHisto->Fit(&energyGauss, "QNS", "", energy-fitRange*energy, energy+fitRange*energy);
	auto fitresultptr = eHisto->Fit(&energyGauss, "QNS", "");

	if( int(fitresultptr) == 0 ){

	  const double mean = energyGauss.GetParameter(0);
	  const double meanErr = fitresultptr.Get()->Error(0);
	  const double variance = energyGauss.GetParameter(1);
	  if( fabs( mean ) < 0.1 ) {

	    recMeans.push_back( mean );
	    recMeanErr.push_back( meanErr );
	    recVariances.push_back( variance );

	    dAngles.push_back( double(iAngle)/double(ANGLE_CONVERSION) );

	  } else {

	    // recMeans.push_back( eHisto->GetMean() );
	    // recVariances.push_back( eHisto->GetRMS() );

	    std::cout << "********************************************************************************"  << std::endl;
	    std::cout << "Fit result ptr status " << int(fitresultptr) << std::endl;
	    std::cout << std::setw(15) << energy << std::setw(15) << iAngle << std::endl;
	    std::cout << std::setw(15) << mean << std::setw(15) << meanErr << std::endl;
	    std::cout << "********************************************************************************"  << std::endl;

	  }

	} else { // we have a very bad fit result, mean is consistent with zero
	  std::cout << "Fit result ptr status "  << int(fitresultptr) << std::endl;
	  std::cout << std::setw(15) << energy
		    << std::setw(15) << iAngle
		    << std::endl;

	  // recMeans.push_back( eHisto->GetMean() );
	  // recVariances.push_back( eHisto->GetRMS() );

	}


	if (first) {
	  first=false;
	  opt="";
	  //eHisto->SetAxisRange(0.0, energy*2,"X");
	}
	else { opt="same"; }
	histCounter++;

	RootUtils::SetAllColors(eHisto, colors.GetColor() );
	canvE.cd();
	eHisto->Draw(opt);

	if( histCounter > 8 && allHistograms ) {
	  saveCanvas(pdf, canvE, energy, Form("Energy Resolution vs Angle: %d GeV", energy ), histCounter, legStartX);
	  histCounter=0;
	  first = true;
	  colors.SetColors();
	}

      }

      if(histCounter > 0 && allHistograms) {
	saveCanvas(pdf, canvE, energy, Form("Energy Resolution vs Angle: %d GeV", energy), histCounter, legStartX);
      }

      if(not dAngles.empty() ) {
	IgnoreError IE;
	TCanvas canvR("cE","cE");
	TGraphErrors resolutions(dAngles.size(), &dAngles[0], &recMeans[0], &zeros[0], &recMeanErr[0] );
	// TGraphErrors resolutions(dAngles.size(), &dAngles[0], &recVariances[0], &zeros[0], &zeros[0] );
	resolutions.Draw("AP");
	gPad->Update();
	resolutions.GetXaxis()->SetTitle( "#theta [rad]" );
	resolutions.GetYaxis()->SetTitle( "(E_{Rec}-E_{MC})/E_{MC}" );
	resolutions.GetYaxis()->SetRangeUser(-0.9e-2, 0.9e-2);
	resolutions.SetMarkerStyle(kOpenStar);
	auto fitresultptr = resolutions.Fit("pol0", "S", "", 0.062, 0.07 );

	if ( int(fitresultptr) == 0 ) {

	  const double mean = fitresultptr.Get()->Parameter(0);
	  const double meanErr = fitresultptr.Get()->Error(0);

	  selectedEnergyMean.push_back( energy * ( mean + 1.0 ) );
	  selectedEnergyDiff.push_back( mean * energy );
	  selectedEnergyVariance.push_back( energy * meanErr ); // variance
	  selectedEnergyEnergy.push_back( energy );

	}
	//canvR.SaveAs("test.eps");
	pdf.AddCanvas(canvR, Form("ERec vs. Angle: %d GeV", energy) );
      }

      if(not dAngles.empty() ) {
	IgnoreError IE;
	TCanvas canvR("cE","cE");
	TGraphErrors resolutions(dAngles.size(), &dAngles[0], &recVariances[0], &zeros[0], &zeros[0] );
	resolutions.Draw("AP");
	gPad->Update();
	resolutions.GetXaxis()->SetTitle( "#theta [rad]" );
	resolutions.GetYaxis()->SetTitle( "#sigma((E_{Rec}-E_{MC}))/E_{MC}" );
	resolutions.GetYaxis()->SetRangeUser(0.0, 0.5);
	resolutions.SetMarkerStyle(kOpenStar);
	auto fitresultptr = resolutions.Fit("pol0", "S", "", 0.062, 0.07 );
	if ( int(fitresultptr) == 0 ) {
	  const double mean = fitresultptr.Get()->Parameter(0);
	  selectedSigmaEoverE.push_back( mean );
	}

	//canvR.SaveAs("test.eps");
	pdf.AddCanvas(canvR, Form("SigmaE vs. Angle: %d GeV", energy) );

      }

    }//all histograms


  }// energy related histograms


  { //Energy Resolution linearity
      if(not selectedEnergyEnergy.empty() ) {
	TCanvas canvER("cER","cER");
	std::vector<double> zeros( selectedEnergyEnergy.size(), 0.0 );
	TGraphErrors resolutions(selectedEnergyEnergy.size(),
				  &selectedEnergyEnergy[0],
				  &selectedEnergyMean[0],
				  &zeros[0],
				  &selectedEnergyVariance[0] );
	resolutions.Draw("AP");
	gPad->Update();
	resolutions.GetXaxis()->SetTitle( "E_{MC} [GeV]" );
	resolutions.GetYaxis()->SetTitle( "E_{Rec} [GeV]" );
	resolutions.Fit("pol1");
	resolutions.SetTitle(" Resolution " );
	const double startY(0.9), legStartX(0.3);
	auto* leg = RootUtils::BuildLegend(canvER, legStartX, startY, legStartX+0.3, startY-(1+2)*legendEntryHeight);
	RootUtils::AddLegendHeader( leg, {"Reconstructed Energy", "#theta = 70 mrad"} );
	IgnoreError IE;
	pdf.AddCanvas(canvER, "Energy Resolution vs. Energy" );
	canvER.SaveAs("ERes.eps");

      }

      if(not selectedEnergyEnergy.empty() ) {
	TCanvas canvER("cER","cER");
	std::vector<double> zeros( selectedEnergyEnergy.size(), 0.0 );
	TGraphErrors resolutions(selectedEnergyEnergy.size(),
				  &selectedEnergyEnergy[0],
				  &selectedEnergyDiff[0],
				  &zeros[0],
				  &selectedEnergyVariance[0] );
	resolutions.Draw("APM");
	resolutions.SetMarkerStyle( kFullStar );
	gPad->Update();
	resolutions.GetXaxis()->SetTitle( "E_{MC} [GeV]" );
	resolutions.GetYaxis()->SetTitle( "E_{Rec} - E_{MC} [GeV]" );
	resolutions.GetYaxis()->SetRangeUser(-0.7, 0.7);
	resolutions.SetTitle(" Resolution " );
	const double startY(0.9), legStartX(0.3);
	auto* leg = RootUtils::BuildLegend(canvER, legStartX, startY, legStartX+0.3, startY-(1+2)*legendEntryHeight);
	RootUtils::AddLegendHeader( leg, {"Reconstructed Energy", "#theta = 70 mrad"} );
	IgnoreError IE;
	pdf.AddCanvas(canvER, "Energy Resolution vs. Energy" );
	canvER.SaveAs("EResDiff.eps");

      }


  }

  drawEnergyResolutionFunction(pdf, selectedEnergyEnergy, selectedSigmaEoverE, fitRes);


  if(allHistograms){
    for (auto const& pIntMapIntResults: histograms ) {
      const int energy = pIntMapIntResults.first;
      RootUtils::Colors colors;
      TCanvas canvT("cT","cT");
      bool first = true; TString opt("");
      int histCounter = 0;
      for (auto const& pIntResults: pIntMapIntResults.second ) {
	//const int iAngle = pIntResults.first;
	Results* results = pIntResults.second;
	if( not results ) continue;
	TH1D* tHisto = results->hs["eRes"];
	if (first) { first=false; opt=""; }
	else { opt="same"; }
	histCounter++;
	RootUtils::SetAllColors(tHisto, colors.GetColor() );
	canvT.cd();
	tHisto->Draw(opt);
	if( histCounter > 5 ) {
	  saveCanvas(pdf, canvT, energy, Form("Energy vs Angle: %d GeV", energy ), histCounter, 0.20);
	  histCounter=0;
	  first = true;
	  colors.SetColors();
	}
      }
      if (histCounter > 0) {
	saveCanvas(pdf, canvT, energy, Form("Energy vs Angle: %d GeV", energy), histCounter, 0.20);
      }
    }
  }//Theta related histograms



  // Draw Theta Related Histograms
  if(allHistograms){
    for (auto const& pIntMapIntResults: histograms ) {
      const int energy = pIntMapIntResults.first;
      RootUtils::Colors colors;
      TCanvas canvT("cT","cT");
      bool first = true; TString opt("");
      int histCounter = 0;
      for (auto const& pIntResults: pIntMapIntResults.second ) {
	//const int iAngle = pIntResults.first;
	Results* results = pIntResults.second;
	if( not results ) continue;
	TH1D* tHisto = results->hs["tRes"];
	if (first) { first=false; opt=""; }
	else { opt="same"; }
	histCounter++;
	RootUtils::SetAllColors(tHisto, colors.GetColor() );
	canvT.cd();
	tHisto->Draw(opt);
	if( histCounter > 5 ) {
	  saveCanvas(pdf, canvT, energy, Form("Theta Resolution vs Angle: %d GeV", energy ), histCounter);
	  histCounter=0;
	  first = true;
	  colors.SetColors();
	}
      }
      if (histCounter > 0) {
	saveCanvas(pdf, canvT, energy, Form("Theta Resolution vs Angle: %d GeV", energy), histCounter);
      }
    }
  }//Theta related histograms


  // Draw Phi Related Histograms
  if(allHistograms){
    for (auto const& pIntMapIntResults: histograms ) {
      RootUtils::Colors colors;
      const int energy = pIntMapIntResults.first;
      TCanvas canv("cP","cP");
      bool first = true; TString opt("");
      int histCounter = 0;
      for (auto const& pIntResults: pIntMapIntResults.second ) {
	//const int iAngle = pIntResults.first;
	Results* results = pIntResults.second;
	if( not results ) continue;
	TH1D* tHisto = results->hs["pRes"];
	if (first) { first=false; opt=""; }
	else { opt="same"; }
	histCounter++;
	RootUtils::SetAllColors(tHisto, colors.GetColor() );
	canv.cd();
	tHisto->Draw(opt);
	if( histCounter > 5 ) {
	  saveCanvas(pdf, canv, energy, Form("Phi Resolution vs Angle: %d GeV", energy ), histCounter);
	  histCounter=0;
	  first = true;
	  colors.SetColors();
	}
      }
      if (histCounter > 0) {
	saveCanvas(pdf, canv, energy, Form("Phi Resolution vs Angle: %d GeV", energy), histCounter);
      }
    }
  }//Phi related histograms

  drawEfficiencies( pdf, efficiencies );

  /// Finish the pdf file with an empty slide to keep the last title
  TCanvas empty("e","e");
  pdf.AddCanvas( empty, "end");

  return 0;
}

void drawEfficiencies( RootUtils::PDFFile& pdf, Efficiencies& efficiencies) {

  for( auto& pIntEfficiency: efficiencies) {

      const int energy = pIntEfficiency.first;
      TCanvas canv("cE","cE");
      pIntEfficiency.second->Draw();

      IgnoreError IE;
      const double startY = 0.4;
      const double legStartX = 0.2;
      RootUtils::BuildLegend(canv, legStartX, startY, legStartX+0.3, startY-(1)*legendEntryHeight);
      //RootUtils::AddLegendHeader( leg, {Form("%d GeV", energy )} );
      TString title = Form("Reconstruction Efficiency %d GeV", energy );
      pdf.AddCanvas( canv, title );
      TString filename=Form("RecoEff%d.eps", energy);
      canv.SaveAs(filename);
    }
  return;
}//end drawEfficiencies


/// draw the sigmaE/E function prop to a/sqrt(E)
void drawEnergyResolutionFunction( RootUtils::PDFFile& pdf,
				   VecD& selectedEnergyEnergy,
				   VecD& selectedSigmaEoverE,
				   TF1& fitRes) {
  if(not selectedEnergyEnergy.empty() ) {

    TCanvas canvER("cER","cER");
    // std::vector<double> zeros( selectedEnergyEnergy.size(), 0.0 );
    TGraph resolutions(selectedEnergyEnergy.size(),
		       &selectedEnergyEnergy[0],
		       &selectedSigmaEoverE[0] );
    resolutions.Draw("AP");
    gPad->Update();
    resolutions.GetXaxis()->SetTitle( "E_{MC} [GeV]" );
    resolutions.GetYaxis()->SetTitle( "#sigma_{E}/E_{MC}" );
    resolutions.Fit(&fitRes);
    resolutions.SetMarkerStyle(kOpenStar);
    resolutions.SetTitle(" Resolution " );
    const double startY(0.9), legStartX(0.5);
    auto* leg = RootUtils::BuildLegend(canvER, legStartX, startY, legStartX+0.3, startY-(1+2)*legendEntryHeight);
    RootUtils::AddLegendHeader( leg, {"Energy Resolution", "a/$sqrt{E}+b"} );
    IgnoreError IE;
    pdf.AddCanvas(canvER, "SigmaE/E" );
    canvER.SaveAs("EResFunc.eps");

  }

}// end drawEnergyResolutionFunction
