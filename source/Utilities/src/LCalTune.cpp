/* To run reconstruction over simulated files:

   parallel --files --joblog Rec.log -j3 --results Rec  RunReco {1} {2} \
   ::: 10 20 30 50 150 200 250 300 400 500 750 1000 1250 1500 \
   ::: `ls --color=never /eos/experiment/clicdp/grid/ilc/user/s/sailer/LumiCalSim/LC3/100/` &

*/

#include "RootUtils.hh"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include "TError.h"

#include <sys/stat.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

class IgnoreError {
public:
  IgnoreError() : originalErrorLevel(gErrorIgnoreLevel) { gErrorIgnoreLevel = kError; }
  ~IgnoreError() { gErrorIgnoreLevel = originalErrorLevel; }

private:
  int originalErrorLevel;
};

struct Results {
  double energy;
  double mradAngle;
  std::map<std::string, TH1D*> hs;
  std::map<std::string, TH2D*> hs2D;

  Results(double e, double a) : energy(e), mradAngle(a), hs(), hs2D() {}
};

typedef std::vector<double> VecD;
typedef std::map<int, TEfficiency*> Efficiencies;
TEfficiency* angleEff(int energy) {
  int    nbins(200);
  double minAngle(0.03), maxAngle(0.120);
  return new TEfficiency(Form("eff%d", energy), Form("%d GeV;#theta [mrad];Efficiency", energy), nbins, minAngle, maxAngle);
}

const double legendEntryHeight(0.0625);

void drawEfficiencies(RootUtils::PDFFile&, Efficiencies& efficiencies);
void drawEnergyResolutionFunction(RootUtils::PDFFile&, VecD& selectedEnergyEnergy, VecD& selectedSigmaEoverE,
                                  VecD& selectedSigmaEoverEErr, TF1&);

constexpr size_t THETALOG_SIZE  = 100;
constexpr size_t THETALOG_START = 30;

void saveCanvas(RootUtils::PDFFile& pdf, TCanvas& canv, int energy, TString title, int nHistos, double legStartX = 0.65) {
  IgnoreError  IE;
  const double startY = 0.9;
  auto* leg = RootUtils::BuildLegend(canv, legStartX, startY, legStartX + 0.3, startY - (nHistos + 1) * legendEntryHeight);
  RootUtils::AddLegendHeader(leg, {Form("%d GeV", energy)});
  RootUtils::AutoSetYRange(canv, 1.1);

  pdf.AddCanvas(canv, title);
}

const int ANGLE_CONVERSION   = 100000;
const int TO_MRAD_CONVERSION = 100;

const double fiducialMinimum = 0.062;
const double fiducialMaximum = 0.077;

Results* makeResults(double energy, double mradAngle, TString filename, TEfficiency* eff) {
  auto results = new Results(energy, mradAngle);

  //Energy ResolutionHistogram
  {
    TString histoName(Form("e%da%d", int(energy), int(mradAngle * ANGLE_CONVERSION)));
    TString histoTitle(Form("%3.1f mrad; E_{Reco}[GeV]; N ", mradAngle * 1000));
    double  minE(0.0), maxE(energy * 1.2);
    int     nBins = 300;
    // if( energy > 330 ) {
    //   minE = 0.85 * energy;
    //   maxE = 1.15 * energy;
    // }
    results->hs["eRes"] = new TH1D(histoName, histoTitle, nBins, minE, maxE);
  }

  {
    TString histoName(Form("relE%da%d", int(energy), int(mradAngle * ANGLE_CONVERSION)));
    TString histoTitle(Form("%3.1f mrad; (E_{Reco} - E_{MC})/E_{MC}; N ", mradAngle * 1000));
    double  maxE           = 0.35 / sqrt(energy / 20.0);
    double  minE           = -maxE;
    int     nBins          = 100;
    results->hs["eResRel"] = new TH1D(histoName, histoTitle, nBins, minE, maxE);
  }

  //Theta ResolutionHistograms
  {
    TString histoName(Form("t%da%d", int(energy), int(mradAngle * ANGLE_CONVERSION)));
    TString histoNameN(Form("tN%da%d", int(energy), int(mradAngle * ANGLE_CONVERSION)));
    TString histoNameE(Form("tE%da%d", int(energy), int(mradAngle * ANGLE_CONVERSION)));
    TString histoTitle(Form("%3.1f mrad; (#theta_{Reco}-#theta_{MC})[rad]; N ", mradAngle * 1000));
    TString histoTitleN(Form("N %3.1f mrad; (#theta_{Reco}-#theta_{MC})[rad]; N ", mradAngle * 1000));
    TString histoTitleE(Form("E %3.1f mrad; (#theta_{Reco}-#theta_{MC})[rad]; N ", mradAngle * 1000));
    double  minT         = -0.0004;
    double  maxT         = -minT;
    int     nBins        = 200;
    results->hs["tRes"]  = new TH1D(histoName, histoTitle, nBins, minT, maxT);
    results->hs["tResN"] = new TH1D(histoNameN, histoTitleN, nBins, minT, maxT);
    results->hs["tResE"] = new TH1D(histoNameE, histoTitleE, nBins, minT, maxT);
    for (int i = 0; i < 100; ++i) {
      TString histoTitleL(Form("L %3.1f mrad %i; (#theta_{Reco}-#theta_{MC})[rad]; N ", mradAngle * 1000, i));
      TString histoNameL(Form("tL%da%d%i", int(energy), int(mradAngle * ANGLE_CONVERSION), i));
      results->hs[Form("tResL%i", i)] = new TH1D(histoNameL, histoTitleL, nBins, minT, maxT);
    }
  }

  //Phi ResolutionHistogram
  {
    TString histoName(Form("p%da%d", int(energy), int(mradAngle * ANGLE_CONVERSION)));
    TString histoTitle(Form("%3.1f mrad; (#phi_{Reco}-#phi_{MC})[rad]; N ", mradAngle * 1000));
    double  minP        = -0.5;
    double  maxP        = -minP;
    int     nBins       = 200;
    results->hs["pRes"] = new TH1D(histoName, histoTitle, nBins, minP, maxP);
  }

  TFile* file = TFile::Open(filename);
  if (not file)
    return NULL;
  TTree* tree;
  file->GetObject("fcalAna", tree);
  if (not tree)
    return NULL;

  double               thetaMC, thetaReco, eMC, eReco, phiMC, phiReco, thetaN, thetaE;
  std::vector<double>* thetaLog = nullptr;
  tree->SetBranchAddress("thetaMC", &thetaMC);
  tree->SetBranchAddress("thetaReco", &thetaReco);
  tree->SetBranchAddress("thetaN", &thetaN);
  tree->SetBranchAddress("thetaE", &thetaE);
  tree->SetBranchAddress("thetaLog", &thetaLog);
  tree->SetBranchAddress("eMC", &eMC);
  tree->SetBranchAddress("eReco", &eReco);
  tree->SetBranchAddress("phiMC", &phiMC);
  tree->SetBranchAddress("phiReco", &phiReco);

  TH1D* eHisto  = results->hs["eRes"];
  TH1D* rHisto  = results->hs["eResRel"];
  TH1D* tHisto  = results->hs["tRes"];
  TH1D* tHistoN = results->hs["tResN"];
  TH1D* tHistoE = results->hs["tResE"];
  //auto tHistoL = results->hs["tResL"];
  TH1D* pHisto = results->hs["pRes"];

  std::vector<TH1D*> thetaLogHisto;
  for (size_t j = 0; j < THETALOG_SIZE; ++j) {
    thetaLogHisto.push_back(results->hs[Form("tResL%i", int(j))]);
  }

  for (int i = 0; i < tree->GetEntries(); ++i) {
    //if( i > 2000) break;
    tree->GetEntry(i);

    eHisto->Fill(eReco);
    rHisto->Fill((eReco - eMC) / eMC);
    tHisto->Fill(thetaReco - thetaMC);
    //tHisto->Fill( thetaLog-thetaMC );
    tHistoN->Fill(thetaN - thetaMC);
    tHistoE->Fill(thetaE - thetaMC);
    if (thetaLog) {
      for (size_t j = 0; j < thetaLog->size(); ++j) {
        thetaLogHisto[j]->Fill(thetaLog->at(j) - thetaMC);
      }
    }
    pHisto->Fill(phiReco - phiMC);

    // wihthin ten percent
    eff->Fill(fabs((eReco - eMC) / eMC) < 0.1, mradAngle);
  }

  file->Close();
  return results;
}

double fitGauss(double* x, double* p) {
  return p[2] / (p[1] * sqrt(2 * M_PI)) * std::exp(-(x[0] - p[0]) * (x[0] - p[0]) / (2.0 * (p[1] * p[1])));
}

double fitResolution(double* x, double* p) {
  return p[0] / sqrt(x[0]);  // + p[1];
}

int main(int, char**) {
  RootUtils::SetStyle();
  RootUtils::PDFFile pdf("LCalTune_1500.pdf");
  const bool         allHistograms = true;
  std::map<int, std::map<int, Results*>> histograms;
  Efficiencies efficiencies;
//  std::vector<int> energies = {10, 20, 30, 50, 100, 150, 500};
//#include "LCalTune_SmallSet.h"
//#include "LCalTune_FullSet.h"
//#include "LCalTune_Compare.h"

  std::vector<int> energies = {1500};

  std::vector<std::string> angles = {"0.0625", "0.0626", "0.0627", "0.0628", "0.0629", "0.063",  "0.0631",
                                     "0.0632", "0.0633", "0.0634", "0.0635", "0.0636", "0.0637", "0.0638",
                                     "0.0639", "0.064",  "0.0641", "0.0642", "0.0643", "0.0644", "0.0645",
                                     "0.0646", "0.0647", "0.0648", "0.0649", "0.065"};

  //#include "LCalTune_HESet.h"

  // energies = {10, 100, 1000};

  //energies = {1500};

  //TF1 thetaGauss("thetaGauss", fitGauss, -0.03, 0.03, 3);
  //TF1 phiGauss("phiGauss", fitGauss, -0.3, 0.3, 3);
  TF1 energyGauss("energyGauss", fitGauss, -0.3, 0.3, 3);
  TF1 fitRes("fitRes", fitResolution, 0, 1550, 1);
  fitRes.SetParameter(0, 0.22);
  //  fitRes.SetParameter(1, 0.0);
  fitRes.SetParName(0, "a");
  //  fitRes.SetParName(1, "b");

  for (auto energy : energies) {
    efficiencies.insert(std::make_pair(energy, angleEff(energy)));

    std::cout << "Running over " << energy << " GeV files " << std::endl;
    for (auto mradAngleStr : angles) {
      int iAngle = int(atof(mradAngleStr.c_str()) * ANGLE_CONVERSION);
      //if( iAngle > 9000) break;
      TString filename = Form("FCA_%d_%s.root", energy, mradAngleStr.c_str());

      // std::cout << "Opening "
      //                << std::setw(50) << filename
      //                << std::setw(10) << energy
      //                << std::setw(10) << iAngle
      //                << std::endl;

      histograms[energy][iAngle] = makeResults(double(energy), atof(mradAngleStr.c_str()), filename, efficiencies[energy]);
    }
  }

  // histograms[30] = makeEnergyResolution( 100.0, "FCAlAna_100.root", "e100");

  VecD selectedEnergyMean, selectedEnergyRelDiff, selectedEnergyVariance, selectedEnergyEnergy;
  VecD selectedEnergyRelVar;
  VecD selectedSigmaEoverE{}, selectedSigmaEoverEErr{};
  // Draw Energy Related Histograms
  {
    for (auto const& pIntMapIntResults : histograms) {
      const int         energy    = pIntMapIntResults.first;
      const double      legStartX = (energy < 310) ? 0.20 : 0.20;
      RootUtils::Colors colors;

      std::vector<double> recEnergies;
      std::vector<double> recMeans;
      std::vector<double> recMeanErr;
      std::vector<double> recVariances;
      std::vector<double> recVariancesErr;
      std::vector<double> dAngles;
      VecD                tMeans, tRMS, pMeans, pRMS, tAngles, pAngles;
      VecD                tMeansN, tRMSN;
      VecD                tMeansE, tRMSE;
      std::vector<VecD>   tMeansL(THETALOG_SIZE), tRMSL(THETALOG_SIZE);
      TCanvas             canvE("cE", "cE");

      bool    first = true;
      TString opt("");
      int     histCounter = 0;
      for (auto const& pIntResults : pIntMapIntResults.second) {
        const int iAngle  = pIntResults.first;
        Results*  results = pIntResults.second;
        if (not results)
          continue;

        TH1D* eHisto = results->hs["eResRel"];
        energyGauss.SetParameter(0, eHisto->GetMean());
        energyGauss.SetParameter(1, eHisto->GetRMS());
        energyGauss.SetParameter(2, eHisto->GetBinContent(eHisto->GetMaximumBin()));

        //const double fitRange= 0.50;
        //auto fitresultptr = eHisto->Fit(&energyGauss, "QNS", "", energy-fitRange*energy, energy+fitRange*energy);
        auto fitresultptr = eHisto->Fit(&energyGauss, "QNS", "");

        if (int(fitresultptr) == 0) {
          const double mean     = energyGauss.GetParameter(0);
          const double meanErr  = fitresultptr.Get()->Error(0);
          const double variance = energyGauss.GetParameter(1);
          if (fabs(mean) < 0.1 && meanErr < 0.005) {
            recMeans.push_back(mean);
            recMeanErr.push_back(meanErr);
            recVariances.push_back(variance);
            recVariancesErr.push_back(fitresultptr.Get()->Error(1));

            dAngles.push_back(double(iAngle) / double(ANGLE_CONVERSION));

          } else {
            // recMeans.push_back( eHisto->GetMean() );
            // recVariances.push_back( eHisto->GetRMS() );

            std::cout << "********************************************************************************" << std::endl;
            std::cout << "Fit result ptr status " << int(fitresultptr) << std::endl;
            std::cout << std::setw(15) << energy << std::setw(15) << iAngle << std::endl;
            std::cout << std::setw(15) << mean << std::setw(15) << meanErr << std::endl;
            std::cout << "********************************************************************************" << std::endl;
          }

        } else {  // we have a very bad fit result, mean is consistent with zero
          std::cout << "Fit result ptr status " << int(fitresultptr) << std::endl;
          std::cout << std::setw(15) << energy << std::setw(15) << iAngle << std::endl;

          // recMeans.push_back( eHisto->GetMean() );
          // recVariances.push_back( eHisto->GetRMS() );
        }

        ///********************************************************************************
        ///********************************************************************************
        ///********************************************************************************

        {  /// theta resolutions
          TH1D*  tHisto = results->hs["tRes"];
          double mean(tHisto->GetMean());
          double rms(tHisto->GetRMS());
          tMeans.push_back(mean * 1000);
          tRMS.push_back(rms * 1000);
          tAngles.push_back(double(iAngle) / double(TO_MRAD_CONVERSION));
        }  // angular resolutions

        {  /// theta resolutions
          TH1D*  tHisto = results->hs["tResN"];
          double mean(tHisto->GetMean());
          double rms(tHisto->GetRMS());
          tMeansN.push_back(mean * 1000);
          tRMSN.push_back(rms * 1000);
        }  // angular resolutions

        {  /// theta resolutions
          TH1D*  tHisto = results->hs["tResE"];
          double mean(tHisto->GetMean());
          double rms(tHisto->GetRMS());
          tMeansE.push_back(mean * 1000);
          tRMSE.push_back(rms * 1000);
        }  // angular resolutions

        {  /// theta resolutions
          for (size_t i = THETALOG_START; i < THETALOG_SIZE; ++i) {
            TH1D*  tHisto = results->hs[Form("tResL%i", int(i))];
            double mean(tHisto->GetMean());
            double rms(tHisto->GetRMS());
            tMeansL[i].push_back(mean * 1000);
            tRMSL[i].push_back(rms * 1000);
          }
        }  // angular resolutions

        {  /// phi resolutions
          TH1D*  pHisto = results->hs["pRes"];
          double mean(pHisto->GetMean());
          double rms(pHisto->GetRMS());
          pMeans.push_back(mean * TMath::RadToDeg());
          pRMS.push_back(rms * TMath::RadToDeg());
          pAngles.push_back(double(iAngle) / double(ANGLE_CONVERSION));
        }  // angular resolutions

        if (allHistograms) {
          if (first) {
            first = false;
            opt   = "";
            //eHisto->SetAxisRange(0.0, energy*2,"X");
          } else {
            opt = "same";
          }
          histCounter++;

          RootUtils::SetAllColors(eHisto, colors.GetColor());
          canvE.cd();
          eHisto->Draw(opt);

          if (histCounter > 8 && allHistograms) {
            saveCanvas(pdf, canvE, energy, Form("Energy Resolution vs Angle: %d GeV", energy), histCounter, legStartX);
            histCounter = 0;
            first       = true;
            colors.SetColors();
          }

        }  //all histograms
      }

      if (histCounter > 0 && allHistograms) {
        saveCanvas(pdf, canvE, energy, Form("Energy Resolution vs Angle: %d GeV", energy), histCounter, legStartX);
      }

      if (not dAngles.empty()) {
        IgnoreError         IE;
        TCanvas             canvR("cE", "cE");
        std::vector<double> zeros(angles.size(), 0.0);
        TGraphErrors        resolutions(dAngles.size(), &dAngles[0], &recMeans[0], &zeros[0], &recMeanErr[0]);
        resolutions.Draw("AP");
        gPad->Update();
        resolutions.SetTitle(Form("%d GeV", energy));
        resolutions.GetXaxis()->SetTitle("#theta [mrad]");
        resolutions.GetYaxis()->SetTitle("(E_{Rec}-E_{MC})/E_{MC}  ");
        resolutions.GetYaxis()->SetTitleOffset(1.31);
        resolutions.GetYaxis()->SetRangeUser(-0.01, 0.005);
        resolutions.SetMarkerStyle(kOpenStar);
        RootUtils::SetAllColors(resolutions);

        const double legStartX_G(0.2);
        const double startY(0.9);
        auto* leg = RootUtils::BuildLegend(canvR, legStartX_G, startY, legStartX_G + 0.3, startY - (2) * legendEntryHeight);
        RootUtils::AddLegendHeader(leg, {"Relative Reconstructed Energy"});

        std::cout << "Fitting energy diff: " << energy << std::endl;
        auto fitresultptr = resolutions.Fit("pol0", "S", "", fiducialMinimum, fiducialMaximum);

        if (int(fitresultptr) == 0) {
          const double mean    = fitresultptr.Get()->Parameter(0);
          const double meanErr = fitresultptr.Get()->Error(0);

          selectedEnergyMean.push_back(energy * (mean + 1.0));
          selectedEnergyVariance.push_back(energy * meanErr);  // variance
          selectedEnergyEnergy.push_back(energy);
          selectedEnergyRelDiff.push_back(mean);
          selectedEnergyRelVar.push_back(meanErr);  // variance
        }
        canvR.SaveAs(Form("ERecVsAngle_%d.eps", energy));
        pdf.AddCanvas(canvR, Form("ERec vs. Angle: %d GeV", energy));
      }

      if (not dAngles.empty()) {
        IgnoreError         IE;
        TCanvas             canvR("cE", "cE");
        std::vector<double> zeros(angles.size(), 0.0);
        TGraphErrors        resolutions(dAngles.size(), &dAngles[0], &recVariances[0], &zeros[0], &recVariancesErr[0]);
        resolutions.Draw("AP");
        gPad->Update();
        std::cout << "Fitting sigmaE energy: " << energy << std::endl;
        resolutions.SetTitle(Form("%d GeV", energy));
        resolutions.GetXaxis()->SetTitle("#theta [mrad]");
        resolutions.GetYaxis()->SetTitle("#sigma(E_{Rec})/E_{MC}             ");
        resolutions.GetYaxis()->SetTitleOffset(1.31);
        resolutions.GetYaxis()->SetRangeUser(0.0, 0.02);
        resolutions.SetMarkerStyle(kOpenStar);
        RootUtils::SetAllColors(resolutions);

        const double legStartX_G(0.2);
        const double startY(0.9);
        auto* leg = RootUtils::BuildLegend(canvR, legStartX_G, startY, legStartX_G + 0.3, startY - (2) * legendEntryHeight);
        RootUtils::AddLegendHeader(leg, {"Energy Resolution"});

        auto fitresultptr = resolutions.Fit("pol0", "S", "", fiducialMinimum, fiducialMaximum);
        if (int(fitresultptr) == 0) {
          const double mean = fitresultptr.Get()->Parameter(0);
          const double err  = fitresultptr.Get()->Error(0);
          selectedSigmaEoverE.push_back(mean);
          selectedSigmaEoverEErr.push_back(err);
        }

        //canvR.SaveAs("test.eps");
        canvR.SaveAs(Form("SigmaEVsAngle_%d.eps", energy));
        pdf.AddCanvas(canvR, Form("SigmaE vs. Angle: %d GeV", energy));
      }

      /// theta vs angle
      if (not tAngles.empty()) {
        IgnoreError         IE;
        TCanvas             canvR("cE", "cE");
        std::vector<double> zeros(angles.size(), 0.0);
        TGraphErrors        bias(tAngles.size(), &tAngles[0], &tMeans[0], &zeros[0], &zeros[0]);
        bias.Draw("AP");
        gPad->Update();
        bias.SetTitle(Form("%d GeV", energy));
        bias.GetXaxis()->SetTitle("#theta [mrad]");
        bias.GetYaxis()->SetTitle("#delta_{#theta} [mrad]  ");
        bias.GetYaxis()->SetTitleOffset(1.31);
        bias.GetYaxis()->SetRangeUser(-0.02, 0.02);
        bias.SetMarkerStyle(kOpenStar);
        RootUtils::SetAllColors(bias, kBlue);

        const double legStartX_G(0.2);
        const double startY(0.9);
        auto* leg = RootUtils::BuildLegend(canvR, legStartX_G, startY, legStartX_G + 0.3, startY - (2) * legendEntryHeight);
        RootUtils::AddLegendHeader(leg, {"Theta Bias"});

        canvR.SaveAs(Form("ThetaVsAngle_%d.eps", energy));
        pdf.AddCanvas(canvR, Form("Theta Bias vs. Angle: %d GeV", energy));
      }

      for (size_t i = THETALOG_START; i < THETALOG_SIZE; ++i) {
        std::vector<double> zeros(angles.size(), 0.0);
        {
          IgnoreError  IE;
          TCanvas      canvR("cE", "cE");
          TGraphErrors resolutions(tAngles.size(), &tAngles[0], &tRMSL[i][0], &zeros[0], &zeros[0]);
          resolutions.Draw("AP");
          gPad->Update();
          resolutions.SetTitle("#sigma_{#theta}");
          resolutions.GetXaxis()->SetTitle("#theta [mrad]");
          resolutions.GetYaxis()->SetTitle("#Delta_{#theta} | #sigma_{#theta} [mrad]  ");
          resolutions.GetYaxis()->SetTitleOffset(1.31);
          resolutions.GetYaxis()->SetRangeUser(-0.05, 0.2);
          resolutions.SetMarkerStyle(kOpenStar);
          RootUtils::SetAllColors(resolutions);

          TGraphErrors bias(tAngles.size(), &tAngles[0], &tMeansL[i][0], &zeros[0], &zeros[0]);
          bias.Draw("P");
          bias.SetTitle("#Delta_{#theta}");
          bias.SetMarkerStyle(kCircle);
          RootUtils::SetAllColors(bias, kBlue);

          const double legStartX_G(0.2);
          const double startY(0.9);
          auto*        leg =
              RootUtils::BuildLegend(canvR, legStartX_G, startY, legStartX_G + 0.3, startY - (4) * legendEntryHeight);
          RootUtils::SetEntryOptions(leg, "p");
          RootUtils::AddLegendHeader(leg, {Form("%d GeV", energy), Form("logConstant = %1.1f", double(i) / 10.0)});

          canvR.SaveAs(Form("ThetaVsAngle_%d.eps", energy));
          pdf.AddCanvas(canvR, Form("%i Theta Res vs. Angle: %d GeV", int(i), energy));
        }
      }

      if (not tAngles.empty()) {
        IgnoreError         IE;
        TCanvas             canvR("cE", "cE");
        std::vector<double> zeros(angles.size(), 0.0);
        TGraphErrors        resolutions(tAngles.size(), &tAngles[0], &tRMS[0], &zeros[0], &zeros[0]);
        resolutions.Draw("AP");
        gPad->Update();
        resolutions.SetTitle("Polar Angle Resolution #sigma_{#theta}");
        resolutions.GetXaxis()->SetTitle("#theta [mrad]");
        resolutions.GetYaxis()->SetTitle("#delta_{#theta} | #sigma_{#theta} [mrad]  ");
        resolutions.GetYaxis()->SetTitleOffset(1.31);
        resolutions.GetYaxis()->SetRangeUser(-0.075, 0.35);
        resolutions.SetMarkerStyle(kOpenStar);

        TGraphErrors bias(tAngles.size(), &tAngles[0], &tMeans[0], &zeros[0], &zeros[0]);
        bias.Draw("P");
        gPad->Update();
        bias.SetTitle("Polar Angle Bias #delta_{#theta}");
        bias.SetMarkerStyle(kCircle);
        RootUtils::SetAllColors(bias, kBlue);

        const double legStartX_G(0.2);
        const double startY(0.9);
        auto* leg = RootUtils::BuildLegend(canvR, legStartX_G, startY, legStartX_G + 0.3, startY - (3) * legendEntryHeight);
        RootUtils::SetEntryOptions(leg, "p");
        RootUtils::AddLegendHeader(leg, {Form("%d GeV, Log Constant: 6.1", energy)});

        canvR.SaveAs(Form("ThetaResVsAngle_%d.eps", energy));
        pdf.AddCanvas(canvR, Form("Theta Res vs. Angle: %d GeV", energy));

        resolutions.GetYaxis()->SetRangeUser(-0.01, 0.05);
        canvR.SaveAs(Form("ThetaResVsAngle_%d_zoom.eps", energy));
      }

      /// phi vs angle
      if (not tAngles.empty()) {
        IgnoreError         IE;
        TCanvas             canvR("cE", "cE");
        std::vector<double> zeros(pAngles.size(), 0.0);
        TGraphErrors        resolutions(pAngles.size(), &pAngles[0], &pMeans[0], &zeros[0], &zeros[0]);
        resolutions.Draw("AP");
        gPad->Update();
        resolutions.SetTitle(Form("%d GeV", energy));
        resolutions.GetXaxis()->SetTitle("#theta [mrad]");
        resolutions.GetYaxis()->SetTitle("#sigma_{#phi} [deg]  ");
        resolutions.GetYaxis()->SetTitleOffset(1.31);
        resolutions.GetYaxis()->SetRangeUser(-0.5, 0.5);
        resolutions.SetMarkerStyle(kOpenStar);

        const double legStartX_G(0.2);
        const double startY(0.9);
        auto* leg = RootUtils::BuildLegend(canvR, legStartX_G, startY, legStartX_G + 0.3, startY - (2) * legendEntryHeight);
        RootUtils::AddLegendHeader(leg, {"Phi Bias"});

        canvR.SaveAs(Form("PhiBiasVsAngle_%d.eps", energy));
        pdf.AddCanvas(canvR, Form("Phi Bias vs. Angle: %d GeV", energy));
      }

      if (not tAngles.empty()) {
        IgnoreError         IE;
        TCanvas             canvR("cE", "cE");
        std::vector<double> zeros(pAngles.size(), 0.0);
        TGraphErrors        resolutions(pAngles.size(), &pAngles[0], &pRMS[0], &zeros[0], &zeros[0]);
        resolutions.Draw("AP");
        gPad->Update();
        resolutions.SetTitle(Form("%d GeV", energy));
        resolutions.GetXaxis()->SetTitle("#theta [mrad]");
        resolutions.GetYaxis()->SetTitle("#sigma_{#phi} [deg]  ");
        resolutions.GetYaxis()->SetTitleOffset(1.31);
        resolutions.GetYaxis()->SetRangeUser(0.0, 5);
        resolutions.SetMarkerStyle(kOpenStar);

        const double legStartX_G(0.2);
        const double startY(0.9);
        auto* leg = RootUtils::BuildLegend(canvR, legStartX_G, startY, legStartX_G + 0.3, startY - (2) * legendEntryHeight);
        RootUtils::AddLegendHeader(leg, {"Phi Resolution"});

        canvR.SaveAs(Form("PhiResVsAngle_%d.eps", energy));
        pdf.AddCanvas(canvR, Form("Phi Res vs. Angle: %d GeV", energy));
      }

    }  //all histograms

  }  // energy related histograms

  {  //Energy Resolution linearity
    if (not selectedEnergyEnergy.empty()) {
      TCanvas             canvER("cER", "cER");
      std::vector<double> zeros(selectedEnergyEnergy.size(), 0.0);
      TGraphErrors resolutions(selectedEnergyEnergy.size(), &selectedEnergyEnergy[0], &selectedEnergyMean[0], &zeros[0],
                               &selectedEnergyVariance[0]);
      resolutions.Draw("AP");
      gPad->Update();
      resolutions.GetXaxis()->SetTitle("E_{MC} [GeV]");
      resolutions.GetYaxis()->SetTitle("E_{Rec} [GeV]");
      resolutions.GetYaxis()->SetTitleOffset(1.31);

      resolutions.Fit("pol1");
      resolutions.SetTitle(" Linearity ");
      const double startY(0.9), legStartX(0.2);
      auto* leg = RootUtils::BuildLegend(canvER, legStartX, startY, legStartX + 0.3, startY - (1 + 2) * legendEntryHeight);
      RootUtils::AddLegendHeader(leg, {"Reconstructed Energy", "#theta #in [62, 77] mrad"});
      IgnoreError IE;
      pdf.AddCanvas(canvER, "Energy Resolution vs. Energy");
      canvER.SaveAs("ERes.eps");
    }

    if (not selectedEnergyEnergy.empty()) {
      TCanvas             canvER("cER", "cER");
      std::vector<double> zeros(selectedEnergyEnergy.size(), 0.0);
      TGraphErrors resolutions(selectedEnergyEnergy.size(), &selectedEnergyEnergy[0], &selectedEnergyRelDiff[0], &zeros[0],
                               &selectedEnergyRelVar[0]);
      resolutions.Draw("APM");
      resolutions.SetMarkerStyle(kFullStar);
      gPad->Update();
      resolutions.GetXaxis()->SetTitle("E_{MC} [GeV]");
      resolutions.GetYaxis()->SetTitle("(E_{Rec} - E_{MC})/E_{MC}");
      resolutions.GetYaxis()->SetRangeUser(-0.016, 0.011);
      resolutions.GetYaxis()->SetTitleOffset(1.31);
      resolutions.SetTitle("Relative Bias");
      const double startY(0.9), legStartX(0.2);
      auto* leg = RootUtils::BuildLegend(canvER, legStartX, startY, legStartX + 0.3, startY - (1 + 2) * legendEntryHeight);
      RootUtils::AddLegendHeader(leg, {"Reconstructed Energy", "#theta #in [62, 77] mrad"});
      IgnoreError IE;
      pdf.AddCanvas(canvER, "Energy Resolution vs. Energy");
      canvER.SaveAs("EResDiff.eps");
    }
  }

  drawEnergyResolutionFunction(pdf, selectedEnergyEnergy, selectedSigmaEoverE, selectedSigmaEoverEErr, fitRes);

  if (allHistograms) {
    for (auto const& pIntMapIntResults : histograms) {
      const int         energy = pIntMapIntResults.first;
      RootUtils::Colors colors;
      TCanvas           canvT("cT", "cT");
      bool              first = true;
      TString           opt("");
      int               histCounter = 0;
      for (auto const& pIntResults : pIntMapIntResults.second) {
        //const int iAngle = pIntResults.first;
        Results* results = pIntResults.second;
        if (not results)
          continue;
        TH1D* tHisto = results->hs["eRes"];
        if (first) {
          first = false;
          opt   = "";
        } else {
          opt = "same";
        }
        histCounter++;
        RootUtils::SetAllColors(tHisto, colors.GetColor());
        canvT.cd();
        tHisto->Draw(opt);
        if (histCounter > 5) {
          saveCanvas(pdf, canvT, energy, Form("Energy vs Angle: %d GeV", energy), histCounter, 0.20);
          histCounter = 0;
          first       = true;
          colors.SetColors();
        }
      }
      if (histCounter > 0) {
        saveCanvas(pdf, canvT, energy, Form("Energy vs Angle: %d GeV", energy), histCounter, 0.20);
      }
    }
  }  //Theta related histograms

  // Draw Theta Related Histograms
  if (allHistograms) {
    std::vector<std::vector<std::string>> titles{
        {"O Theta Resolution vs Angle: %d GeV", "tRes"},
        // {"N Theta Resolution vs Angle: %d GeV", "tResN"},
        // {"E Theta Resolution vs Angle: %d GeV", "tResE"},
        // {"L Theta Resolution vs Angle: %d GeV 0", "tResL0"},
        // {"L Theta Resolution vs Angle: %d GeV 1", "tResL1"},
        // {"L Theta Resolution vs Angle: %d GeV 2", "tResL2"},
        // {"L Theta Resolution vs Angle: %d GeV 3", "tResL3"},
        // {"L Theta Resolution vs Angle: %d GeV 4", "tResL4"},
        // {"L Theta Resolution vs Angle: %d GeV 5", "tResL5"},
        // {"L Theta Resolution vs Angle: %d GeV 6", "tResL6"},
        // {"L Theta Resolution vs Angle: %d GeV 7", "tResL7"},
        // {"L Theta Resolution vs Angle: %d GeV 8", "tResL8"},
        // {"L Theta Resolution vs Angle: %d GeV 9", "tResL9"},
        // {"L Theta Resolution vs Angle: %d GeV 10", "tResL10"},
        // {"L Theta Resolution vs Angle: %d GeV 11", "tResL11"},
        // {"L Theta Resolution vs Angle: %d GeV 12", "tResL12"},
        // {"L Theta Resolution vs Angle: %d GeV 13", "tResL13"},
        // {"L Theta Resolution vs Angle: %d GeV 14", "tResL14"},
        // {"L Theta Resolution vs Angle: %d GeV 15", "tResL15"},
        // {"L Theta Resolution vs Angle: %d GeV 16", "tResL16"},
        // {"L Theta Resolution vs Angle: %d GeV 17", "tResL17"},
        // {"L Theta Resolution vs Angle: %d GeV 18", "tResL18"},
        // {"L Theta Resolution vs Angle: %d GeV 19", "tResL19"},
        // {"L Theta Resolution vs Angle: %d GeV 20", "tResL20"},
        // {"L Theta Resolution vs Angle: %d GeV 21", "tResL21"},
        // {"L Theta Resolution vs Angle: %d GeV 22", "tResL22"},
        // {"L Theta Resolution vs Angle: %d GeV 23", "tResL23"},
        // {"L Theta Resolution vs Angle: %d GeV 24", "tResL24"},
        // {"L Theta Resolution vs Angle: %d GeV 25", "tResL25"},
        // {"L Theta Resolution vs Angle: %d GeV 26", "tResL26"},
        // {"L Theta Resolution vs Angle: %d GeV 27", "tResL27"},
        // {"L Theta Resolution vs Angle: %d GeV 28", "tResL28"},
        // {"L Theta Resolution vs Angle: %d GeV 29", "tResL29"},

        {"L Theta Resolution vs Angle: %d GeV 30", "tResL30"},
        {"L Theta Resolution vs Angle: %d GeV 31", "tResL31"},
        {"L Theta Resolution vs Angle: %d GeV 32", "tResL32"},
        {"L Theta Resolution vs Angle: %d GeV 33", "tResL33"},
        {"L Theta Resolution vs Angle: %d GeV 34", "tResL34"},
        {"L Theta Resolution vs Angle: %d GeV 35", "tResL35"},
        {"L Theta Resolution vs Angle: %d GeV 36", "tResL36"},
        {"L Theta Resolution vs Angle: %d GeV 37", "tResL37"},
        {"L Theta Resolution vs Angle: %d GeV 38", "tResL38"},
        {"L Theta Resolution vs Angle: %d GeV 39", "tResL39"},
        {"L Theta Resolution vs Angle: %d GeV 40", "tResL40"},
        {"L Theta Resolution vs Angle: %d GeV 41", "tResL41"},
        {"L Theta Resolution vs Angle: %d GeV 42", "tResL42"},
        {"L Theta Resolution vs Angle: %d GeV 43", "tResL43"},
        {"L Theta Resolution vs Angle: %d GeV 44", "tResL44"},
        {"L Theta Resolution vs Angle: %d GeV 45", "tResL45"},
        {"L Theta Resolution vs Angle: %d GeV 46", "tResL46"},
        {"L Theta Resolution vs Angle: %d GeV 47", "tResL47"},
        {"L Theta Resolution vs Angle: %d GeV 48", "tResL48"},
        {"L Theta Resolution vs Angle: %d GeV 49", "tResL49"},
        {"L Theta Resolution vs Angle: %d GeV 50", "tResL50"},
        {"L Theta Resolution vs Angle: %d GeV 51", "tResL51"},
        {"L Theta Resolution vs Angle: %d GeV 52", "tResL52"},
        {"L Theta Resolution vs Angle: %d GeV 53", "tResL53"},
        {"L Theta Resolution vs Angle: %d GeV 54", "tResL54"},
        {"L Theta Resolution vs Angle: %d GeV 55", "tResL55"},
        {"L Theta Resolution vs Angle: %d GeV 56", "tResL56"},
        {"L Theta Resolution vs Angle: %d GeV 57", "tResL57"},
        {"L Theta Resolution vs Angle: %d GeV 58", "tResL58"},
        {"L Theta Resolution vs Angle: %d GeV 59", "tResL59"},
        {"L Theta Resolution vs Angle: %d GeV 61", "tResL61"},
        {"L Theta Resolution vs Angle: %d GeV 62", "tResL62"},
        {"L Theta Resolution vs Angle: %d GeV 63", "tResL63"},
        {"L Theta Resolution vs Angle: %d GeV 64", "tResL64"},
        {"L Theta Resolution vs Angle: %d GeV 65", "tResL65"},
        {"L Theta Resolution vs Angle: %d GeV 66", "tResL66"},
        {"L Theta Resolution vs Angle: %d GeV 67", "tResL67"},
        {"L Theta Resolution vs Angle: %d GeV 68", "tResL68"},
        {"L Theta Resolution vs Angle: %d GeV 69", "tResL69"},
        {"L Theta Resolution vs Angle: %d GeV 70", "tResL70"},
        {"L Theta Resolution vs Angle: %d GeV 71", "tResL71"},
        {"L Theta Resolution vs Angle: %d GeV 72", "tResL72"},
        {"L Theta Resolution vs Angle: %d GeV 73", "tResL73"},
        {"L Theta Resolution vs Angle: %d GeV 74", "tResL74"},
        {"L Theta Resolution vs Angle: %d GeV 75", "tResL75"},
        {"L Theta Resolution vs Angle: %d GeV 76", "tResL76"},
        {"L Theta Resolution vs Angle: %d GeV 77", "tResL77"},
        {"L Theta Resolution vs Angle: %d GeV 78", "tResL78"},
        {"L Theta Resolution vs Angle: %d GeV 79", "tResL79"},
        {"L Theta Resolution vs Angle: %d GeV 80", "tResL80"},
        {"L Theta Resolution vs Angle: %d GeV 81", "tResL81"},
        {"L Theta Resolution vs Angle: %d GeV 82", "tResL82"},
        {"L Theta Resolution vs Angle: %d GeV 83", "tResL83"},
        {"L Theta Resolution vs Angle: %d GeV 84", "tResL84"},
        {"L Theta Resolution vs Angle: %d GeV 85", "tResL85"},
        {"L Theta Resolution vs Angle: %d GeV 86", "tResL86"},
        {"L Theta Resolution vs Angle: %d GeV  87", "tResL87"},
        {"L Theta Resolution vs Angle: %d GeV  88", "tResL88"},
        {"L Theta Resolution vs Angle: %d GeV  89", "tResL89"},
        {"L Theta Resolution vs Angle: %d GeV  90", "tResL90"},
        {"L Theta Resolution vs Angle: %d GeV  91", "tResL91"},
        {"L Theta Resolution vs Angle: %d GeV  92", "tResL92"},
        {"L Theta Resolution vs Angle: %d GeV  93", "tResL93"},
        {"L Theta Resolution vs Angle: %d GeV  94", "tResL94"},
        {"L Theta Resolution vs Angle: %d GeV  95", "tResL95"},
        {"L Theta Resolution vs Angle: %d GeV  96", "tResL96"},
        {"L Theta Resolution vs Angle: %d GeV  97", "tResL97"},
        {"L Theta Resolution vs Angle: %d GeV  98", "tResL98"},
        {"L Theta Resolution vs Angle: %d GeV  99", "tResL99"},
        //        {"L Theta Resolution vs Angle: %d GeV 100", "tResL100"},

    };

    for (auto const& pIntMapIntResults : histograms) {
      for (auto const& entries : titles) {
        std::string const& tTitle = entries[0];
        std::string const& hTitle = entries[1];
        const int          energy = pIntMapIntResults.first;
        RootUtils::Colors  colors;
        TCanvas            canvT("cT", "cT");
        bool               first = true;
        TString            opt("");
        int                histCounter = 0;
        for (auto const& pIntResults : pIntMapIntResults.second) {
          //const int iAngle = pIntResults.first;
          Results* results = pIntResults.second;
          if (not results)
            continue;
          TH1D* tHisto = results->hs[hTitle];
          if (first) {
            first = false;
            opt   = "";
          } else {
            opt = "same";
          }
          histCounter++;
          RootUtils::SetAllColors(tHisto, colors.GetColor());
          canvT.cd();
          tHisto->Draw(opt);
          if (histCounter > 100) {
            saveCanvas(pdf, canvT, energy, Form(tTitle.c_str(), energy), histCounter);
            histCounter = 0;
            first       = true;
            colors.SetColors();
          }
        }
        if (histCounter > 0) {
          saveCanvas(pdf, canvT, energy, Form(tTitle.c_str(), energy), histCounter);
        }
      }
    }
  }  //Theta related histograms

  // Draw Phi Related Histograms
  if (allHistograms) {
    for (auto const& pIntMapIntResults : histograms) {
      RootUtils::Colors colors;
      const int         energy = pIntMapIntResults.first;
      TCanvas           canv("cP", "cP");
      bool              first = true;
      TString           opt("");
      int               histCounter = 0;
      for (auto const& pIntResults : pIntMapIntResults.second) {
        //const int iAngle = pIntResults.first;
        Results* results = pIntResults.second;
        if (not results)
          continue;
        TH1D* tHisto = results->hs["pRes"];
        if (first) {
          first = false;
          opt   = "";
        } else {
          opt = "same";
        }
        histCounter++;
        RootUtils::SetAllColors(tHisto, colors.GetColor());
        canv.cd();
        tHisto->Draw(opt);
        if (histCounter > 5) {
          saveCanvas(pdf, canv, energy, Form("Phi Resolution vs Angle: %d GeV", energy), histCounter);
          histCounter = 0;
          first       = true;
          colors.SetColors();
        }
      }
      if (histCounter > 0) {
        saveCanvas(pdf, canv, energy, Form("Phi Resolution vs Angle: %d GeV", energy), histCounter);
      }
    }
  }  //Phi related histograms

  drawEfficiencies(pdf, efficiencies);

  /// Finish the pdf file with an empty slide to keep the last title
  TCanvas empty("e", "e");
  pdf.AddCanvas(empty, "end");

  return 0;
}

void drawEfficiencies(RootUtils::PDFFile& pdf, Efficiencies& efficiencies) {
  for (auto& pIntEfficiency : efficiencies) {
    const int energy = pIntEfficiency.first;
    TCanvas   canv("cE", "cE");
    pIntEfficiency.second->Draw();

    IgnoreError  IE;
    const double startY    = 0.4;
    const double legStartX = 0.2;
    RootUtils::BuildLegend(canv, legStartX, startY, legStartX + 0.3, startY - (1) * legendEntryHeight);
    //RootUtils::AddLegendHeader( leg, {Form("%d GeV", energy )} );
    TString title = Form("Reconstruction Efficiency %d GeV", energy);
    pdf.AddCanvas(canv, title);
    TString filename = Form("RecoEff%d.eps", energy);
    canv.SaveAs(filename);
  }
  return;
}  //end drawEfficiencies

/// draw the sigmaE/E function prop to a/sqrt(E)
void drawEnergyResolutionFunction(RootUtils::PDFFile& pdf, VecD& selectedEnergyEnergy, VecD& selectedSigmaEoverE,
                                  VecD& selectedSigmaEoverEErr, TF1& fitRes) {
  if (not selectedEnergyEnergy.empty()) {
    TCanvas             canvER("cER", "cER");
    std::vector<double> zeros(selectedEnergyEnergy.size(), 0.0);
    TGraphErrors resolutions(selectedEnergyEnergy.size(), &selectedEnergyEnergy[0], &selectedSigmaEoverE[0], &zeros[0],
                             &selectedSigmaEoverEErr[0]);
    resolutions.Draw("AP");
    gPad->Update();
    resolutions.GetXaxis()->SetTitle("E_{MC} [GeV]");
    resolutions.GetYaxis()->SetTitle("#sigma_{E}/E_{MC}");
    resolutions.Fit(&fitRes);
    resolutions.SetMarkerStyle(kOpenStar);
    resolutions.SetTitle("Selected Angles");
    const double startY(0.9), legStartX(0.3);
    auto*        leg = RootUtils::BuildLegend(canvER, legStartX, startY, legStartX + 0.3, startY - (3) * legendEntryHeight);
    RootUtils::AddLegendHeader(leg, {"Energy Resolution", "a/#sqrt{E}"});
    IgnoreError IE;
    pdf.AddCanvas(canvER, "SigmaE/E");
    canvER.SaveAs("EResFunc.eps");
  }

}  // end drawEnergyResolutionFunction
