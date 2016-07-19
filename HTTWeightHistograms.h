#ifndef HTTWeightHistograms_h
#define HTTWeightHistograms_h

// Original Author:  Artur Kalinowski i ja :D
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include <string>
#include "AnalysisHistograms.h"

class THStack;


class HTTWeightHistograms: public AnalysisHistograms {
   public:

  HTTWeightHistograms(std::string fileName="Histos.root", int opt=0);

  HTTWeightHistograms(TFileDirectory *myDir);

  HTTWeightHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HTTWeightHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  virtual bool fill1DHistogram(const std::string &name, float val, float weight=1.0);

  float getLumi();

  ///Return sample normalisation based only on
  ///luminosity and cross section.
  ///MC to DATA scaling factors should be applied
  ///on top of this normalisation.
  float getSampleNormalisation(std::string sampleName);

  ///Calculate scaling factor for the WJets MC
  ///SCaling factor is estimated in high Mt region.
  ///Other backgrounds are subtracted, basing on MC
  ///QCD contribution is neglected.
  std::pair<float,float> getWNormalisation(std::string selName, std::string SubSelName);

   private:
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  ///Return histogram for sum of all DY decay modes.
  TH1F *get1D_DY_Histogram(const std::string& name);

  ///Return histogram for sum of all WJet HT bins
  TH1F *get1D_WJet_Histogram(const std::string& name);


  ///Parameter SubSelName was used to plot histograms for different, more sophisticated analysis.For simple asymmetry analysis use "Plus".
  ///Prepare histograms: Plus-Minus; Plus+Minus
  double MakeDiff(TH1F *hTTbar, TH1F* hDYJets, TH1F* hSoup, TH1F* hWJets, TH1F* hQCD, TH1F *hTTbarS, TH1F* hDYJetsS, TH1F* hSoupS, TH1F* hWJetsS, TH1F* hQCDS, std::string varName, std::string selName, std::string SubSelName);

  ///Prepare asymmetry histograms and W+Jet background estimation
  std::pair<std::pair<TH1*,TH1*>,TH1*> PlotAsymm(std::string varName, std::string selName, std::string SubSelName);

  ///Plots W+Jet estimation
  double * WJetEstimation(std::string varName, std::string selName, std::string SubSelName);

  void PlotTwoHistograms(std::string varName,std::string nazwa, TH1* th1,std::string opis1, TH1* th2, std::string opis2);
  void PlotOneHistogram(std::string varName,std::string nazwa, TH1* th1,std::string opis1);

  std::pair<float,float> getQCDOStoSS(std::string selName, std::string SubSelName);
  TH1F* getQCDbackground(std::string varName, std::string selName, std::string SubSelName);

  ///Plot stack histograms for Plus-Minus
  THStack*  PlotDiff(std::string varName, std::string selName, std::string SubSelName);
};

#endif
