#ifndef HTTHistogramsTT_h
#define HTTHistogramsTT_h

// Original Author:  Michal Bluj
//         Created:  pon, 21 lis 2016, 13:33:31 CEST
//
//
#include <string>

#include "AnalysisHistograms.h"

class THStack;

class HTTHistogramsTT: public AnalysisHistograms {
 public:

  HTTHistogramsTT(TDirectory *myDir);

  HTTHistogramsTT(TDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HTTHistogramsTT();

  void finalizeHistograms(int nRuns, float weight=1.0);

  std::string getTemplateName(const std::string& name);
  
  float getLumi();

  ///Return sample normalisation based only on
  ///luminosity and cross section.
  ///MC to DATA scaling factors should be applied
  ///on top of this normalisation.
  float getSampleNormalisation(std::string sampleName);

  ///Estimate QCD background using the Loose/Tight method.
  TH1F* getQCDbackground(unsigned int iCategory, std::string varName);

  ///Calculate QCD OS/SS ratio using non isolated events.
  std::pair<float,float> getQCDLooseToTight(unsigned int iCategory);

 private:

  std::pair<float,float> wselOSCorrection;
  std::pair<float,float> wselSSCorrection;
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  //Plot stacked histograms for each contributing process.
  ///iCategory - selection category to be plotted.
  ///selName - secondary type of selection (OS/SS/mt) used for background estimation
  //varName - name of variable to be plotted,
  THStack* plotStack(unsigned int iCategory,
		     std::string varName, 
		     std::string selName = "OS");

  void plotnPCA(const std::string & type);

  void plotVerticesPulls(const std::string & hName);

  void plotProfiles(const std::string & hName,
		    const std::string & sysType);
  
  void plotPhiDecayPlanes(const std::string& name);

  void plot_HAZ_Histograms(const std::string & hName,
			   const std::string & sysType);

  void plotCPhistograms(int nRuns, float weight);

  ///Return histogram for sum of all DY decay modes, and jet bins
  TH1F *get1D_DYJet_Histogram(const std::string& name);

  ///Return histogram for sum of all jet bins
  TH1F *get1D_WJet_Histogram(const std::string& name);

  ///Return histogram for sum VV processes
  TH1F *get1D_VV_Histogram(const std::string& name);

  ///Return histogram for sum single top processes
  TH1F *get1D_ST_Histogram(const std::string& name);

  ///Return sum of DY histograms. Sum can run over
  ///decay modes, jet bins, or both.
  TH1F *get1D_DYSum(const std::string& name, bool sumDecayModes, bool sumJetBins);

  ///Return histogram for sum of all W/Z jet bins
  ///The results is scaled to 1/LO_xsection.
  TH1F *get1D_VJetSum(const std::string& name);

  ///Return histogram from nJets sample normalised by
  ///preselection/number of analysed events
  TH1F *getNormalised_NJet_Histogram(const std::string& hName);

  //Plot a single histogram. One has to provide the full
  //histogram name, e.g. including h1D prefix.
  void plotSingleHistogram(std::string hName);

  //float tauTauDYScale, lepTauDYScale, jetTauDYScale;//MB not used
  float ttScale;

};

#endif
