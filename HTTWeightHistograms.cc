#include <iostream>
#include <cmath>

#include "commonUtils.h"
#include "HTTWeightHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TMath.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THStack.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
float HTTWeightHistograms::getLumi(){

  //pileupCalc.py -i 15_12_2015.json --inputLumiJSON pileup_latest.txt --calcMode observed --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 MyDataPileupHistogram.root

  //brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i lumiSummary.json
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 19    | 36   | 10256 | 10256 | 573248145.792     | 552672886.226    |

  //return 552672886.226e-6;//pb-1 data for NTUPLES_23_11_2015

  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i 15_12_2015.json
  //+-------+------+-------+-------+-------------------+------------------+
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 47    | 115  | 30707 | 30707 | 2135679014.929    | 2066764067.818   |
  //+-------+------+-------+-------+-------------------+------------------+

  //return 2066764067.818e-6;//pb-1 data for NTUPLES_15_12_2015


  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i 18_12_2015.json
  //+-------+------+-------+-------+-------------------+------------------+
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 47    | 116  | 32567 | 32567 | 2282566714.405    | 2207823354.548   |
  //+-------+------+-------+-------+-------------------+------------------+

  //return 2207823354.548e-6;//pb-1 data for NTUPLES_18_12_2015
  
  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i 12_02_2016_CD.json
  //+-------+------+-------+-------+-------------------+------------------+
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 47    | 116  | 32702 | 32702 | 2289599379.230    | 2214575536.754   |
  //+-------+------+-------+-------+-------------------+------------------+
  return 2214575536.754e-6;//pb-1 data for NTUPLES_12_02_2016
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
float HTTWeightHistograms::getSampleNormalisation(std::string sampleName){

  std::string hName = "h1DStats"+sampleName;
  TH1F *hStats = 0;
  if(sampleName=="DYJets") hStats = get1D_DY_Histogram(hName.c_str());
  else if(sampleName=="WJets") hStats = get1D_WJet_Histogram(hName.c_str());
  else hStats = get1DHistogram(hName.c_str());

  if(!hStats) return 0;

  float genPresEff = 1.0;
  float recoPresEff = hStats->GetBinContent(3)/hStats->GetBinContent(2);
  float presEff = genPresEff*recoPresEff;
  float kFactor = 1.0;

  float crossSection = 1.0;
  int nEventsAnalysed = hStats->GetBinContent(1);

  ///Cross sections taken from
  if(sampleName=="DYJetsLowM"){
    //https://cmsweb.cern.ch/das/request?input=mcm%20prepid=SMP-RunIISpring15MiniAODv2-00016
    crossSection = 71600;
  }
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
  if(sampleName=="DYJets"){
    //xsection for 3xZ->mu mu M50 in [pb]  
    crossSection = 3*2008.4; 
  }
  if(sampleName=="WJets"){
    //xsection for 3xW->mu nu in [pb]
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
    crossSection = 3*20508.9;
  }

  if(sampleName=="TTbar"){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
    crossSection = 831.76*0.95; 
  }
  if(sampleName=="H"){
    ///mH = 125, gg fussion only
    //https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx
    //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
    crossSection = 44.14*6.32E-02; 
  }
  
  float weight = crossSection*presEff/nEventsAnalysed;
  if(presEff<0 || fabs(fabs(crossSection)-1.0)<1e-5) weight = 1.0;

  std::cout<<"Sample name: "<<sampleName<<" ";
  std::cout<<"Xsection: "<<crossSection<<" [pb] "<<" ";
  std::cout<<"Events analyzed: "<<nEventsAnalysed<<" ";
  //std::cout<<"Gen preselection efficiency: "<<genPresEff<<" ";
  std::cout<<"Reco preselection efficiency: "<<recoPresEff<<" ";
  //std::cout<<"External scaling: "<<kFactor<<" ";
  std::cout<<"Final weight: "<<weight<<std::endl;

  return weight;  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::HTTWeightHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::HTTWeightHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::HTTWeightHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){

  selectionFlavours_ = flavours;

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::~HTTWeightHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTWeightHistograms::get1D_DY_Histogram(const std::string& name){
  
  TString hName = name;
  hName.ReplaceAll("DYJets","DYJetsMuTau");
  TH1F *hDYJetsMuTau = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("DYJets","DYJetsMuMu");
  TH1F *hDYJetsMuMu = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("DYJets","DYJetsEE");
  TH1F *hDYJetsEE = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("DYJets","DYJetsOther");
  TH1F *hDYJetsOther = get1DHistogram(hName.Data());

  if(!hDYJetsMuTau) return 0;
  if(hDYJetsMuMu) hDYJetsMuTau->Add(hDYJetsMuMu);
  if(hDYJetsEE) hDYJetsMuTau->Add(hDYJetsEE);
  if(hDYJetsOther) hDYJetsMuTau->Add(hDYJetsOther);
  hDYJetsMuTau->SetName(name.c_str());

  return hDYJetsMuTau;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTWeightHistograms::get1D_WJet_Histogram(const std::string& name){

  TString hName = name;
  hName.ReplaceAll("WJets","WJetsHT0");
  TH1F *hWJets0to100 = get1DHistogram(hName.Data());
  
  hName = name;
  hName.ReplaceAll("WJets","WJetsHT100to200");
  TH1F *hWJets100to200 = get1DHistogram(hName.Data());
    
  hName = name;
  hName.ReplaceAll("WJets","WJetsHT200to400");
  TH1F *hWJets200to400 = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("WJets","WJetsHT400to600");
  TH1F *hWJets400to600 = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("WJets","WJetsHT600toInf");
  TH1F *hWJets600toInf = get1DHistogram(hName.Data());

  TH1F *hWJets = 0;
  if(hWJets0to100) hWJets = (TH1F*)hWJets0to100->Clone(name.c_str());
  else if(hWJets100to200) hWJets = (TH1F*)hWJets100to200->Clone(name.c_str());
  else if(hWJets200to400) hWJets = (TH1F*)hWJets200to400->Clone(name.c_str());

  if(!hWJets) return 0;
  
  hWJets->Clear();
  if(hWJets0to100) hWJets->Add(hWJets0to100);
  if(hWJets100to200) hWJets->Add(hWJets100to200);
  if(hWJets200to400) hWJets->Add(hWJets200to400);
  if(hWJets400to600) hWJets->Add(hWJets400to600);
  if(hWJets600toInf) hWJets->Add(hWJets600toInf);

  return hWJets;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool HTTWeightHistograms::fill1DHistogram(const std::string& name, float val, float weight){

  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill1DHistogram(name,val,weight)){
    if(name.find("h1DMass")!=std::string::npos) hTemplateName = "h1DMassTemplate";
    if(name.find("h1DStats")!=std::string::npos) hTemplateName = "h1DStatsTemplate";
    if(name.find("h1DEta")!=std::string::npos) hTemplateName = "h1DEtaTemplate";
    if(name.find("h1DIso")!=std::string::npos) hTemplateName = "h1DIsoTemplate";
    
    if(get1DHistogram(hTemplateName,true)->GetXaxis()->IsVariableBinSize()){
      Float_t* binsArray = new Float_t[this->get1DHistogram(hTemplateName,true)->GetNbinsX()+1];
      for(unsigned int iBin=0;iBin<=this->get1DHistogram(hTemplateName,true)->GetNbinsX();++iBin){
	binsArray[iBin] = this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXbins()->At(iBin);
      }
      this->add1DHistogram(name,"",this->get1DHistogram(hTemplateName,true)->GetNbinsX(),
			   binsArray, file_);
      this->get1DHistogram(name,true)->SetDirectory(this->get1DHistogram(hTemplateName,true)->GetDirectory());
      delete binsArray;      
    }
    else{
      this->add1DHistogram(name,"",
			   this->get1DHistogram(hTemplateName,true)->GetNbinsX(),
			   this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmin(),
			   this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmax(),
			   file_);
      this->get1DHistogram(name,true)->SetDirectory(this->get1DHistogram(hTemplateName,true)->GetDirectory());
    }   
    return AnalysisHistograms::fill1DHistogram(name,val,weight);
  }
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   std::cout<<"defineHistograms Adding histogram: "<<file_<<" "<<file_->fullPath()<<std::endl;

   add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
   add1DHistogram("h1DMassTemplate",";SVFit mass [GeV/c^{2}]; Events",50,0,200,file_);
   add1DHistogram("h1DEtaTemplate",";#eta; Events",24,-2.4,2.4,file_);
   float bins[31] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5};
   add1DHistogram("h1DIsoTemplate",";Isolation; Events",30,bins,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

  TFile f("AsymmWeights.root","recreate");

   // control histograms
  WJetEstimation("MassTrans","","Plus");
  WJetEstimation("MassTrans","All","Plus");
  WJetEstimation("MassTrans","qcdselSS","Plus");
  WJetEstimation("EtaMuon","","Plus");
  WJetEstimation("EtaMuon","All","Plus");
  WJetEstimation("EtaMuon","qcdselSS","Plus");

  PlotDiff("MassVis","","Plus");
  PlotDiff("MassVis","qcdselSS","Plus");
  PlotDiff("MassVis","All","Plus");

  PlotDiff("MassTrans","","Plus");
  PlotDiff("MassTrans","qcdselSS","Plus");
  PlotDiff("MassTrans","All","Plus");

   // asymmetry histograms 
  std::pair<std::pair<TH1*,TH1*>, TH1*> WEstimation = PlotAsymm("EtaMuon","","Plus");
  TH1F * hAsymmEtaMC = (TH1F*) WEstimation.second;
  hAsymmEtaMC -> SetName("EtaAsymmetryOS");
  hAsymmEtaMC -> Write(); 

  WEstimation = PlotAsymm("EtaMuon","qcdselSS","Plus");
  TH1F * hAsymmEtaMCSS = (TH1F*) WEstimation.second;
  hAsymmEtaMCSS -> SetName("EtaAsymmetrySS");
  hAsymmEtaMCSS -> Write(); 

  WEstimation = PlotAsymm("MassTrans","qcdselSS","Plus");
  TH1F * hAsymmMassTransMCSS = (TH1F*) WEstimation.second;
  hAsymmMassTransMCSS -> SetName("MassTransAsymmetrySS");
  hAsymmMassTransMCSS -> Write(); 

  WEstimation = PlotAsymm("EtaMuon","All","Plus");
  TH1F * hAsymmEtaMCAll = (TH1F*) WEstimation.second;
  hAsymmEtaMCAll -> SetName("EtaAsymmetryAll");
  hAsymmEtaMCAll -> Write(); 

  WEstimation = PlotAsymm("MassTrans","All","Plus");
  TH1F * hAsymmMassTransMCAll = (TH1F*) WEstimation.second;
  hAsymmMassTransMCAll -> SetName("MassTransAsymmetryAll");
  hAsymmMassTransMCAll -> Write(); 

  WEstimation = PlotAsymm("MassTrans","","Plus");
  TH1F * hAsymmMassTransMC = (TH1F*) WEstimation.second;
  hAsymmMassTransMC -> SetName("MassTransAsymmetryOS");
  hAsymmMassTransMC -> Write(); 

  f.Close();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double HTTWeightHistograms::MakeDiff(TH1F *hTTbar, TH1F* hDYJets, TH1F* hSoup, TH1F* hWJets, TH1F* hQCD, TH1F *hTTbarS, TH1F* hDYJetsS, TH1F* hSoupS, TH1F* hWJetsS, TH1F* hQCDS, std::string varName, std::string selName, std::string SubSelName){

  std::string selName2="";

  if(selName.find("All")!=std::string::npos) {
	selName="";
	selName2="All";
	}

  std::string hName = "h1D"+varName;
  hTTbar -> Reset();
  hDYJets-> Reset();
  hSoup  -> Reset();
  hWJets -> Reset();
  hQCD   -> Reset();

  hTTbarS -> Reset();
  hDYJetsS-> Reset();
  hSoupS  -> Reset();
  hWJetsS -> Reset();
  hQCDS   -> Reset();

// added histograms
  TH1F *hTTbar1 = get1DHistogram((hName+"TTbar"+selName+"Plus").c_str());
  TH1F *hSoup1 = get1DHistogram((hName+"Data"+selName+"Plus").c_str());
  TH1F *hWJets11 = get1D_WJet_Histogram((hName+"WJets"+selName+"Plus").c_str());
  TH1F *hDYJets1 = get1D_DY_Histogram((hName+"DYJets"+selName+"Plus").c_str());
  TH1F *hDYJetsLowM1 = get1DHistogram((hName+"DYJetsLowM"+selName+"Plus").c_str());

  if(!hDYJets1){
    hDYJets1 = (TH1F*)hSoup1->Clone((hName+"hDYJets1"+selName).c_str()); hDYJets1->Reset();
  }
  if(!hDYJetsLowM1){
    hDYJetsLowM1 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM1"+selName).c_str()); hDYJetsLowM1->Reset();
  }

  TH1F *hTTbar2 = get1DHistogram((hName+"TTbar"+selName+"Minus").c_str()); 
  TH1F *hSoup2 = get1DHistogram((hName+"Data"+selName+"Minus").c_str());
  TH1F *hWJets21 = get1D_WJet_Histogram((hName+"WJets"+selName+"Minus").c_str());
  TH1F *hDYJets2 = get1D_DY_Histogram((hName+"DYJets"+selName+"Minus").c_str());
  TH1F *hDYJetsLowM2 = get1DHistogram((hName+"DYJetsLowM"+selName+"Minus").c_str());

  if(!hDYJets2){
    hDYJets2 = (TH1F*)hSoup1->Clone((hName+"hDYJets2"+selName).c_str()); hDYJets2->Reset();
  }
  if(!hDYJetsLowM2){
    hDYJetsLowM2 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM2"+selName).c_str()); hDYJetsLowM2->Reset();
  }

// sum "OS" and "SS" histograms for "All" asymmetry

  if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
  TH1F *hSoup12 = get1DHistogram((hName+"Data"+"qcdselSS"+"Plus").c_str());
  TH1F *hSoup22 = get1DHistogram((hName+"Data"+"qcdselSS"+"Minus").c_str());
  hSoup1->Add(hSoup12,1);
  hSoup2->Add(hSoup22,1);
  }

  float lumi = getLumi();
  std::string sampleName = "WJets";
  float weight=1;
  float scale =1;
  weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
  hWJets11->Scale(scale);
        scale = weight*lumi;
  hWJets21->Scale(scale);

  if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
  TH1F *hWJets12 = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+"Plus").c_str());
  TH1F *hWJets22 = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+"Minus").c_str());
    scale = weight*lumi;
    hWJets12-> Scale(scale);
    hWJets11 -> Add(hWJets12,1);
    hWJets22-> Scale(scale);
    hWJets21 -> Add(hWJets22,1);
  }

  sampleName = "DYJetsLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowM1->Scale(scale);
  hDYJetsLowM2->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets1->Scale(scale);
  hDYJets2->Scale(scale);

 if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
  TH1F *hDYJets12 = get1D_DY_Histogram((hName+"DYJets"+selName+"qcdselSSPlus").c_str());
  TH1F *hDYJetsLowM12 = get1DHistogram((hName+"DYJetsLowM"+selName+"qcdselSSPlus").c_str());

  if(!hDYJets12){
    hDYJets12 = (TH1F*)hSoup1->Clone((hName+"hDYJets12"+selName).c_str()); hDYJets12->Reset();
  }
  if(!hDYJetsLowM12){
    hDYJetsLowM12 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM12"+selName).c_str()); hDYJetsLowM12->Reset();
  } 

  TH1F *hDYJets22 = get1D_DY_Histogram((hName+"DYJets"+selName+"qcdselSSMinus").c_str());
  TH1F *hDYJetsLowM22 = get1DHistogram((hName+"DYJetsLowM"+selName+"qcdselSSMinus").c_str());

  if(!hDYJets22){
    hDYJets22 = (TH1F*)hSoup1->Clone((hName+"hDYJets22"+selName).c_str()); hDYJets22->Reset();
  }
  if(!hDYJetsLowM22){
    hDYJetsLowM22 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM22"+selName).c_str()); hDYJetsLowM22->Reset();
  } 

  std::string sampleName = "DYJetsLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowM12->Scale(scale);
  hDYJetsLowM22->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets12->Scale(scale);
  hDYJets22->Scale(scale);

  hDYJetsLowM1->Add(hDYJetsLowM12,1);
  hDYJets1->Add(hDYJets12,1);
  hDYJetsLowM2->Add(hDYJetsLowM22,1);
  hDYJets2->Add(hDYJets22,1);
  }

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar1->Scale(scale);
  hTTbar2->Scale(scale);

 if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
  TH1F *hTTbar12 = get1DHistogram((hName+"TTbar"+"qcdselSS"+"Plus").c_str());
  TH1F *hTTbar22 = get1DHistogram((hName+"TTbar"+"qcdselSS"+"Minus").c_str());
  hTTbar12->Scale(scale);
  hTTbar22->Scale(scale);
  hTTbar1->Add(hTTbar12,1);
  hTTbar2->Add(hTTbar22,1);
  }

// QCD background
  TH1F * hQCD1 = (TH1F*) getQCDbackground(varName,selName,"Plus");
  TH1F * hQCD2 = (TH1F*) getQCDbackground(varName,selName,"Minus");
  
  if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
	TH1F * hQCD12 = (TH1F*) getQCDbackground(varName,"qcdselSS","Plus");
	TH1F * hQCD22 = (TH1F*) getQCDbackground(varName,"qcdselSS","Minus");

	hQCD1->Add(hQCD12,1);
	hQCD2->Add(hQCD22,1);
	}

// make difference hist (Minus)
  hTTbar -> Add(hTTbar1,1);
  hTTbar -> Add(hTTbar2,-1);
  hDYJets-> Add(hDYJetsLowM1,1);
  hDYJets-> Add(hDYJetsLowM2,-1);
  hDYJets-> Add(hDYJets1,1);
  hDYJets-> Add(hDYJets2,-1);
  hSoup  -> Add(hSoup1,1);
  hSoup  -> Add(hSoup2,-1);
  hWJets -> Add(hWJets11,1);
  hWJets -> Add(hWJets21,-1);  
  hQCD   -> Add(hQCD1,1);
  hQCD   -> Add(hQCD2,-1);

// make sum hist (Plus)
  hTTbarS -> Add(hTTbar1,1);
  hTTbarS -> Add(hTTbar2,1);
  hDYJetsS-> Add(hDYJetsLowM1,1);
  hDYJetsS-> Add(hDYJetsLowM2,1);
  hDYJetsS-> Add(hDYJets1,1);
  hDYJetsS-> Add(hDYJets2,1);
  hSoupS  -> Add(hSoup1,1);
  hSoupS  -> Add(hSoup2,1);
  hWJetsS -> Add(hWJets11,1);
  hWJetsS -> Add(hWJets21,1);
  hQCDS   -> Add(hQCD1,1);
  hQCDS   -> Add(hQCD2,1);

return 0;
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
std::pair<std::pair<TH1*,TH1*>,TH1*>  HTTWeightHistograms::PlotAsymm(std::string varName, std::string selName, std::string SubSelName){

  std::string hName = "h1D"+varName;
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());
  TH1F *hDYJets = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());   
  TH1F *hSoup = get1DHistogram((hName+"Data"+"Diff"+"Minus").c_str());
  TH1F *hWJets = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());
  TH1F *hQCD = get1DHistogram((hName+"QCD"+"Diff"+"Minus").c_str());

  TH1F *hTTbarS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str());
  TH1F *hDYJetsS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str()); // I need to copy unexisting histogram
  TH1F *hSoupS = get1DHistogram((hName+"Data"+"Diff"+"Plus").c_str());
  TH1F *hWJetsS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str());
  TH1F *hQCDS = get1DHistogram((hName+"QCD"+"Diff"+"Minus").c_str());

  MakeDiff(hTTbar, hDYJets, hSoup, hWJets, hQCD, hTTbarS, hDYJetsS, hSoupS, hWJetsS, hQCDS, varName, selName,"Plus");

  int rebinFactor = 1;  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hTTbar->Rebin(rebinFactor);
  hDYJets->Rebin(rebinFactor);
  hQCD -> Rebin(rebinFactor);

  hSoupS->Rebin(rebinFactor);
  hWJetsS->Rebin(rebinFactor);
  hTTbarS->Rebin(rebinFactor);
  hDYJetsS->Rebin(rebinFactor);
  hQCDS  -> Rebin(rebinFactor);

  float lumi = getLumi();

  TH1F *test = (TH1F*)hSoup->Clone(hName.c_str());
  test  -> Reset();
  test  -> Add(hSoup,1);

// Delta from the data
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);
  hSoup->Add(hQCD,-1);

  hName = hName + selName + "WEstimation";

  TH1F *W = (TH1F*)hSoup->Clone(hName.c_str());
  W  -> Reset();
  W  -> Add(hSoup,1);

  std::string hNameMC = hName +"MC" ;
  TH1F* WMC  = (TH1F*)hSoup->Clone(hNameMC.c_str());
  WMC  -> Reset();
  WMC  -> Add(hWJets,1);

// Asymmetry from Data

  hSoupS->Add(hDYJetsS,-1);
  hSoupS->Add(hTTbarS,-1);
  hSoupS->Add(hQCDS,-1);

  hSoup->Divide(hSoupS);

// Asymmetry from MC
  hWJets -> Divide(hWJetsS);

// WJet background estimated from the data
  W -> Divide(hWJets);

// WJet background estimated from the MC
  WMC -> Divide(hWJets);

// drawing
  
  hNameMC = "AsymmMC_"+varName+"_"+ selName;
  PlotOneHistogram(varName,hNameMC,hWJets,"Asymm MC");

// return WJet background

  std::pair <TH1*, TH1*> WEstimation;
  WEstimation = make_pair(W,WMC);

  return std::make_pair( WEstimation, hWJets);
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
THStack*  HTTWeightHistograms::PlotDiff(std::string varName, std::string selName, std::string SubSelName){

  std::string hName = "h1D"+varName;
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());
  TH1F *hDYJets = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());   
  TH1F *hSoup = get1DHistogram((hName+"Data"+"Diff"+"Minus").c_str());
  TH1F *hWJets = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());
  TH1F *hQCD = get1DHistogram((hName+"QCD"+"Diff"+"Minus").c_str());

  TH1F *hTTbarS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str());
  TH1F *hDYJetsS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str()); // I need to copy unexisting histogram
  TH1F *hSoupS = get1DHistogram((hName+"Data"+"Diff"+"Plus").c_str());
  TH1F *hWJetsS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str());
  TH1F *hQCDS = get1DHistogram((hName+"QCD"+"Diff"+"Minus").c_str());

  MakeDiff(hTTbar, hDYJets, hSoup, hWJets, hQCD, hTTbarS, hDYJetsS, hSoupS, hWJetsS, hQCDS, varName, selName,"Plus");

  float lumi = getLumi();

  std::string selName2="";

  if(selName.find("All")!=std::string::npos) {
	selName="";
	selName2="All";
	}

// corrected W+Jet background
  hWJets -> Reset();
  TH1F *hWJets11 = get1D_WJet_Histogram((hName+"WJets"+selName+"Plus").c_str());
  TH1F *hWJets21 = get1D_WJet_Histogram((hName+"WJets"+selName+"Minus").c_str());

  std::string sampleName = "WJets";
  std::string WselType = "wselOS";
  if(selName.find("SS")!=std::string::npos) WselType = "wselSS";
  pair<float,float> dataToMCScaleP = getWNormalisation( WselType,"Plus");
  pair<float,float> dataToMCScaleM = getWNormalisation( WselType,"Minus");
  float weight = getSampleNormalisation(sampleName);
  float scale=1;
        scale = weight*lumi*dataToMCScaleP.first;
  	hWJets11->Scale(scale);
        scale = weight*lumi*dataToMCScaleM.first;
  	hWJets21->Scale(scale);

  if(selName2=="All"){
  TH1F *hWJets12 =  get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+"Plus").c_str());
  TH1F *hWJets22 =  get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+"Minus").c_str());
       WselType = "wselSS";
    std::pair<float,float> dataToMCScaleP = getWNormalisation(WselType,"Plus");
    std::pair<float,float> dataTOMCScaleM = getWNormalisation(WselType,"Minus");
       scale = weight*lumi*dataToMCScaleP.first;
    hWJets11-> Scale(scale);
    hWJets11 -> Add(hWJets12,1);
       scale = weight*lumi*dataTOMCScaleM.first;
    hWJets21-> Scale(scale);
    hWJets21 -> Add(hWJets22,1);
  }

  hWJets -> Add(hWJets11,1);
  hWJets -> Add(hWJets21,-1);

// plot stack

  int rebinFactor = 1;  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hTTbar->Rebin(rebinFactor);
  hDYJets->Rebin(rebinFactor);
  hQCD -> Rebin(rebinFactor);

  hSoupS->Rebin(rebinFactor);
  hWJetsS->Rebin(rebinFactor);
  hTTbarS->Rebin(rebinFactor);
  hDYJetsS->Rebin(rebinFactor);
  hQCDS  -> Rebin(rebinFactor);

// plot stack
  hSoup->SetLineColor(1);
  hSoup->SetFillColor(1);
  hSoup->SetMarkerStyle(20);

  hWJets->SetFillColor(kRed+2);
  hTTbar->SetFillColor(kBlue+2);
  hDYJets->SetFillColor(kOrange-4);
  hQCD->SetFillColor(kMagenta-10);

  hSoup->SetLineWidth(1);

  THStack *hs = new THStack("hs","Stacked histograms");      
  /////////
  hs->Add(hWJets,"hist");
  hs->Add(hQCD,"hist");
  hs->Add(hDYJets,"hist");
  hs->Add(hTTbar,"hist");

  ////////
  TH1F *hMCSum = (TH1F*)hDYJets->Clone("hMCSum");
  hMCSum->Reset();
  hMCSum->Add(hDYJets);
  hMCSum->Add(hWJets);
  hMCSum->Add(hTTbar);
  hMCSum->Add(hQCD);

  if(!selName.size()) selName = "baseline";
  cout<<"Event count summary for selecion name: "<<selName<<std::endl;
  std::cout<<"Data: "<<hSoup->Integral(0,hSoup->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC: "<<hMCSum->Integral(0,hMCSum->GetNbinsX()+1)<<std::endl;  
  std::cout<<"MC W->l: "<<hWJets->Integral(0,hWJets->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC TTbar: "<<hTTbar->Integral(0,hTTbar->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->ll: "<<hDYJets->Integral(0,hDYJets->GetNbinsX()+1)<<std::endl;  
  std::cout<<"QCD: "<<hQCD->Integral(0,hQCD->GetNbinsX()+1)<<std::endl; 
  std::cout<<"Correction factors:"<<std::endl;
  std::cout<<"----------------------------------------"<<std::endl;

  TCanvas *c1 = getDefaultCanvas();
  c1->SetName("c1");
  c1->SetTitle("HTauTau analysis");
  c1->Divide(2);

  TPad *pad1 = (TPad*)c1->GetPad(1);
  TPad *pad2 = (TPad*)c1->GetPad(2);
  pad1->SetPad(0.01,0.29,0.99,0.99);
  pad2->SetPad(0.01,0.01,0.99,0.29);
  pad1->SetRightMargin(0.23);
  pad2->SetRightMargin(0.23);
  pad2->SetFillStyle(4000);
  ///
  pad1->Draw();
  pad1->cd();

  if(!selName.size()) selName = "baseline";
  hs->SetMaximum(4400);
  hs->Draw("hist");
  hs -> SetTitle("");
  hs->GetXaxis()->SetTitle(varName.c_str());
  hs->GetYaxis()->SetTitleOffset(1.4);
  hMCSum->SetFillColor(5);
  /////////
  float highEnd = 150;
  float lowEnd = -150;

  int binHigh = hs->GetXaxis()->FindBin(highEnd);  
  int binLow = hs->GetXaxis()->FindBin(lowEnd);

  if(hs->GetXaxis()->GetXmax()<highEnd) binHigh = hs->GetXaxis()->GetNbins();
  if(hs->GetXaxis()->GetXmin()>lowEnd) lowEnd = 1;

  hs->GetXaxis()->SetRange(binLow,binHigh);

  char yTitle[200];
 // sprintf(yTitle,"Events/%2.1f",hSoup->GetXaxis()->GetBinWidth(1));
  hs->GetYaxis()->SetTitle("Number of events");
  
  float max = hs->GetMaximum();
  if(hSoup->GetMaximum()>max) max = hSoup->GetMaximum();

  hs->GetHistogram()->SetTitleOffset(1.0);
  hs->SetMaximum(1.3*max);
  hs->SetMinimum(0.0);

  hSoup->Draw("same");
  TH1F *hMCSumref = (TH1F*)hMCSum->Clone("hMCSumref");
  hMCSumref->SetLineColor(12);
//  hMCSumref->Draw("same");

  TH1F *hEmpty = new TH1F("hEmpty","",1,0,1);
  hEmpty->SetLineColor(10);
  hEmpty->SetFillColor(10);

  TLegend *leg = new TLegend(0.79,0.32,0.99,0.82,NULL,"brNDC");
  setupLegend(leg);
  leg->AddEntry(hSoup,"Data","lep");
  leg->AddEntry(hDYJets,"Z#rightarrow ll","f");
  leg->AddEntry(hWJets,"W#rightarrow l #nu","f");
  leg->AddEntry(hTTbar,"TTbar","f");
  leg->AddEntry(hQCD,"QCD","f");
  leg->SetHeader(Form("#int L = %.2f fb^{-1}",lumi/1000));
  leg->Draw();

  float x = 0.6*(hs->GetXaxis()->GetXmax() - 
		 hs->GetXaxis()->GetXmin()) +
    hs->GetXaxis()->GetXmin(); 

  float y = 0.8*(max - 
		 hs->GetMinimum()) +
                 hs->GetMinimum(); 
  c1->cd();
  pad2->Draw();
  pad2->cd();

  hMCSum->GetXaxis()->SetRange(binLow,binHigh);
  hMCSum->SetTitle("");
  hMCSum->SetXTitle("");
  hMCSum->SetYTitle("#frac{N_{obs} - N_{exp}}{#sqrt{N_{obs}}}");
  hMCSum->GetXaxis()->SetLabelSize(0.09);
  hMCSum->GetYaxis()->SetLabelSize(0.09);
  hMCSum->GetYaxis()->SetTitleSize(0.09);
  hMCSum->GetYaxis()->SetTitleOffset(0.5);
  hMCSum->Add(hSoup,-1);
  for(int i=0;i<hMCSum->GetNbinsX()+1;++i){
    if(hSoup->GetBinContent(i)>0) hMCSum->SetBinContent(i,-hMCSum->GetBinContent(i)/sqrt(hSoup->GetBinContent(i)));
    else  hMCSum->SetBinContent(i,0);
    hMCSum->SetBinError(i,0);
  }
  hMCSum->SetLineWidth(3);
  hMCSum->SetMinimum(-5);
  hMCSum->SetMaximum(5);
  hMCSum->SetStats(kFALSE);
  hMCSum->Draw("hist");
  TLine *aLine = new TLine(hMCSum->GetXaxis()->GetXmin(),0.0,highEnd,0.0);
  aLine->SetLineColor(1);
  aLine->SetLineWidth(2);
  aLine->Draw();

  string plotName;
  selName="Diff_Stack_"+varName+"_"+selName+selName2+"_"+SubSelName;
  hName=hName+SubSelName;
  if(hName.find_last_of("/")<string::npos) plotName = "fig_pngA/" + hName.substr(hName.find_last_of("/")) + ".png";    
  else plotName = "fig_pngA/hTree_"+hName+Form("_%s",selName.c_str())+".png";
  c1->Print(plotName.c_str());

  if(hName.find_last_of("/")<string::npos) plotName = "fig_C/" + hName.substr(hName.find_last_of("/")) + ".C";    
  else plotName = "fig_C/hTree_"+hName+Form("_%s",selName.c_str())+".C";
  c1->Print(plotName.c_str()); 

  pad1->SetLogy(1);
  if(hName.find_last_of("/")<string::npos) plotName = "fig_pngA/" + hName.substr(hName.find_last_of("/")) + "_LogY.png";    
  else plotName = "fig_pngA/hTree_"+hName+Form("_%s",selName.c_str())+"_LogY.png";
  c1->Print(plotName.c_str()); 

  std::cout<<"-------------------------------------------------------------"<<std::endl;

  return hs;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double * HTTWeightHistograms::WJetEstimation(std::string varName, std::string selName, std::string SubSelName){

  std::string hName = "h1D"+varName;

  std::pair<std::pair<TH1*,TH1*>, TH1*> Asymm = PlotAsymm(varName,selName,"Plus");
  std::pair<TH1*, TH1*> WEstimation = Asymm.first;

  TH1F * hWJetsAsymm = (TH1F*) WEstimation.first;
  TH1F * hWJetsAsymmMC = (TH1F*) WEstimation.second;

  std::string selName2="";
  if(selName.find("All")!=std::string::npos) {
	selName="";
	selName2="All";
	}

  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+selName).c_str());
  TH1F *hWJets2 = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS").c_str());

  TH1F *hWJetsMC = get1D_WJet_Histogram((hName+"WJets"+selName).c_str());
  TH1F *hWJets2MC = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS").c_str());

  int rebinFactor = 1;  
  hWJets->Rebin(rebinFactor); 
  hWJets2->Rebin(rebinFactor);

// different scalling for OS and SS WJets
  float lumi = getLumi();

  std::string WselType = "wselOS";
   if(selName.find("SS")!=std::string::npos) WselType = "wselSS";
   pair<float,float> dataToMCScale = getWNormalisation(WselType,"");

   TH1F* Werror  = (TH1F*)hWJets->Clone("Werror");
   Werror -> Reset();
  
 	for(int i=0; i< Werror->GetNbinsX()+2; i++){
 	Werror -> SetBinContent(i,dataToMCScale.first);
 	Werror -> SetBinError(i,dataToMCScale.second);
 	}

   float weight = getSampleNormalisation("WJets");
   float scale = weight*lumi;
   hWJets -> Scale(scale);
   hWJetsMC -> Scale(scale);
   hWJets -> Multiply(Werror);

   if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
    	WselType = "wselSS";
    	dataToMCScale = getWNormalisation(WselType,"");
 	hWJets2 -> Scale(scale);
	hWJets2MC -> Scale(scale);

	Werror->Reset();
 	for(int i=0; i< Werror->GetNbinsX()+2; i++){
 	Werror -> SetBinContent(i,dataToMCScale.first);
 	Werror -> SetBinError(i,dataToMCScale.second);
 	}

	hWJets2 -> Multiply(Werror);
     	hWJets -> Add(hWJets2,1);
	hWJetsMC -> Add(hWJets2MC,1);
   }

// calculate contribution to the background for mT<40GeV; classical and asymmetry

  double  StatError40;
  double  StatError40MC;
  if(varName.find("MassTrans")!=std::string::npos){
	//int iBin = hWJetsAsymm->FindBin(40);
	double intWJet = hWJets -> IntegralAndError(0,10,StatError40MC,"");
	double intWJetAsymm = hWJetsAsymm -> IntegralAndError(0,10,StatError40,"");
	cout<<"10 binów"<<varName<<"_"<<selName<<"_"<<selName2<<"_"<<SubSelName<<"_"<<": intWJet "<<intWJet<<"pm"<<StatError40MC<<" intWJetAsymm "<<intWJetAsymm<<"pm"<<StatError40<<std::endl;
  }

  if(varName.find("Eta")!=std::string::npos){
	//int iBin = hWJetsAsymm->FindBin(40);
	double intWJet = hWJets -> IntegralAndError(0,hWJets->GetNbinsX()+1,StatError40MC,"");
	double intWJetAsymm = hWJetsAsymm -> IntegralAndError(0,hWJetsAsymm->GetNbinsX()+1,StatError40,"");
	cout<<"10 binów"<<varName<<"_"<<selName<<"_"<<selName2<<"_"<<SubSelName<<"_"<<": intWJet "<<intWJet<<"pm"<<StatError40MC<<" intWJetAsymm "<<intWJetAsymm<<"pm"<<StatError40<<std::endl;
  }

// drawing
  hName="WJetBack_"+varName+"_"+selName+"_"+selName2;
  PlotTwoHistograms(varName,hName, hWJetsAsymm,"Asymm Data", hWJets, "MC scaling");

  hName="WJetBackMC_"+varName+"_"+selName+"_"+selName2;
  PlotTwoHistograms(varName,hName, hWJetsAsymmMC,"Asymm MC", hWJetsMC, "MC scaling");

return 0;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTWeightHistograms::getWNormalisation(std::string selName, std::string SubSelName){

  std::cout<<"Calling method: "<<__func__<<selName<<"_"<<SubSelName<<std::endl;

  std::string hName = "h1DMassTrans";
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+selName+SubSelName).c_str());
  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+selName+SubSelName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+selName+SubSelName).c_str());
 
  TH1F *hTT = get1DHistogram((hName+"TTbar"+selName+SubSelName).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+selName+SubSelName).c_str());

  float lumi = getLumi();

  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hWJets->Clone((hName+"hDYJetsLowM"+selName).c_str()); hDYJetsLowM->Reset();
  }
  if(!hTT){
    hTT = (TH1F*)hWJets->Clone((hName+"hTTbar"+selName+SubSelName).c_str()); hTT->Reset();
  } 
		 
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "WJets";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hWJets->Scale(scale);

  sampleName = "DYJetsLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowM->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTT->Scale(scale);

  // Prepare a histogram with data minus backgrounds: DYJets, hTT, Other
  TH1F* datamtlo = (TH1F*)hSoup->Clone("datamtlo");
  datamtlo->Add(hDYJets,-1);
  datamtlo->Add(hDYJetsLowM,-1);
  datamtlo->Add(hTT,-1);

  float inthWJets=hWJets->Integral(0,hWJets->GetNbinsX()+1);
  float intdata=datamtlo->Integral(0,datamtlo->GetNbinsX()+1);

  // Calculate weight
  weight=intdata/inthWJets;

  float dweight;
  float inthSoup = hSoup->Integral(0,hSoup->GetNbinsX()+1);
  float inthDYJetsLowM = hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1);
  float inthDYJets = hDYJets->Integral(0,hDYJets->GetNbinsX()+1);
  float inthTT = hTT->Integral(0,hTT->GetNbinsX()+1);
  float inthOther = 0;//hOther->Integral(0,hOther->GetNbinsX()+1);
  dweight=((inthSoup+inthDYJets+inthDYJetsLowM+inthTT+inthOther)/inthWJets/inthWJets+intdata*intdata/(inthWJets*inthWJets*inthWJets));
  dweight=sqrt(dweight);
  cout<<"Selecion name: "<<selName<<std::endl;
  cout<<"DATA: "<<inthSoup<<" DATA - MC(!WJets): "<<intdata<<" MC WJets "<<inthWJets
      <<" DYJets: "<<inthDYJets<<" DYJetsLowM: "<<inthDYJetsLowM
      <<" TTbar: "<<inthTT<<" Other: "<<inthOther<<endl;
  cout<<"WJets scale:"<<weight<<" dweight "<<dweight<<endl;
  return std::make_pair(weight, dweight);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::PlotTwoHistograms(std::string varName,std::string nazwa, TH1* th1,std::string opis1, TH1* th2, std::string opis2){

  TCanvas* c = new TCanvas("WEstim","Asymm",460,500);
  TPad* p1 = new TPad("p1","p1",0.06,0.03,0.98,0.98,0); p1->Draw();
  p1->SetLeftMargin(0.15);
  p1->Draw();
  p1->cd();
 
  th1->Draw("hist");
  th1->SetStats(kFALSE);
  th2->SetLineColor(2);
  th2->Draw("same");

  th1->SetFillColor(kWhite);

  float max = th1->GetMaximum();
  th1->SetMaximum(1.3*max);

  float th1int = th1->Integral(0,th1->GetNbinsX()+1);

  if(varName.find("Eta")!=std::string::npos && (nazwa.find("SS")!=std::string::npos || nazwa.find("All")!=std::string::npos || nazwa.find("MC")!=std::string::npos)){

  	TLegend *leg2 = new TLegend(0.5,0.12,0.7,0.32,NULL,"brNDC");
  	setupLegend(leg2);
  	leg2->SetTextSize(0.03);
  	leg2->AddEntry(th1,opis1.c_str(),"l");
  	leg2->AddEntry(th2,opis2.c_str(),"l");
  	leg2->SetHeader(Form("#int N = %.2f",th1int));
	leg2->Draw();
  }  else {
  	TLegend *leg = new TLegend(0.6,0.6,0.90,0.88,NULL,"brNDC");
  	setupLegend(leg);
  	leg->SetTextSize(0.03);
  	leg->AddEntry(th1,opis1.c_str(),"l");
  	leg->AddEntry(th2,opis2.c_str(),"l");
  	leg->SetHeader(Form("#int N = %.2f",th1int));
  	leg->Draw();
  }

  th1->GetYaxis()->SetTitle("Number of events");
  th2->GetYaxis()->SetTitle("Number of events");
  th1->GetYaxis()->SetTitleOffset(2.3);

  if(varName.find("SVfit")!=std::string::npos){
    th1->GetXaxis()->SetTitle("SVFit mass [GeV]");
    th2->GetXaxis()->SetTitle("SVFit mass [GeV]");
    }

  if(varName.find("MassTrans")!=std::string::npos){
    th1->GetXaxis()->SetTitle("M_{T} [GeV]");
    th2->GetXaxis()->SetTitle("M_{T} [GeV]");
    }

  if(varName.find("NPV")!=std::string::npos){
    th1->GetXaxis()->SetTitle("NPV");
    th2->GetXaxis()->SetTitle("NPV");
    }

  if(varName.find("Eta")!=std::string::npos){
    th1->GetXaxis()->SetTitle("#eta [units]");
    th2->GetXaxis()->SetTitle("#eta [units]");
    }

  if(varName.find("MassVis")!=std::string::npos){
    th1->GetXaxis()->SetTitle("M_{VIS} [GeV]");
    th2->GetXaxis()->SetTitle("M_{VIS} [GeV]");
    }

  c->Print(TString::Format("fig_pngA/%s.png",nazwa.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::PlotOneHistogram(std::string varName,std::string nazwa, TH1* th1,std::string opis1){

  TCanvas* c = new TCanvas("WEstim","Asymm",460,500);
  TPad* p1 = new TPad("p1","p1",0.06,0.03,0.98,0.98,0); p1->Draw();
  p1->SetLeftMargin(0.1);
  p1->Draw();
  p1->cd();

  float max = th1->GetMaximum();
  th1->SetMaximum(1.3*max);
 
  th1->Draw("hist");
  th1->SetStats(kFALSE);

  TLegend *leg = new TLegend(0.4,0.8,0.65,0.88,NULL,"brNDC");
  setupLegend(leg);
  leg->SetTextSize(0.03);
  leg->AddEntry(th1,opis1.c_str(),"l");
  leg->Draw();

  th1->GetYaxis()->SetTitle("Number of events");
  th1->GetYaxis()->SetTitleOffset(1.7);

  if(varName.find("SVfit")!=std::string::npos){
    th1->GetXaxis()->SetTitle("SVFit mass [GeV]");
    }

  if(varName.find("MassTrans")!=std::string::npos){
    th1->GetXaxis()->SetTitle("M_{T} [GeV]");
    }

  if(varName.find("NPV")!=std::string::npos){
    th1->GetXaxis()->SetTitle("NPV");
    }

  if(varName.find("Eta")!=std::string::npos){
    th1->GetXaxis()->SetTitle("#eta [units]");
    }

  if(varName.find("MassVis")!=std::string::npos){
    th1->GetXaxis()->SetTitle("M_{VIS} [GeV]");
    }

  c->Print(TString::Format("fig_pngA/%s.png",nazwa.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTWeightHistograms::getQCDOStoSS(std::string selName, std::string SubSelName){

  std::cout<<"Calling method: "<<__func__<<std::endl;
  if(selName.find("SS")!=std::string::npos) return  std::make_pair(1.0,0.0);

  std::string hName = "h1DIso";

  // SS selection
  TH1F *hWJetsSS = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hDYJetsLowMSS = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hDYJetsSS = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hTTbarSS = get1DHistogram((hName+"TTbar"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hSoupSS = get1DHistogram((hName+"Data"+"qcdselSS"+selName+SubSelName).c_str());

  if(!hDYJetsLowMSS){
    hDYJetsLowMSS = (TH1F*)hSoupSS->Clone((hName+"DYJetsLowM"+"qcdselSS"+selName+SubSelName).c_str()); hDYJetsLowMSS->Reset();
  } 
  if(!hTTbarSS){
    hTTbarSS = (TH1F*)hWJetsSS->Clone((hName+"hTTbar"+"qcdselSS"+selName+SubSelName).c_str()); hTTbarSS->Reset();
  } 

  // OS selection
  TH1F *hWJetsOS = get1D_WJet_Histogram((hName+"WJets"+"qcdselOS"+selName+SubSelName).c_str());
  TH1F *hDYJetsLowMOS = get1DHistogram((hName+"DYJetsLowM"+"qcdselOS"+selName+SubSelName).c_str());
  TH1F *hDYJetsOS = get1D_DY_Histogram((hName+"DYJets"+"qcdselOS"+selName+SubSelName).c_str());
  TH1F *hTTbarOS = get1DHistogram((hName+"TTbar"+"qcdselOS"+selName+SubSelName).c_str());
  TH1F *hSoupOS = get1DHistogram((hName+"Data"+"qcdselOS"+selName+SubSelName).c_str());

  if(!hDYJetsLowMOS){
    hDYJetsLowMOS = (TH1F*)hSoupOS->Clone((hName+"DYJetsLowM"+"qcdselOS"+selName+SubSelName).c_str()); hDYJetsLowMOS->Reset();
  }
  if(!hTTbarOS){
    hTTbarOS = (TH1F*)hWJetsOS->Clone((hName+"hTTbar"+"qcdselOS"+selName+SubSelName).c_str()); hTTbarOS->Reset();
  } 

 float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DYJetsLowM";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJetsLowMOS->Scale(scale);
  hDYJetsLowMSS->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsOS->Scale(scale);
  hDYJetsSS->Scale(scale);
  
  sampleName = "WJets";
  weight = getSampleNormalisation(sampleName);
  scale= weight*lumi;
  	hWJetsOS->Scale(getWNormalisation("wselOS",SubSelName).first*scale);
  	hWJetsSS->Scale(getWNormalisation("wselSS",SubSelName).first*scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbarOS->Scale(scale);
  hTTbarSS->Scale(scale);
 
  ///Subtract backgrounds (not QCD)
  hSoupSS->Add(hWJetsSS,-1);
  hSoupSS->Add(hDYJetsSS,-1);  
  hSoupSS->Add(hDYJetsLowMSS,-1);
  hSoupSS->Add(hTTbarSS,-1);
  
  hSoupOS->Add(hWJetsOS,-1);
  hSoupOS->Add(hDYJetsOS,-1);  
  hSoupOS->Add(hDYJetsLowMOS,-1);
  hSoupOS->Add(hTTbarOS,-1);
  
  hSoupOS->Divide(hSoupSS);

  //funtion fitting
  TF1 *line=new TF1("line","[0]",0,2);
  line->SetParameter(0,1);
  TCanvas* c = new TCanvas("QCD_OStoSS","QCD_OStoSS",460,500);
    TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->Clear();
   leg->Draw();

  hSoupOS->SetLineWidth(3);
  hSoupOS->GetYaxis()->SetTitleOffset(1.4);
  hSoupOS->GetYaxis()->SetTitle("OS/SS");
  hSoupOS->GetXaxis()->SetTitle("muon relative isolation");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(11);
  hSoupOS->SetStats(kFALSE);
  hSoupOS->Draw();
  hSoupOS->Fit("line","","",0.3,0.5);

  hName = hName + "_" + selName +"_"+ SubSelName;
  c->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());

  float param, dparam;
  param=line->GetParameter(0);
  dparam=line->GetParError(0);

  std::cout<<"QCD OS/SS ratio: "<<param<<" +- "<<dparam<<std::endl;

  return std::make_pair(param, dparam);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTWeightHistograms::getQCDbackground(std::string varName, std::string selName, std::string SubSelName){

  std::cout<<"Calling method: "<<__func__<<std::endl;
  
  ///Not very clear and elegant. AK
  ///Need this to avoid resursive control region labels like
  ///qcdselSSqcdselOS
  int ifscale=0; //when SS there should be no scaling
  std::string OS="";
  if(selName.find("SS")==std::string::npos) {ifscale=1; OS="OS";}
  if(selName.find("qcdsel")!=std::string::npos) selName = "";

  std::string hName = "h1D" + varName;
  // SS selection
  TH1F *hWJetsSS = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hDYJetsLowMSS = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hDYJetsSS = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hTTbarSS = get1DHistogram((hName+"TTbar"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hSoupSS = get1DHistogram((hName+"Data"+"qcdselSS"+selName+SubSelName).c_str());

  if(!hSoupSS) return 0;

  if(!hDYJetsLowMSS){
    hDYJetsLowMSS = (TH1F*)hSoupSS->Clone((hName+"DYJetsLowM"+"qcdselSS"+selName+SubSelName).c_str()); hDYJetsLowMSS->Reset();
  } 
  if(!hTTbarSS){
    hTTbarSS = (TH1F*)hSoupSS->Clone((hName+"hTTbar"+"qcdselSS"+selName+SubSelName).c_str()); hTTbarSS->Reset();
  }

  //////////////////////////////////////////////////////////////////////
 float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DYJetsLowM";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJetsLowMSS->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsSS->Scale(scale);
  
  sampleName = "WJets";
  weight = getSampleNormalisation(sampleName);
  scale= weight*lumi;
  	hWJetsSS->Scale(getWNormalisation("wselSS",SubSelName).first*scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbarSS->Scale(scale);

  hSoupSS->SetName(("h1DQCDEstimate"+varName+SubSelName).c_str());

  hSoupSS->Add(hWJetsSS,-1);
  hSoupSS->Add(hDYJetsSS,-1);
  hSoupSS->Add(hDYJetsLowMSS,-1);
  hSoupSS->Add(hTTbarSS,-1);

  scale =1;
  if(ifscale==1){
  scale = getQCDOStoSS(selName,SubSelName).first;
  hSoupSS->Scale(scale);
  }

  ///Clean up the QCD shape, and remove fluctuations around 0 counts.
//  for(unsigned int iBinX=0;iBinX<=hSoupSS->GetNbinsX();++iBinX){
//    if(hSoupSS->GetBinContent(iBinX)<3.0) hSoupSS->SetBinContent(iBinX,0);
//  }

	stringstream ss;
	ss << scale;
	string str = ss.str();

  std::string name="QCDBack_"+varName+"_"+OS+"_"+selName+"_"+SubSelName+"_scale="+str;
  PlotOneHistogram(varName,name, hSoupSS,"QCD back");

  return hSoupSS;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
