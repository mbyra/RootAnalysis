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
  TH1F *hStats = get1DHistogram(hName.c_str());
  if(sampleName=="DYJets") hStats = get1D_DY_Histogram(hName.c_str());
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

  if(sampleName=="WJetsHT100to200"){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2#Backgrounds_DY_jets_AN2
    crossSection = 1345;
  }

  if(sampleName=="WJetsHT200to400"){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2#Backgrounds_DY_jets_AN2
    crossSection = 359.7;
  }

  if(sampleName=="WJetsHT400to600"){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2#Backgrounds_DY_jets_AN2
    crossSection = 48.91;
  }

  if(sampleName=="WJetsHT600toInf"){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2#Backgrounds_DY_jets_AN2
    crossSection = 18.77;
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
  TH1F *hWJets = get1DHistogram(hName.Data());
  return hWJets; //TEST

  if(!hWJets) return hWJets;

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

  if(hWJets200to400) hWJets->Add(hWJets100to200);
  if(hWJets400to600) hWJets->Add(hWJets200to400);
  if(hWJets600toInf) hWJets->Add(hWJets400to600);
  hWJets->SetName(name.c_str());

  return hWJets;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

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
   add1DHistogram("h1DMassTemplate",";SVFit mass [GeV/c^{2}]; Events",25,0,200,file_); //temporary change: 25 instead of 50
   add1DHistogram("h1DEtaTemplate",";#eta; Events",23,-2.3,2.3,file_);
   float bins[31] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5};
   add1DHistogram("h1DIsoTemplate",";Isolation; Events",30,bins,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

//  PlotAsymm("EtaMuon","","Plus"); 

  TFile f("AsymmWeights.root","recreate");

  std::pair<std::pair<TH1*,TH1*>, TH1*> WEstimation = PlotAsymm("EtaMuon","","Plus");
  TH1F * hAsymmEtaMC = (TH1F*) WEstimation.second;
  hAsymmEtaMC -> SetName("EtaAsymmetryOS");
  hAsymmEtaMC -> Write(); 

  WEstimation = PlotAsymm("EtaMuon","qcdselSS","Plus");
  TH1F * hAsymmEtaMCSS = (TH1F*) WEstimation.second;
  hAsymmEtaMCSS -> SetName("EtaAsymmetrySS");
  hAsymmEtaMCSS -> Write(); 

  WEstimation = PlotAsymm("MassTrans","","Plus");
  TH1F * hAsymmMassTransMC = (TH1F*) WEstimation.second;
  hAsymmMassTransMC -> SetName("MassTransAsymmetryOS");
  hAsymmMassTransMC -> Write(); 

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

  WJetEstimation("EtaMuon","","Plus");
  WJetEstimation("EtaMuon","All","Plus");
  WJetEstimation("EtaMuon","qcdselSS","Plus");
  WJetEstimation("MassVis","All","Plus");

  std::pair<float,float> ErrorOS = WJetEstimation("MassTrans","","Plus");
  cout<<"OS: StatError40 "<<ErrorOS.first<<" StatError40MC "<<ErrorOS.second<<std::endl;

  std::pair<float,float> ErrorAll = WJetEstimation("MassTrans","All","Plus");
  cout<<"All: StatError40 "<<ErrorAll.first<<" StatError40MC "<<ErrorAll.second<<std::endl;

  std::pair<float,float> ErrorSS = WJetEstimation("MassTrans","qcdselSS","Plus");
  cout<<"SS: StatError40 "<<ErrorSS.first<<" StatError40MC "<<ErrorSS.second<<std::endl;

  f.Close();

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double HTTWeightHistograms::MakeDiff(TH1F *hTTbar, TH1F* hDYJets, TH1F* hSoup, TH1F* hWJets, TH1F* hQCD, TH1F *hTTbarS, TH1F* hDYJetsS, TH1F* hSoupS, TH1F* hWJetsS, TH1F* hQCDS, std::string varName, std::string selName, std::string SubSelName){

 //note: I left QCD just in case. If we will find "a brand new" method for QCD it will be good to subtract it also from Delta in asymmetry.
  std::cout<<"--- Drawing THStack for variable: "<<varName
	   <<" selection: "<<selName<<std::endl;

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
  TH1F *hWJets11 = get1DHistogram((hName+"WJets"+selName+"Plus").c_str());
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
  TH1F *hWJets21 = get1DHistogram((hName+"WJets"+selName+"Minus").c_str());
  TH1F *hDYJets2 = get1D_DY_Histogram((hName+"DYJets"+selName+"Minus").c_str());
  TH1F *hDYJetsLowM2 = get1DHistogram((hName+"DYJetsLowM"+selName+"Minus").c_str());

  if(!hDYJets2){
    hDYJets2 = (TH1F*)hSoup1->Clone((hName+"hDYJets2"+selName).c_str()); hDYJets2->Reset();
  }
  if(!hDYJetsLowM2){
    hDYJetsLowM2 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM2"+selName).c_str()); hDYJetsLowM2->Reset();
  }

// sum OS and SS histograms for "All" asymmetry

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
  TH1F *hWJets12 = get1DHistogram((hName+"WJets"+"qcdselSS"+"Plus").c_str());
  TH1F *hWJets22 = get1DHistogram((hName+"WJets"+"qcdselSS"+"Minus").c_str());
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

// plot stack diff for Plus / Minus

  hSoup1->SetLineColor(1);
  hSoup1->SetFillColor(1);
  hSoup1->SetMarkerStyle(20);

  hWJets11->SetFillColor(kRed+2);
  hTTbar1->SetFillColor(kBlue+2);
  hDYJets1->SetFillColor(kOrange-3);
  hDYJetsLowM1->SetFillColor(kOrange-7);

  THStack *hs1 = new THStack("hs1","Stacked histograms Plus");          
  hs1->Add(hTTbar1,"hist");
  hs1->Add(hWJets11,"hist");
  hs1->Add(hDYJetsLowM1,"hist");
  hs1->Add(hDYJets1,"hist");

  TCanvas* c = new TCanvas("Diff","Plus",460,500);
  std::string nazwa = "Stack_Plus_"+varName+"_" + selName +"_"+selName2;
  hs1->Draw("hist");
  hSoup1->Draw("same");

  c->Print(TString::Format("fig_pngA/%s.png",nazwa.c_str()).Data());

// plot stack diff for Minus

  hSoup2->SetLineColor(1);
  hSoup2->SetFillColor(1);
  hSoup2->SetMarkerStyle(20);

  hWJets21->SetFillColor(kRed+2);
  hTTbar2->SetFillColor(kBlue+2);
  hDYJets2->SetFillColor(kOrange-3);
  hDYJetsLowM2->SetFillColor(kOrange-7);

  THStack *hs2 = new THStack("hs1","Stacked histograms Plus");          
  hs2->Add(hTTbar2,"hist");
  hs2->Add(hWJets21,"hist");
  hs2->Add(hDYJetsLowM2,"hist");
  hs2->Add(hDYJets2,"hist");

  TCanvas* c2 = new TCanvas("Diff","Plus",460,500);
  nazwa = "Stack_Minus_"+varName+"_" + selName +"_"+selName2;
  hs2->Draw("hist");
  hSoup2->Draw("same");

  c2->Print(TString::Format("fig_pngA/%s.png",nazwa.c_str()).Data());

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

return 0;
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
std::pair<std::pair<TH1*,TH1*>,TH1*>  HTTWeightHistograms::PlotAsymm(std::string varName, std::string selName, std::string SubSelName){

  std::string hName = "h1D"+varName;
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());
  TH1F *hDYJets = get1DHistogram((hName+"TTbar"+"Diff"+"Minus").c_str());   
  TH1F *hSoup = get1DHistogram((hName+"Data"+"Diff"+"Minus").c_str());
  TH1F *hWJets = get1DHistogram((hName+"WJets"+"Diff"+"Minus").c_str());
  TH1F *hQCD = get1DHistogram((hName+"QCD"+"Diff"+"Minus").c_str());

  TH1F *hTTbarS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str());
  TH1F *hDYJetsS = get1DHistogram((hName+"TTbar"+"Diff"+"Plus").c_str()); // I need to copy unexisting histogram
  TH1F *hSoupS = get1DHistogram((hName+"Data"+"Diff"+"Plus").c_str());
  TH1F *hWJetsS = get1DHistogram((hName+"WJets"+"Diff"+"Plus").c_str());
  TH1F *hQCDS = get1DHistogram((hName+"QCD"+"Diff"+"Minus").c_str());

  MakeDiff(hTTbar, hDYJets, hSoup, hWJets, hQCD, hTTbarS, hDYJetsS, hSoupS, hWJetsS, hQCDS, varName, selName,"Plus");

  int rebinFactor = 1;  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hTTbar->Rebin(rebinFactor);
  hDYJets->Rebin(rebinFactor);

  hSoupS->Rebin(rebinFactor);
  hWJetsS->Rebin(rebinFactor);
  hTTbarS->Rebin(rebinFactor);
  hDYJetsS->Rebin(rebinFactor);

  float lumi = getLumi();

// Delta from the data
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);
// I delete it to have W+Jet estimation method without QCD  hSoup->Add(hQCD,-1);

//  TH1F* W  = (TH1F*)hSoup->Clone(("Asymm"+varName+selName).c_str);

  hName = hName + "_"+selName + "_WEstimation";

  TH1F *W = (TH1F*)hSoup->Clone(hName.c_str());
  W  -> Reset();
  W  -> Sumw2();
  W  -> Add(hSoup,1);

cout<<"hSoup Delta Dane Np-Nn: --------------------------------"<<endl;
W -> Print("all");
cout<<"hWJets Delta MC Np-Nn: --------------------------------"<<endl;
hWJets -> Print("all");

  std::string hNameMC = hName +"MC" ;
  TH1F* WMC  = (TH1F*)hSoup->Clone(hNameMC.c_str());
  WMC  -> Reset();
  WMC  -> Add(hWJets,1);

// here Asymmetry from Data

  hSoupS->Add(hDYJetsS,-1);
  hSoupS->Add(hTTbarS,-1);
//  hSoupS->Add(hQCDS,-1);

cout<<"hWJets MONTE CARLO: NMC--------------------------------"<<endl;
hWJetsS -> Print("all");
cout<<"hSoup DATA: ND--------------------------------"<<endl;
hSoupS -> Print("all");

  hSoup->Divide(hSoupS);

// here Asymmetry from MC
  hWJets -> Divide(hWJetsS);

cout<<"Asymetia MC: AMC--------------------------------"<<endl;
hWJets -> Print("all");

// WJet Bac. Estimation from the Data
  W -> Divide(hWJets);

// WJet Bac. Estimation from the MC
  WMC -> Divide(hWJets);

// drawing
  TCanvas *c2 = getDefaultCanvas();
  c2->SetName("c2");
  c2->SetTitle("W asymmetry");

  hName = "Asymm_"+ varName + "_"+selName;

  if(varName.find("SVfit")!=std::string::npos){
    hSoup->GetXaxis()->SetTitle("SVFit mass [GeV/c^{2}]");
    hSoup->SetTitle("asymmetry SVFit mass [GeV/c^{2}]");}

  if(varName.find("MassTrans")!=std::string::npos){
    hSoup->GetXaxis()->SetTitle("Mt [GeV/c^{2}]");
    hSoup->SetTitle("asymmetry Mt [GeV/c^{2}]");}

  if(varName.find("NPV")!=std::string::npos){
    hSoup->GetXaxis()->SetTitle("NPV");
    hSoup->SetTitle("asymmetry NPV");}

  if(varName.find("Eta")!=std::string::npos){
    hSoup->GetXaxis()->SetTitle("Eta");
    hSoup->SetTitle("asymmetry Eta");}

  hSoup->SetStats(kFALSE);

// drawing asymmetries
//  hSoup->Draw("hist");
  hWJets->SetLineColor(1);
  hWJets->SetLineWidth(1.2);
//  hWJets->SetMinimum(0.0);
//  hWJets->SetMaximum(0.5);
  hWJets->Draw("hist");
  
  c2->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());

// return WJetBack estimations

  std::pair <TH1*, TH1*> WEstimation;
  WEstimation = make_pair(W,WMC);

  return std::make_pair( WEstimation, hWJets);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<double,double> HTTWeightHistograms::WJetEstimation(std::string varName, std::string selName, std::string SubSelName){

  std::string hName = "h1D"+varName;

  std::pair<std::pair<TH1*,TH1*>, TH1*> Asymm = PlotAsymm(varName,selName,"Plus");
  std::pair<TH1*, TH1*> WEstimation = Asymm.first;

  TH1F * hWJetsAsymm = (TH1F*) WEstimation.first;
  TH1F * hWJetsAsymmMC = (TH1F*) WEstimation.second;

  std::string selName2="";
  if(selName=="All") {
	selName="";
	selName2="All";
	}

  TH1F *hWJets = get1DHistogram((hName+"WJets"+selName).c_str());
  hWJets->Sumw2();
  TH1F *hWJets2 = get1DHistogram((hName+"WJets"+"qcdselSS").c_str());
  hWJets2->Sumw2();

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
  hWJets->Scale(scale);
  hWJets-> Multiply(Werror);

  if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
    WselType = "wselSS";
    dataToMCScale = getWNormalisation(WselType,"");
    float scale = weight*lumi*dataToMCScale.first;
    hWJets2-> Scale(scale);
    hWJets -> Add(hWJets2,1);
  }

// drawing
  TCanvas *c2 = getDefaultCanvas();
  c2->SetName("c2");
  c2->SetTitle("WJets Background estimation");

  if(varName.find("SVfit")!=std::string::npos){
    hWJetsAsymm->GetXaxis()->SetTitle("SVFit mass [GeV/c^{2}]");
    hWJetsAsymm->SetTitle("WJet SVFit mass [GeV/c^{2}]");}

  if(varName.find("MassTrans")!=std::string::npos){
    hWJetsAsymm->GetXaxis()->SetTitle("Mt [GeV/c^{2}]");
    hWJetsAsymm->SetTitle("WJet Mt [GeV/c^{2}]");}

  if(varName.find("NPV")!=std::string::npos){
    hWJetsAsymm->GetXaxis()->SetTitle("NPV");
    hWJetsAsymm->SetTitle("WJet NPV");}

  if(varName.find("Eta")!=std::string::npos){
    hWJetsAsymm->GetXaxis()->SetTitle("Eta");
    hWJetsAsymm->SetTitle("WJet Eta");}

  hWJetsAsymm -> SetStats(kFALSE);

  hWJetsAsymm -> Draw("hist");

// estimation of statistical error

  double  StatError40;
  double  StatError40MC;
  if(varName.find("MassTrans")!=std::string::npos){
	//int iBin = hWJetsAsymm->FindBin(40);
	double intWJet = hWJets -> IntegralAndError(0,5,StatError40,"");
	double intWJetAsymm = hWJetsAsymm -> IntegralAndError(0,5,StatError40MC,"");
	cout<<"5 binów: intWJet "<<intWJet<<"pm"<<intWJet<<" intWJetAsymm "<<intWJetAsymm<<"pm"<<intWJetAsymm<<std::endl;
  }

//  hWJetsAsymm->SetMinimum(-2000);
//  hWJetsAsymm->SetMaximum(8000);
  hWJets -> SetLineColor(1);
  hWJets->SetLineWidth(1.2);
  hWJets->Draw("same");

 TLegend *leg = new TLegend(0.6,0.6,0.99,0.9,NULL,"brNDC");
  setupLegend(leg);
  leg->SetTextSize(0.03);
  leg->AddEntry(hWJets,"Old","l");
  leg->AddEntry(hWJetsAsymm,"Data Asymm","l");

  float inthWJetsAsymm = hWJetsAsymm->Integral(0,hWJetsAsymm->GetNbinsX()+1);
  leg->SetHeader(Form("#int N = %.2f",inthWJetsAsymm));
//  leg->AddEntry((TObject*)0, Form("#int WAsymm = %.3f ",inthWJetsAsymm), "");
//  leg->AddEntry((TObject*)0, Form("#int WOld = %.3f ",inthWJets), "");
  leg->Draw();
  
  hName="WJetBack"+varName+selName+selName2;
  c2->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());

  TCanvas *c1 = getDefaultCanvas();
  c1->SetName("c1");
  c1->SetTitle("WJets Background estimation");
  hWJetsAsymmMC -> SetStats(kFALSE);

  hWJetsAsymmMC -> Draw("hist");
//  hWJetsAsymmMC->SetMinimum(-2000);
//  hWJetsAsymmMC->SetMaximum(8000);
  hWJets -> SetLineColor(1);
  hWJets->SetLineWidth(1.2);
  hWJets->Draw("same");

  leg->Clear();
  leg->SetTextSize(0.03);
  leg->AddEntry(hWJets,"Old","l");
  leg->AddEntry(hWJetsAsymmMC,"Data AsymmMC","l");
//  leg->AddEntry((TObject*)0, Form("#int WAsymm = %.3f ",inthWJetsAsymm), "");
//  leg->AddEntry((TObject*)0, Form("#int WOld = %.3f ",inthWJets), "");
  leg->Draw();
  
  hName="WJetBackMC"+varName+selName+selName2;
  c1->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());

  return std::make_pair(StatError40, StatError40MC);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTWeightHistograms::getWNormalisation(std::string selName, std::string SubSelName){

  std::cout<<"Calling method: "<<__func__<<std::endl;

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

  // Create a histogram with data minus backgrounds: DYJets, hTT, Other
  TH1F* datamtlo = (TH1F*)hSoup->Clone("datamtlo");
  datamtlo->Add(hDYJets,-1);
  datamtlo->Add(hDYJetsLowM,-1);
  datamtlo->Add(hTT,-1);
  //datamtlo->Add(hQCD,-1);

  float inthWJets=hWJets->Integral(0,hWJets->GetNbinsX()+1);
  float intdata=datamtlo->Integral(0,datamtlo->GetNbinsX()+1);


  double errorW, errorData;
  hWJets->IntegralAndError(0,hWJets->GetNbinsX()+1,errorW,"");
  datamtlo->IntegralAndError(0,datamtlo->GetNbinsX()+1,errorData,"");

  // Calculate weight
  weight=intdata/inthWJets;
//---------test
// TCanvas* c = new TCanvas("QCD_OStoSS","QCD_OStoSS",460,500);

//  hWJets->Scale(weight);
//  hWJets->Draw("hist");
//  datamtlo->Draw("same");

//   hName="TloW"+selName+SubSelName;

//  c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
//  c->Print(TString::Format("fig_C/%s.C",hName.c_str()).Data());
//---------test

  float dweight;
  float inthSoup = hSoup->Integral(0,hSoup->GetNbinsX()+1);
  float inthDYJetsLowM = hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1);
  float inthDYJets = hDYJets->Integral(0,hDYJets->GetNbinsX()+1);
  float inthTT = hTT->Integral(0,hTT->GetNbinsX()+1);
  float inthOther = 0;//hOther->Integral(0,hOther->GetNbinsX()+1);
  dweight=((inthSoup+inthDYJets+inthDYJetsLowM+inthTT+inthOther)/inthWJets/inthWJets+intdata*intdata/(inthWJets*inthWJets*inthWJets));
  dweight=sqrt(dweight);
  cout<<"Selecion name: "<<selName<<std::endl;
  cout<<"DATA: "<<inthSoup<<" DATA - MC(!WJets): "<<intdata<<" error "<<errorData<<" MC WJets "<<inthWJets<<" error "<<errorW
      <<" DYJets: "<<inthDYJets<<" DYJetsLowM: "<<inthDYJetsLowM
      <<" TTbar: "<<inthTT<<" Other: "<<inthOther<<endl;
  cout<<"WJets scale:"<<weight<<" dweight "<<dweight<<endl;
  return std::make_pair(weight, dweight);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
