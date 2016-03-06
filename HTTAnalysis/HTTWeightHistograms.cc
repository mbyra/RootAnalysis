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
   add1DHistogram("h1DMassTemplate",";SVFit mass [GeV/c^{2}]; Events",50,0,200,file_);
   add1DHistogram("h1DEtaTemplate",";#eta; Events",24,-2.4,2.4,file_);
   float bins[31] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5};
   add1DHistogram("h1DIsoTemplate",";Isolation; Events",30,bins,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1* HTTWeightHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

//  PlotAsymm("EtaMuon","","Plus"); 

  WJetEstimation("EtaMuon","","Plus");

  std::pair<TH1*, TH1*> WEstimation = PlotAsymm("EtaMuon","","Plus");
  TH1F * hAsymm = (TH1F*) WEstimation.second;

  return hAsymm;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double HTTWeightHistograms::MakeDiff(TH1F *hTTbar, TH1F* hDYJets, TH1F* hSoup, TH1F* hWJets, TH1F* hQCD, TH1F *hTTbarS, TH1F* hDYJetsS, TH1F* hSoupS, TH1F* hWJetsS, TH1F* hQCDS, std::string varName, std::string selName, std::string SubSelName){

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
    scale = weight*lumi;
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
  std::cout<<"------------------------------------------"<<std::endl;
  std::cout<<"--------dodalam TT-----------------"<<std::endl;
  TH1F *hTTbar12 = get1DHistogram((hName+"TTbar"+"qcdselSS"+"Plus").c_str());
  TH1F *hTTbar22 = get1DHistogram((hName+"TTbar"+"qcdselSS"+"Minus").c_str());
  hTTbar12->Scale(scale);
  hTTbar22->Scale(scale);
  hTTbar1->Add(hTTbar12,1);
  hTTbar2->Add(hTTbar22,1);
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
std::pair<TH1*,TH1*>  HTTWeightHistograms::PlotAsymm(std::string varName, std::string selName, std::string SubSelName){

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

  int rebinFactor = 2;  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hTTbar->Rebin(rebinFactor);
  hDYJets->Rebin(rebinFactor);
  hQCD->Rebin(rebinFactor);

  hSoupS->Rebin(rebinFactor);
  hWJetsS->Rebin(rebinFactor);
  hTTbarS->Rebin(rebinFactor);
  hDYJetsS->Rebin(rebinFactor);
  hQCDS->Rebin(rebinFactor);

  float lumi = getLumi();

// Delta from the data
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);
// I delete it to have W+Jet estimation method without QCD  hSoup->Add(hQCD,-1);

//  TH1F* W  = (TH1F*)hSoup->Clone(("Asymm"+varName+selName).c_str);

  TH1F *W = get1DHistogram("h1DEtaAsymm");
  W  -> Reset();
  W  -> Add(hSoup,1);

  TH1F* WMC  = (TH1F*)hSoup->Clone("WMC");
  WMC  -> Reset();
  WMC  -> Add(hWJets,1);

// here Asymmetry from Data

  hSoupS->Add(hDYJetsS,1);
  hSoupS->Add(hTTbarS,1);
//  hSoupS->Add(hQCDS,1);

  hSoup->Divide(hSoupS);

// here Asymmetry from MC
  hWJets -> Divide(hWJetsS);

// WJet Bac. Estimation from the Data
  W -> Divide(hWJets);

// WJet Bac. Estimation from the MC
  WMC -> Divide(hWJets);

// drawing
  TCanvas *c2 = getDefaultCanvas();
  c2->SetName("c2");
  c2->SetTitle("W asymmetry");

  if(varName.find("SVfit")!=std::string::npos){
    hSoup->GetXaxis()->SetTitle("SVFit mass [GeV/c^{2}]");
    hSoup->SetTitle("asymmetry SVFit mass [GeV/c^{2}]");}

  if(varName.find("Mt")!=std::string::npos){
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
  hSoup->Draw("hist");
  hWJets->SetLineColor(1);
  hWJets->SetLineWidth(1.2);
//  hWJets->SetMinimum(0.0);
//  hWJets->SetMaximum(0.5);
  hWJets->Draw("same");
  
  hName="Asymmetry"+varName+selName;
  c2->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());

// return WJetBack estimations
  return std::make_pair(W, hWJets);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double * HTTWeightHistograms::WJetEstimation(std::string varName, std::string selName, std::string SubSelName){

  std::string hName = "h1D"+varName;

  std::pair<TH1*, TH1*> WEstimation = PlotAsymm(varName,selName,"Plus");
  TH1F * hWJetsAsymm = (TH1F*) WEstimation.first;
//  TH1F * hWJetsAsymmMC = (TH1F*) WEstimation.second;

  std::string selName2="";
  if(selName=="All") {
	selName="";
	selName2="All";
	}

  TH1F *hWJets = get1DHistogram((hName+"WJets"+selName).c_str());
//  TH1F *hWJets = get1DHistogram((hName+"WJets"+selName+SubSelName).c_str());
  TH1F *hWJets2 = get1DHistogram((hName+"WJets"+"qcdselSS").c_str());

  int rebinFactor = 2;  
  hWJets->Rebin(rebinFactor); 
  hWJets2->Rebin(rebinFactor);

// different scalling for OS and SS WJets
  float lumi = getLumi();

  std::string WselType = "wselOS";
  if(selName.find("SS")!=std::string::npos) WselType = "wselSS";
  pair<float,float> dataToMCScale = getWNormalisation(WselType,"");
  float weight = getSampleNormalisation("WJets");
  float scale = weight*lumi*dataToMCScale.first;
  hWJets->Scale(scale);

  if((SubSelName=="Plus" || SubSelName=="Minus") && selName2=="All"){
  std::cout<<"------------------------------------------"<<std::endl;
  std::cout<<"--------dodalam Wreferencje-----------------"<<std::endl;
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

  if(varName.find("Mt")!=std::string::npos){
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
//  leg->AddEntry((TObject*)0, Form("#int WAsymm = %.3f ",inthWJetsAsymm), "");
//  leg->AddEntry((TObject*)0, Form("#int WOld = %.3f ",inthWJets), "");
  leg->Draw();
  
  hName="WJetBack"+varName+selName+selName2;
  c2->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());

//  TCanvas *c1 = getDefaultCanvas();
//  c1->SetName("c1");
//  c1->SetTitle("WJets Background estimation");
//  hWJetsAsymmMC -> SetStats(kFALSE);

//  hWJetsAsymmMC -> Draw("hist");
//  hWJetsAsymmMC->SetMinimum(-2000);
//  hWJetsAsymmMC->SetMaximum(8000);
//  hWJets -> SetLineColor(1);
//  hWJets->SetLineWidth(1.2);
//  hWJets->Draw("same");

//  leg->Clear();
//  leg->SetTextSize(0.03);
//  leg->AddEntry(hWJets,"Old","l");
//  leg->AddEntry(hWJetsAsymmMC,"Data Asymm","l");
//  leg->AddEntry((TObject*)0, Form("#int WAsymm = %.3f ",inthWJetsAsymm), "");
//  leg->AddEntry((TObject*)0, Form("#int WOld = %.3f ",inthWJets), "");
//  leg->Draw();
  
//  hName="WJetBackMC"+varName+selName+selName2;
//  c1->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());

return 0;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTWeightHistograms::getQCDOStoSS(std::string selName, std::string SubSelName){

  std::cout<<"Calling method: "<<__func__<<std::endl;
  if(selName.find("SS")!=std::string::npos) return  std::make_pair(1.0,0.0);

  std::string hName = "h1DIso";

  // SS selection
  TH1F *hWJetsSS = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+SubSelName).c_str());
  TH1F *hDYJetsLowMSS = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+SubSelName).c_str());
  TH1F *hDYJetsSS = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+SubSelName).c_str()); 
  TH1F *hTTSS = get1DHistogram((hName+"TTbar"+"qcdselSS"+SubSelName).c_str());
  TH1F *hSoupSS = get1DHistogram((hName+"Data"+"qcdselSS"+SubSelName).c_str());
  TH1F *hSoupSSb = get1DHistogram((hName+"Data"+"qcdselSS"+SubSelName).c_str());

  // OS selection
  TH1F *hWJetsOS = get1D_WJet_Histogram((hName+"WJets"+"qcdselOS"+SubSelName).c_str());
  TH1F *hDYJetsLowMOS = get1DHistogram((hName+"DYJetsLowM"+"qcdselOS"+SubSelName).c_str());
  TH1F *hDYJetsOS = get1D_DY_Histogram((hName+"DYJets"+"qcdselOS"+SubSelName).c_str());
  TH1F *hTTOS = get1DHistogram((hName+"TTbar"+"qcdselOS"+SubSelName).c_str());
  TH1F *hSoupOS = get1DHistogram((hName+"Data"+"qcdselOS"+SubSelName).c_str());
  TH1F *hSoupOSb = get1DHistogram((hName+"Data"+"qcdselOS"+SubSelName).c_str());

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
  scale = getSampleNormalisation(sampleName); 
  hWJetsOS->Scale(scale*getWNormalisation("wselOS",SubSelName).first);
  hWJetsSS->Scale(scale*getWNormalisation("wselSS",SubSelName).first);

// UWAGA UWAGA UWAGA
// tu mysi byc czytanie histogramow w zaleznosci od ISO!!! 
// na razie robię podmiankę i daję normalne WJety
  TH1F *hWJetsOSAsymm = (TH1F*)hWJetsOS->Clone((hName+"hWJetsOSAsymm"+SubSelName).c_str());
  TH1F *hWJetsSSAsymm = (TH1F*)hWJetsSS->Clone((hName+"hWJetsSSAsymm"+SubSelName).c_str());
//  std::pair<TH1*, TH1*> WEstimationOS = PlotAsymm(varName,"","Plus");
//  TH1F * hWJetsOSAsymm = (TH1F*) WEstimationOS.first;
//  std::pair<TH1*, TH1*> WEstimationSS = PlotAsymm(varName,"qcdselSS","Plus");
//  TH1F * hWJetsSSAsymm = (TH1F*) WEstimationSS.first;
  
  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTOS->Scale(scale);
  hTTSS->Scale(scale);
 
  ///Subtract backgrounds other than QCD using MC
  hSoupSS->Add(hDYJetsLowMSS,-1);
  hSoupSS->Add(hDYJetsSS,-1);
  hSoupSS->Add(hTTSS,-1);

  TH1F *hSoupSSAsymm = (TH1F*)hSoupSS->Clone((hName+"hSoupSSAsymm"+SubSelName).c_str());
  hSoupSSAsymm->Reset();
  hSoupSSAsymm->Add(hSoupSS,1);
  hSoupSS->Add(hWJetsSS,-1);
  hSoupSSAsymm->Add(hWJetsSSAsymm,-1);
  
  hSoupOS->Add(hDYJetsLowMOS,-1);
  hSoupOS->Add(hDYJetsOS,-1);

  hSoupOS->Add(hTTOS,-1);

  TH1F *hSoupOSAsymm = (TH1F*)hSoupOS->Clone((hName+"hSoupOSAsymm"+SubSelName).c_str());
  hSoupOSAsymm->Reset();
  hSoupOSAsymm->Add(hSoupOS,1);
  hSoupOS->Add(hWJetsOS,-1);
  hSoupOSAsymm->Add(hWJetsOSAsymm,-1);

  hSoupOS->Divide(hSoupSS);
  hSoupOSAsymm->Divide(hSoupSSAsymm);

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
  hSoupOS->Draw();
  hSoupOS->Fit("line","","",0.2,0.3);
  c->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());
  c->Print(TString::Format("fig_C/%s.C",hName.c_str()).Data());

  float param, dparam;
  param=line->GetParameter(0);
  dparam=line->GetParError(0);

  //funtion fitting for Asymmetry
  TF1 *lineA=new TF1("lineA","[0]",0,2);
  lineA->SetParameter(0,1);
//  TCanvas* c = new TCanvas("QCD_OStoSS","QCD_OStoSS",460,500);
//    TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->Clear();
   leg->Draw();

  hSoupOSAsymm->SetLineWidth(3);
  hSoupOSAsymm->GetYaxis()->SetTitleOffset(1.4);
  hSoupOSAsymm->GetYaxis()->SetTitle("OS/SS");
  hSoupOSAsymm->GetXaxis()->SetTitle("muon relative isolation");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(11);
  hSoupOSAsymm->Draw();
  hSoupOSAsymm->Fit("lineA","","",0.2,0.4);

  hName=hName+"_Asymm";
  c->Print(TString::Format("fig_pngA/%s.png",hName.c_str()).Data());
  c->Print(TString::Format("fig_C/%s.C",hName.c_str()).Data());

  float paramA, dparamA;
  paramA=line->GetParameter(0);
  dparamA=line->GetParError(0);

//

  std::cout<<"QCD OS/SS ratio: "<<param<<" +- "<<dparam<<std::endl;
  std::cout<<"QCD OS/SS ratio Asymmetry: "<<paramA<<" +- "<<dparamA<<std::endl;

  return std::make_pair(param, paramA);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<TH1*,TH1*> HTTWeightHistograms::getQCDbackground(std::string varName, std::string selName, std::string SubSelName){

  std::cout<<"Calling method: "<<__func__<<std::endl;
  
  ///Not very clear and elegant. AK
  ///Need this to avoid resursive control region labels like
  ///qcdselSSqcdselOS
  int ifscale=0; //when SS there should be no scaling
  if(selName.find("SS")==std::string::npos) ifscale=1;
  if(selName.find("qcdsel")!=std::string::npos) selName = "";

  std::string hName = "h1D" + varName;
  // SS selection

  std::pair<TH1*, TH1*> WEstimationSS = PlotAsymm(varName,"qcdselSS","Plus");
  TH1F * hWJetsAsymm = (TH1F*) WEstimationSS.first;

  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"qcdselSS"+selName+SubSelName).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+"qcdselSS"+selName+SubSelName).c_str());


  ///Protection against null pointers
  ///Null pointers happen when sample was not read, or there were no
  ///events passing particular selection.
//  if(!hSoup) return 0;
  if(!hWJets){
    hWJets = (TH1F*)hSoup->Clone((hName+"WJets"+"qcdselSS"+selName+SubSelName).c_str()); hWJets->Reset();
  }

  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hSoup->Clone((hName+"hDYJetsLowM"+"qcdselSS"+selName+SubSelName).c_str()); hDYJetsLowM->Reset();
  }
  if(!hDYJets){
    hDYJets = (TH1F*)hSoup->Clone((hName+"hDYJets"+"qcdselSS"+selName+SubSelName).c_str()); hDYJets->Reset();
  }
  if(!hTTbar){
    hTTbar = (TH1F*)hSoup->Clone((hName+"hTTbar"+"qcdselSS"+selName+SubSelName).c_str()); hTTbar->Reset();
  }
  //////////////////////////////////////////////////////////////////////
  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections

  std::string sampleName = "DYJetsLowM";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJetsLowM->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "WJets";
  float dataToMCScale = getWNormalisation("wselSS",SubSelName).first;
  scale = getSampleNormalisation(sampleName)*lumi*dataToMCScale;
  hWJets->Scale(scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar->Scale(scale);


  TH1F *hSoupAsymm = (TH1F*)hSoup->Clone((hName+"hSoupAsymm"+SubSelName).c_str());
  hSoupAsymm->Reset();

  hSoup->SetName(("h1DQCDEstimate"+varName+SubSelName).c_str());
  hSoup->Add(hWJets,-1);
  hSoup->Add(hDYJetsLowM,-1);
  hSoup->Add(hDYJets,-1);

  hSoup->Add(hTTbar,-1);

  hSoupAsymm->Add(hSoup,1);

//  TH1F *hSoupAsymm = (TH1F*)hSoup->Clone((hName+"hSoupAsymm"+SubSelName).c_str());
  hSoup->Add(hWJets,-1);
  hSoupAsymm->Add(hWJetsAsymm,-1); // tu bylo hWJetsAsymm

  if(ifscale==1){
  hSoup->Scale(getQCDOStoSS(selName,SubSelName).first);
  hSoupAsymm->Scale(getQCDOStoSS(selName,SubSelName).second);
  }

  ///Clean up the QCD shape, and remove fluctuations around 0 counts.
  for(unsigned int iBinX=0;iBinX<=hSoup->GetNbinsX();++iBinX){
    if(hSoup->GetBinContent(iBinX)<3.0) hSoup->SetBinContent(iBinX,0);
  }

  for(unsigned int iBinX=0;iBinX<=hSoupAsymm->GetNbinsX();++iBinX){
    if(hSoupAsymm->GetBinContent(iBinX)<3.0) hSoupAsymm->SetBinContent(iBinX,0);
  }

  TCanvas* c = new TCanvas("QCD_OStoSS","QCD_OStoSS",460,500);

   TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->Clear();
   leg->Draw();
   hSoup->SetStats(kFALSE);
   hSoup->SetTitle(("QCDBack"+varName+"; selection: "+selName+SubSelName).c_str());
   hSoup->GetXaxis()->SetTitle(varName.c_str());

   hSoup->Draw("hist");
   hSoupAsymm->SetLineColor(2);
   hSoupAsymm->Draw("same");
  std::string name="QCDBack"+varName+selName+SubSelName;
  c->Print(TString::Format("fig_pngA/%s.png",name.c_str()).Data());
  c->Print(TString::Format("fig_C/%s.C",name.c_str()).Data());

  return std::make_pair(hSoup, hSoupAsymm);
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
  cout<<"DATA: "<<inthSoup<<" DATA - MC(!WJets): "<<intdata<<" MC WJets "<<inthWJets
      <<" DYJets: "<<inthDYJets<<" DYJetsLowM: "<<inthDYJetsLowM
      <<" TTbar: "<<inthTT<<" Other: "<<inthOther<<endl;
  cout<<"WJets scale:"<<weight<<" dweight "<<dweight<<endl;
  return std::make_pair(weight, dweight);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
