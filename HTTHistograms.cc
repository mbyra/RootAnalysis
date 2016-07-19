#include <iostream>
#include <cmath>

#include "commonUtils.h"
#include "HTTHistograms.h"
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
float HTTHistograms::getLumi(){

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
float HTTHistograms::getSampleNormalisation(std::string sampleName){

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
HTTHistograms::HTTHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){

  selectionFlavours_ = flavours;

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::~HTTHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1D_DY_Histogram(const std::string& name){
  
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
TH1F *HTTHistograms::get1D_WJet_Histogram(const std::string& name){

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
bool HTTHistograms::fill1DHistogram(const std::string& name, float val, float weight){

  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill1DHistogram(name,val,weight)){
    if(name.find("h1DNPV")!=std::string::npos) hTemplateName = "h1DNPVTemplate";
    if(name.find("h1DMass")!=std::string::npos) hTemplateName = "h1DMassTemplate";
    if(name.find("h1DStats")!=std::string::npos) hTemplateName = "h1DStatsTemplate";
    if(name.find("h1DPt")!=std::string::npos) hTemplateName = "h1DPtTemplate";
    if(name.find("h1DEta")!=std::string::npos) hTemplateName = "h1DEtaTemplate";
    if(name.find("h1DIso")!=std::string::npos) hTemplateName = "h1DIsoTemplate";
    if(name.find("h1DPhi")!=std::string::npos) hTemplateName = "h1DPhiTemplate";
    if(name.find("h1DCosPhi")!=std::string::npos) hTemplateName = "h1DCosPhiTemplate";
    if(name.find("h1DCSVBtag")!=std::string::npos) hTemplateName = "h1DCSVBtagTemplate";
    if(name.find("h1DID")!=std::string::npos) hTemplateName = "h1DIDTemplate";
    //std::cout<<"fill1DHistogram Adding histogram: "<<name<<" "<<file_<<" "<<file_->fullPath()<<std::endl;
    
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
void HTTHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   std::cout<<"defineHistograms Adding histogram: "<<file_<<" "<<file_->fullPath()<<std::endl;

   add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
   add1DHistogram("h1DNPVTemplate",";Number of PV; Events",61,-0.5,60.5,file_);
   add1DHistogram("h1DMassTemplate",";SVFit mass [GeV/c^{2}]; Events",50,0,200,file_); 
   add1DHistogram("h1DPtTemplate",";p_{T}; Events",20,0,100,file_);
   add1DHistogram("h1DEtaTemplate",";#eta; Events",24,-2.4,2.4,file_);
   add1DHistogram("h1DPhiTemplate",";#phi; Events",8,-M_PI,M_PI,file_);
   add1DHistogram("h1DCosPhiTemplate",";cos(#phi); Events",10,-1.0,1.0,file_);
   add1DHistogram("h1DCSVBtagTemplate",";CSV btag; Events",20,0,1,file_);
   float bins[31] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5};
   add1DHistogram("h1DIsoTemplate",";Isolation; Events",30,bins,file_);
   add1DHistogram("h1DIDTemplate",";ID; Events",30,-0.5,15.5,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

  plotPhiDecayPlanes("Phi_nVectorsData");
  plotPhiDecayPlanes("Phi_nVectorsH");
  plotPhiDecayPlanes("Phi_nVectorsA");
  plotPhiDecayPlanes("Phi_nVectorsDYJetsMuTau");
  plotPhiDecayPlanes("CosPhi_CosPositiveH");
  plotPhiDecayPlanes("CosPhi_CosNegativeH");
    
  ///Control regions plots
  plotStack("Iso","qcdselOS");
  plotStack("Iso","qcdselSS");
  plotStack("StatsDecayMode","");
  
  plotStack("MassVis","qcdselSS");
  plotStack("StatsNJets30","qcdselSS");
  plotStack("CSVBtagLeadingJet","qcdselSS");
  
  plotStack("MassTrans","");
  plotStack("MassTrans","qcdselSS");
 
  plotStack("EtaMuon","");
  plotStack("EtaMuon","qcdselSS");

  WJetAsymm("MassTrans","","MT");
  WJetAsymm("MassTrans","","Eta");
  WJetAsymm("MassTrans","All","MT");
  WJetAsymm("MassTrans","All","Eta");
  WJetAsymm("MassTrans","qcdselSS","MT");
  WJetAsymm("MassTrans","qcdselSS","Eta");

  WJetAsymm("EtaMuon","","MT");
  WJetAsymm("EtaMuon","","Eta");
  WJetAsymm("EtaMuon","All","MT");
  WJetAsymm("EtaMuon","All","Eta");
  WJetAsymm("EtaMuon","qcdselSS","MT");
  WJetAsymm("EtaMuon","qcdselSS","Eta");

  WJetAsymm("MassVis","","MT");
  WJetAsymm("MassVis","","Eta");
  WJetAsymm("MassVis","All","MT");
  WJetAsymm("MassVis","All","Eta");
  WJetAsymm("MassVis","qcdselSS","MT");
  WJetAsymm("MassVis","qcdselSS","Eta");

  plotStack("MassTrans","wselOS");  
  plotStack("MassTrans","wselSS");

  plotStack("MassTrans","ttselOS");  
  plotStack("MassTrans","ttselSS");
  plotStack("IsoMuon","ttselOS");

  plotStack("PtMET","mumuselOS");
  plotStack("PtMET","mumuselSS");

  ///Baseline selection plots
  plotStack("MassSV","");
  plotStack("MassVis","");  
    plotStack("MassVis","qcdselSS");
  plotStack("MassVis","AsymmEta");  
  plotStack("MassVis","AsymmMT"); 

  plotStack("PtMuon","");
  plotStack("EtaMuon","");
  plotStack("IsoMuon","");
  
  plotStack("PtTau","");  
  plotStack("EtaTau","");
  plotStack("IDTau","");
  plotStack("StatsDecayMode","");
  plotStack("PtTauLeadingTk","");
 
  plotStack("PhiMuon","");
  plotStack("PhiTau","");

  plotStack("PtMET","");  

  plotStack("StatsNJets30","");
  
  plotStack("PtLeadingJet","");
  plotStack("EtaLeadingJet","");
  plotStack("CSVBtagLeadingJet","");
  
  plotStack("PtLeadingBJet","");
  plotStack("EtaLeadingBJet","");

  plotStack("NPV","");


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotPhiDecayPlanes(const std::string & name){

  TCanvas aCanvas(TString::Format("PhiDecayPlanes_%s",name.c_str()),
		  TString::Format("PhiDecayPlanes_%s",name.c_str()),
		  460,500);
  
  TLegend l(0.15,0.15,0.35,0.37,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TString hName = "h1D"+name+"RefitPV";
  TH1F* h1DRefitPV = this->get1DHistogram(hName.Data());

  hName = "h1D"+name+"AODPV";
  TH1F* h1DAODPV = this->get1DHistogram(hName.Data());
  
  hName = "h1D"+name+"GenPV";
  TH1F* h1DGenPV = this->get1DHistogram(hName.Data());

  hName = "h1D"+name+"GenNoOfflineSel";
  TH1F* h1DGen = this->get1DHistogram(hName.Data());

  if(h1DGen){
    h1DGen->SetLineWidth(3);
    h1DGen->Scale(1.0/h1DGen->Integral(0,h1DGen->GetNbinsX()+1));
    h1DGen->SetLineColor(1);
  }

  if(h1DAODPV){
    h1DAODPV->SetLineWidth(3);
    h1DAODPV->Scale(1.0/h1DAODPV->Integral(0,h1DAODPV->GetNbinsX()+1));
    h1DAODPV->SetLineColor(2);
  }

  if(h1DGenPV){
    h1DGenPV->SetLineWidth(3);
    h1DGenPV->Scale(1.0/h1DGenPV->Integral(0,h1DGenPV->GetNbinsX()+1));
    h1DGenPV->SetLineColor(3);
  }
  
  if(h1DRefitPV){
    h1DRefitPV->SetLineWidth(3);
    h1DRefitPV->SetLineColor(4);    
    h1DRefitPV->Scale(1.0/h1DRefitPV->Integral(0,h1DRefitPV->GetNbinsX()+1));
    h1DRefitPV->SetXTitle("#phi^{*}");
    h1DRefitPV->SetYTitle("Events");
    h1DRefitPV->SetTitle(name.c_str());
    h1DRefitPV->GetYaxis()->SetTitleOffset(1.4);
    h1DRefitPV->SetStats(kFALSE);
    //h1DRefitPV->GetXaxis()->SetRangeUser(0,M_PI);

    h1DRefitPV->SetMaximum(1.02);
    //h1DRefitPV->SetMaximum(0.4);
    //h1DRefitPV->SetMinimum(0.10);
    h1DRefitPV->Draw("HISTO");    
    l.AddEntry(h1DRefitPV,"reco PCA with refit. PV");
    if(h1DGenPV){
      h1DGenPV->Draw("HISTO same");
      l.AddEntry(h1DGenPV,"reco PCA with gen. PV");
    }
    
    if(h1DAODPV){
      h1DAODPV->Draw("HISTO same");
      l.AddEntry(h1DAODPV,"reco PCA with AOD PV");
    }
    if(h1DGen){
      h1DGen->Draw("HISTO same");
      l.AddEntry(h1DGen,"#splitline{PCA with gen. particles}{no offline selection}");
    }
    l.Draw();
    aCanvas.Print(TString::Format("fig_png/%s.png",name.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
THStack*  HTTHistograms::plotStack(std::string varName, std::string selName){

  std::cout<<"--- Drawing THStack for variable: "<<varName
	   <<" selection: "<<selName<<std::endl;

  std::string hName = "h1D"+varName;

  TH1F *hHiggs = get1DHistogram((hName+"H"+selName).c_str());

  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+selName).c_str());

  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+selName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+selName).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+selName).c_str());
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+selName).c_str());

  pair<float,float> qcdOStoSS = getQCDOStoSS(selName,"");

  TH1F *hQCD = (TH1F*)getQCDbackground(varName,selName);

  ///Protection against null pointers
  ///Null pointers happen when sample was not read, or there were no
  ///events passing particular selection.

  if(!hSoup) return 0;
  
  if(!hQCD){
    hQCD = (TH1F*)hSoup->Clone((hName+"QCD"+selName).c_str()); hQCD->Reset();
  }
  if(!hWJets){
    hWJets = (TH1F*)hSoup->Clone((hName+"WJets"+selName).c_str()); hWJets->Reset();
  }
  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hSoup->Clone((hName+"hDYJetsLowM"+selName).c_str()); hDYJetsLowM->Reset();
  }  
  if(!hTTbar){
    hTTbar = (TH1F*)hSoup->Clone((hName+"hTTbar"+selName).c_str()); hTTbar->Reset();
  }
  if(!hHiggs){
    hHiggs = (TH1F*)hSoup->Clone((hName+"hH"+selName).c_str()); hHiggs->Reset();
  } 
  /////////////////////////////////////////////////////////////////

  float lumi = getLumi(); 
  std::string sampleName = "WJets";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight * lumi;

	std::string WselType = "wselOS";
	if(selName.find("SS")!=std::string::npos) WselType = "wselSS";
	pair<float,float> dataToMCScale = getWNormalisation(WselType,"");

 	hWJets->Scale(scale);
  	hWJets->Scale(dataToMCScale.first);
 

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
  hTTbar->Scale(scale);

  sampleName = "H";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hHiggs->Scale(scale);
  //////////////////////////////////////////////////////
  hSoup->SetLineColor(1);
  hSoup->SetFillColor(1);
  hSoup->SetMarkerStyle(20);

  hWJets->SetFillColor(kRed+2);
  hTTbar->SetFillColor(kBlue+2);
  hDYJets->SetFillColor(kOrange-1);
  hDYJetsLowM->SetFillColor(kOrange-7);

  hQCD->SetFillColor(kMagenta-10);
  hHiggs->SetFillColor(kCyan+4);

  hSoup->SetLineWidth(1);
  int rebinFactor = 1;  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hTTbar->Rebin(rebinFactor);
  hDYJets->Rebin(rebinFactor);
  hDYJetsLowM->Rebin(rebinFactor);
  hHiggs->Rebin(rebinFactor);

  THStack *hs = new THStack("hs","Stacked histograms");      
  /////////
  hs->Add(hHiggs,"hist");    
  hs->Add(hQCD,"hist");
  hs->Add(hTTbar,"hist");
  hs->Add(hWJets,"hist");
  hs->Add(hDYJetsLowM,"hist");
  hs->Add(hDYJets,"hist");
  ////////
  TH1F *hMCSum = (TH1F*)hWJets->Clone("hMCSum");
  hMCSum->Reset();
  hMCSum->Add(hDYJetsLowM);
  hMCSum->Add(hDYJets);
  hMCSum->Add(hWJets);
  hMCSum->Add(hTTbar);
  hMCSum->Add(hQCD);
  hMCSum->Add(hHiggs);

  if(!selName.size()) selName = "baseline";
  cout<<"Event count summary for selecion name: "<<selName<<std::endl;
  std::cout<<"Data: "<<hSoup->Integral(0,hSoup->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC: "<<hMCSum->Integral(0,hMCSum->GetNbinsX()+1)<<std::endl;  
  std::cout<<"MC W->l: "<<hWJets->Integral(0,hWJets->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC TTbar: "<<hTTbar->Integral(0,hTTbar->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->: "<<hDYJets->Integral(0,hDYJets->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->ll(m<50): "<<hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC H->tau tau: "<<hHiggs->Integral(0,hHiggs->GetNbinsX()+1)<<std::endl;  
  std::cout<<"QCD: "<<hQCD->Integral(0,hQCD->GetNbinsX()+1)<<std::endl; 
  std::cout<<"Correction factors:"<<std::endl;
  std::cout<<"QCD SS to OS: "<<qcdOStoSS.first<<" +- "<<qcdOStoSS.second<<std::endl;
  std::cout<<"W DATA to MC: "<<dataToMCScale.first<<" +- "<<dataToMCScale.second<<std::endl;
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
//  hs->SetTitle(("Variable: "+varName+" selection: "+selName).c_str());
  hs -> SetTitle("");
  hs->SetMaximum(4400);
  hs->Draw("hist");
  hs->GetXaxis()->SetTitle(varName.c_str());

  if(varName.find("SVfit")!=std::string::npos){
    hs->GetXaxis()->SetTitle("SVFit mass [GeV]");
    }

  if(varName.find("MassTrans")!=std::string::npos){
    hs->GetXaxis()->SetTitle("M_{T} [GeV]");
    }

  if(varName.find("NPV")!=std::string::npos){
    hs->GetXaxis()->SetTitle("NPV");
    }

  if(varName.find("Eta")!=std::string::npos){
    hs->GetXaxis()->SetTitle("#eta [units]");
    }

  if(varName.find("MassVis")!=std::string::npos){
    hs->GetXaxis()->SetTitle("M_{VIS} [GeV]");
    }

  hs->GetYaxis()->SetTitleOffset(1.4);
  hMCSum->SetFillColor(5);
  /////////
  float highEnd = 150;
  float lowEnd = -150;

  int binHigh = hs->GetXaxis()->FindBin(highEnd);  
  int binLow = hs->GetXaxis()->FindBin(lowEnd);

  if(hs->GetXaxis()->GetXmax()<highEnd || hs->GetXaxis()->GetXmax()>300) binHigh = hs->GetXaxis()->GetNbins();
  if(hs->GetXaxis()->GetXmin()>lowEnd) lowEnd = 1;

  hs->GetXaxis()->SetRange(binLow,binHigh);
  highEnd =  hs->GetXaxis()->GetBinUpEdge(binHigh);

  char yTitle[200];
  sprintf(yTitle,"Number of events");
  hs->GetYaxis()->SetTitle(yTitle);
  
  float max = hs->GetMaximum();
  if(hSoup->GetMaximum()>max) max = hSoup->GetMaximum();

  hs->GetHistogram()->SetTitleOffset(1.0);
  hs->SetMaximum(1.1*max);
  hs->SetMinimum(0.1);

  hSoup->Draw("same");
  TH1F *hEmpty = new TH1F("hEmpty","",1,0,1);
  hEmpty->SetLineColor(10);
  hEmpty->SetFillColor(10);

  TLegend *leg = new TLegend(0.79,0.32,0.99,0.82,NULL,"brNDC");
  setupLegend(leg);
  leg->AddEntry(hSoup,"Data","lep");
  leg->AddEntry(hDYJets,"Z#rightarrow ll","f");
  leg->AddEntry(hDYJetsLowM,"Z#rightarrow ll(m<50)","f");
  leg->AddEntry(hWJets,"W#rightarrow l #nu","f");
  leg->AddEntry(hTTbar,"TTbar","f");
  leg->AddEntry(hQCD,"QCD","f");
  leg->AddEntry(hHiggs,"H#rightarrow #tau #tau","f");
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
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/"))+ ".png";    
  else plotName = "fig_png/hTree_"+hName+Form("_%s",selName.c_str())+".png";
  c1->Print(plotName.c_str());

  if(hName.find_last_of("/")<string::npos) plotName = "fig_C/" + hName.substr(hName.find_last_of("/"))+ ".C";    
  else plotName = "fig_C/hTree_"+hName+Form("_%s",selName.c_str())+".C";
  c1->Print(plotName.c_str()); 

  pad1->SetLogy(1);
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/"))+ "_LogY.png";    
  else plotName = "fig_png/hTree_"+hName+Form("_%s",selName.c_str())+"_LogY.png";
  c1->Print(plotName.c_str()); 

  std::cout<<"-------------------------------------------------------------"<<std::endl;

  return hs;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotSingleHistogram(std::string hName){
  
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",			   
			   460,500);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  
  TH1F* h1D = this->get1DHistogram(hName.c_str());
  
  if(h1D){
    h1D->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw();
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getQCDOStoSS(std::string selName,std::string subSelName){

  std::cout<<"Calling method: "<<__func__<<std::endl;
  if(selName.find("SS")!=std::string::npos) return  std::make_pair(1.0,0.0);

  std::string hName = "h1DIso";

    selName = "";

  // SS selection
  TH1F *hDYJetsLowMSS = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+selName+subSelName).c_str());
  TH1F *hDYJetsSS = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+selName+subSelName).c_str());  
  TH1F *hTTSS = get1DHistogram((hName+"TTbar"+"qcdselSS"+selName+subSelName).c_str());
  TH1F *hSoupSS = get1DHistogram((hName+"Data"+"qcdselSS"+selName+subSelName).c_str());

  if(!hDYJetsLowMSS){
    hDYJetsLowMSS = (TH1F*)hSoupSS->Clone((hName+"DYJetsLowM"+"qcdselSS"+selName+subSelName).c_str()); hDYJetsLowMSS->Reset();
  } 
  if(!hTTSS){
    hTTSS = (TH1F*)hSoupSS->Clone((hName+"hTTbar"+"qcdselSS"+selName+subSelName).c_str()); hTTSS->Reset();
  } 

  // OS selection
  TH1F *hDYJetsLowMOS = get1DHistogram((hName+"DYJetsLowM"+"qcdselOS"+selName+subSelName).c_str());
  TH1F *hDYJetsOS = get1D_DY_Histogram((hName+"DYJets"+"qcdselOS"+selName+subSelName).c_str());
  TH1F *hTTOS = get1DHistogram((hName+"TTbar"+"qcdselOS"+selName+subSelName).c_str());
  TH1F *hSoupOS = get1DHistogram((hName+"Data"+"qcdselOS"+selName+subSelName).c_str());

  if(!hDYJetsLowMOS){
    hDYJetsLowMOS = (TH1F*)hSoupOS->Clone((hName+"DYJetsLowM"+"qcdselOS"+selName+subSelName).c_str()); hDYJetsLowMOS->Reset();
  }

  if(!hTTOS){
    hTTOS = (TH1F*)hSoupOS->Clone((hName+"hTTbar"+"qcdselOS"+selName+subSelName).c_str()); hTTOS->Reset();
  } 

  float lumi = getLumi();

  // W+Jet
  TH1F *hWJetsOS = get1D_WJet_Histogram((hName+"WJets"+"qcdselOS"+selName+subSelName).c_str());
  TH1F *hWJetsSS = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+selName+subSelName).c_str());

  std::string sampleName = "WJets";
  float weight = getSampleNormalisation(sampleName);
  float scale= weight*lumi;
  	hWJetsOS->Scale(getWNormalisation("wselOS",subSelName).first*scale);
  	hWJetsSS->Scale(getWNormalisation("wselSS",subSelName).first*scale);

  ///Normalise MC histograms according to cross sections
  sampleName = "DYJetsLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowMOS->Scale(scale);
  hDYJetsLowMSS->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsOS->Scale(scale);
  hDYJetsSS->Scale(scale);
  
  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTOS->Scale(scale);
  hTTSS->Scale(scale);
 
  ///Subtract backgrounds other than QCD using MC
  hSoupSS->Add(hWJetsSS,-1);
  hSoupSS->Add(hDYJetsLowMSS,-1);
  hSoupSS->Add(hDYJetsSS,-1);
  hSoupSS->Add(hTTSS,-1);
  
  hSoupOS->Add(hWJetsOS,-1);
  hSoupOS->Add(hDYJetsLowMOS,-1);
  hSoupOS->Add(hDYJetsOS,-1);
  hSoupOS->Add(hTTOS,-1);
  
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

  hName = hName + "_" + selName +"_"+ subSelName;
  c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());

  float param, dparam;
  param=line->GetParameter(0);
  dparam=line->GetParError(0);

  std::cout<<"QCD OS/SS ratio: "<<param<<" +- "<<dparam<<std::endl;

  return std::make_pair(param, dparam);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::getQCDbackground(std::string varName, std::string selName){

  std::cout<<"Calling method: "<<__func__<<std::endl;

  float qcdScale = getQCDOStoSS(selName,"").first;
  
  ///Not very clear and elegant. AK
  ///Need this to avoid resursive control region labels like
  ///qcdselSSqcdselOS
  std::string imie = "";
  if(selName.find("qcdsel")!=std::string::npos) {selName = ""; imie="qcdselSS";}

  std::string hName = "h1D" + varName;
  // SS selection
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+selName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+selName).c_str());
  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+selName).c_str());
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"qcdselSS"+selName).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+"qcdselSS"+selName).c_str());

  ///Protection against null pointers
  ///Null pointers happen when sample was not read, or there were no
  ///events passing particular selection.
  if(!hSoup) return 0;
  if(!hWJets){
    hWJets = (TH1F*)hSoup->Clone((hName+"WJets"+"qcdselSS"+selName).c_str()); hWJets->Reset();
  }
  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hSoup->Clone((hName+"hDYJetsLowM"+"qcdselSS"+selName).c_str()); hDYJetsLowM->Reset();
  }
  if(!hDYJets){
    hDYJets = (TH1F*)hSoup->Clone((hName+"hDYJets"+"qcdselSS"+selName).c_str()); hDYJets->Reset();
  }
  if(!hTTbar){
    hTTbar = (TH1F*)hSoup->Clone((hName+"hTTbar"+"qcdselSS"+selName).c_str()); hTTbar->Reset();
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
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hWJets->Scale(getWNormalisation("wselSS","").first * scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar->Scale(scale);

  hSoup->SetName(("h1DQCDEstimate"+varName).c_str());

  hSoup->Add(hWJets,-1);
  hSoup->Add(hDYJetsLowM,-1);
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);

  ///Clean up the QCD shape, and remove fluctuations around 0 counts.
//  for(unsigned int iBinX=0;iBinX<=hSoup->GetNbinsX();++iBinX){
//    if(hSoup->GetBinContent(iBinX)<3.0) hSoup->SetBinContent(iBinX,0);
// }
  
  hSoup->Scale(qcdScale);

  TCanvas* c = new TCanvas("WEstim","QCD_OStoSS",460,500);
  std::string nazwa = "QCDBack_"+varName+"_" + selName+"_"+imie;
  hSoup->Draw();
  c->Print(TString::Format("fig_png/%s.png",nazwa.c_str()).Data());

  return hSoup;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getWNormalisation(std::string selName,std::string subSelName){

  std::cout<<"Calling method: "<<__func__<<std::endl;

  std::string hName = "h1DMassTrans";
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+selName+subSelName).c_str());
  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+selName+subSelName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+selName+subSelName).c_str());
  TH1F *hTT = get1DHistogram((hName+"TTbar"+selName+subSelName).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+selName+subSelName).c_str());
  float lumi = getLumi();

  if(!hDYJets){
    hDYJets = (TH1F*)hWJets->Clone((hName+"hDYJets"+selName+subSelName).c_str()); hDYJets->Reset();
  }
  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hWJets->Clone((hName+"hDYJetsLowM"+selName+subSelName).c_str()); hDYJetsLowM->Reset();
  }
  if(!hTT){
    hTT = (TH1F*)hWJets->Clone((hName+"hTTbar"+selName+subSelName).c_str()); hTT->Reset();
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

  float inthWJets=hWJets->Integral(0,hWJets->GetNbinsX()+1);
  float intdata=datamtlo->Integral(0,datamtlo->GetNbinsX()+1);

  // Calculate weight
  weight=intdata/inthWJets;
  float dweight;
  float inthSoup = hSoup->Integral(0,hSoup->GetNbinsX()+1);
  float inthDYJetsLowM = hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1);
  float inthDYJets = hDYJets->Integral(0,hDYJets->GetNbinsX()+1);
  float inthTT = hTT->Integral(0,hTT->GetNbinsX()+1);
  float inthOther = 0;  //hOther->Integral(0,hOther->GetNbinsX()+1);
  dweight=((inthSoup+inthDYJets+inthTT+inthOther)/inthWJets/inthWJets+intdata*intdata/(inthWJets*inthWJets*inthWJets));
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
TH1* HTTHistograms::WJetAsymm(std::string varName, std::string selName, std::string subSelName){

// This method estimate W+Jet background using W charge asymmetry.
// SubSelName: Eta, MT.
// OS / SS / All samples

  std::cout<<"--- Drawing WJetAsymm for variable: "<<varName
	   <<" selection: "<<selName<<std::endl;

  std::string selName2="";

  if(selName.find("All")!=std::string::npos) {
	selName="";
	selName2="All";
	}

  std::string hName = "h1D"+varName;

// added histograms
  TH1F *hTTbar1 = get1DHistogram((hName+"TTbar"+selName+"Asymm"+subSelName+"Plus").c_str());
  TH1F *hSoup1 = get1DHistogram((hName+"Data"+selName+"Asymm"+subSelName+"Plus").c_str());
  TH1F *hDYJets1 = get1D_DY_Histogram((hName+"DYJets"+selName+"Asymm"+subSelName+"Plus").c_str());
  TH1F *hDYJetsLowM1 = get1DHistogram((hName+"DYJetsLowM"+selName+"Asymm"+subSelName+"Plus").c_str());

  if(!hDYJets1){
    hDYJets1 = (TH1F*)hSoup1->Clone((hName+"hDYJets1"+selName+"Asymm"+subSelName+"Plus").c_str()); hDYJets1->Reset();
  }
  if(!hDYJetsLowM1){
    hDYJetsLowM1 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM1"+selName+"Asymm"+subSelName+"Plus").c_str()); hDYJetsLowM1->Reset();
  }

  if(!hTTbar1){
    hTTbar1 = (TH1F*)hSoup1->Clone((hName+"hTTbar1"+"qcdselSS").c_str()); hTTbar1->Reset();
  }

  TH1F *hTTbar2 = get1DHistogram((hName+"TTbar"+selName+"Asymm"+subSelName+"Minus").c_str()); 
  TH1F *hSoup2 = get1DHistogram((hName+"Data"+selName+"Asymm"+subSelName+"Minus").c_str());
  TH1F *hDYJets2 = get1D_DY_Histogram((hName+"DYJets"+selName+"Asymm"+subSelName+"Minus").c_str());
  TH1F *hDYJetsLowM2 = get1DHistogram((hName+"DYJetsLowM"+selName+"Asymm"+subSelName+"Minus").c_str());

  if(!hDYJets2){
    hDYJets2 = (TH1F*)hSoup1->Clone((hName+"hDYJets2"+selName+"Asymm"+subSelName+"Minus").c_str()); hDYJets2->Reset();
  }
  if(!hDYJetsLowM2){
    hDYJetsLowM2 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM2"+selName+"Asymm"+subSelName+"Minus").c_str()); hDYJetsLowM2->Reset();
  }

  if(!hTTbar2){
    hTTbar2 = (TH1F*)hSoup2->Clone((hName+"hTTbar2"+"qcdselSS").c_str()); hTTbar2->Reset();
  }

// sum OS and SS histograms for "All" asymmetry

  if(selName2=="All"){
  TH1F *hSoup12 = get1DHistogram((hName+"Data"+"qcdselSS"+"Asymm"+subSelName+"Plus").c_str());
  TH1F *hSoup22 = get1DHistogram((hName+"Data"+"qcdselSS"+"Asymm"+subSelName+"Minus").c_str());
  hSoup1->Add(hSoup12,1);
  hSoup2->Add(hSoup22,1);
  }

  float lumi = getLumi();

  std::string sampleName = "DYJetsLowM";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJetsLowM1->Scale(scale);
  hDYJetsLowM2->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets1->Scale(scale);
  hDYJets2->Scale(scale);

 if(selName2=="All"){
  TH1F *hDYJets12 = get1D_DY_Histogram((hName+"DYJets"+"qcdselSSAsymm"+subSelName+"Plus").c_str());
  TH1F *hDYJetsLowM12 = get1DHistogram((hName+"DYJetsLowM"+"qcdselSSAsymm"+subSelName+"Plus").c_str());

  if(!hDYJets12){
    hDYJets12 = (TH1F*)hSoup1->Clone((hName+"hDYJets12"+selName).c_str()); hDYJets12->Reset();
  }
  if(!hDYJetsLowM12){
    hDYJetsLowM12 = (TH1F*)hSoup1->Clone((hName+"hDYJetsLowM12"+selName).c_str()); hDYJetsLowM12->Reset();
  } 

  TH1F *hDYJets22 = get1D_DY_Histogram((hName+"DYJets"+"qcdselSSAsymm"+subSelName+"Minus").c_str());
  TH1F *hDYJetsLowM22 = get1DHistogram((hName+"DYJetsLowM"+"qcdselSSAsymm"+subSelName+"Minus").c_str());

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

   if(selName2=="All"){
	  TH1F *hTTbar12 = get1DHistogram((hName+"TTbar"+"qcdselSS"+"Asymm"+subSelName+"Plus").c_str());
  	TH1F *hTTbar22 = get1DHistogram((hName+"TTbar"+"qcdselSS"+"Asymm"+subSelName+"Minus").c_str());
  	hTTbar12->Scale(scale);
  	hTTbar22->Scale(scale);
  	hTTbar1->Add(hTTbar12,1);
  	hTTbar2->Add(hTTbar22,1);
   }

//qcd background

  TH1F * hQCD1 = (TH1F*) getQCDbackgroundAsymm(varName,selName,"Asymm"+subSelName,"Plus");
  TH1F * hQCD2 = (TH1F*) getQCDbackgroundAsymm(varName,selName,"Asymm"+subSelName,"Minus");
  
  if(selName2=="All"){
   	TH1F * hQCD12 = (TH1F*) getQCDbackgroundAsymm(varName,"qcdselSS","Asymm"+subSelName,"Plus");
   	TH1F * hQCD22 = (TH1F*) getQCDbackgroundAsymm(varName,"qcdselSS","Asymm"+subSelName,"Minus");

 	hQCD1->Add(hQCD12,1);
 	hQCD2->Add(hQCD22,1);
   }

  hDYJets1 -> Add(hDYJetsLowM1,1);
  hDYJets2 -> Add(hDYJetsLowM2,1);

// difference histograms (Minus)

  hTTbar1 -> Add(hTTbar2,-1);
  hDYJets1-> Add(hDYJets2,-1);
  hSoup1  -> Add(hSoup2,-1);
  hQCD1 -> Add(hQCD2,-1);

  std::string nazwa44 = "Delta_wazone_"+varName+"_" + selName +"_"+selName2+"_"+subSelName;
  PlotOneHistogram(varName,nazwa44,hSoup1,"delta_wazone");

  hSoup1 -> Add(hTTbar1,-1);
  hSoup1 -> Add(hDYJets1,-1);
  hSoup1 -> Add(hQCD1,-1);

// MC asymmetry

  TH1F *hWJetsP1 = get1D_WJet_Histogram((hName+"WJets"+selName+"Asymm"+subSelName+"Plus").c_str());
  TH1F *hWJetsN1 = get1D_WJet_Histogram((hName+"WJets"+selName+"Asymm"+subSelName+"Minus").c_str());
   sampleName = "WJets";
   weight = getSampleNormalisation(sampleName);
   scale = weight*lumi;
   hWJetsP1->Scale(scale);
   hWJetsN1->Scale(scale);

	 if(selName2=="All"){
	  TH1F *hWJetsP2 = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+"Asymm"+subSelName+"Plus").c_str());
	  TH1F *hWJetsN2 = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+"Asymm"+subSelName+"Minus").c_str());
	  hWJetsP2->Scale(scale);
	  hWJetsN2->Scale(scale);
	  hWJetsN1->Add(hWJetsN2,1);
	  hWJetsP1->Add(hWJetsP2,1);
	  }

  hWJetsP1->Add(hWJetsN1,-1);

// reference W+Jet background

	TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+selName).c_str()); 
	TH1F *hWJets2 = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS").c_str()); 

  std::string WselType = "wselOS";
  if(selName.find("SS")!=std::string::npos) WselType = "wselSS";
  pair<float,float> dataToMCScale = getWNormalisation(WselType,"");
  weight = getSampleNormalisation("WJets");
  scale = weight*lumi*dataToMCScale.first;
  hWJets->Scale(scale);

  if(selName2=="All"){
    WselType = "wselSS";
    dataToMCScale = getWNormalisation(WselType,"");
    float scale = weight*lumi*dataToMCScale.first;
    hWJets2-> Scale(scale);
    hWJets -> Add(hWJets2,1);
  }

// reference W+Jet background, not scaled

	TH1F *hWJetsMC = get1D_WJet_Histogram((hName+"WJets"+selName).c_str()); 
	TH1F *hWJets2MC = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS").c_str()); 

  weight = getSampleNormalisation("WJets");
  scale = weight*lumi;
  hWJetsMC->Scale(scale);

  if(selName2=="All"){
    hWJets2MC-> Scale(scale);
    hWJetsMC -> Add(hWJets2MC,1);
  }

// plots
  std::string nazwa = "WEstimQCD_"+varName+"_" + selName +"_"+selName2+"_"+subSelName;
  PlotTwoHistograms(varName,nazwa,hSoup1,"Asymm Data",hWJets,"MC scaling");

  nazwa = "WEstimMC_"+varName+"_" + selName +"_"+selName2+"_"+subSelName;
  PlotTwoHistograms(varName,nazwa,hWJetsP1,"Asymm MC",hWJetsMC,"MC scaling");

return hSoup1;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::PlotOneHistogram(std::string varName,std::string nazwa, TH1* th1,std::string opis1){

  TCanvas* c = new TCanvas("WEstim","Asymm",460,500);
  TPad* p1 = new TPad("p1","p1",0.06,0.03,0.98,0.98,0); p1->Draw();
  p1->SetLeftMargin(0.1);
  p1->Draw();
  p1->cd();

  float max = th1->GetMaximum();
  th1->SetMaximum(1.3*max);
 
  th1->Draw("hist");
  th1->SetStats(kFALSE);

  TLegend *leg = new TLegend(0.55,0.8,0.8,0.88,NULL,"brNDC");
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

  c->Print(TString::Format("fig_png/%s.png",nazwa.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::PlotTwoHistograms(std::string varName,std::string nazwa, TH1* th1,std::string opis1, TH1* th2, std::string opis2){
// Sometimes parameters should have different values when perfect histograms are needed.

  TCanvas* c = new TCanvas("WEstim","Asymm",460,500);
  TPad* p1 = new TPad("p1","p1",0.05,0.03,0.98,0.98,0); p1->Draw();
  p1->SetLeftMargin(0.12);
  p1->Draw();
  p1->cd();
 
  th1->Draw("hist");
  th1->SetStats(kFALSE);
  th2->SetLineColor(2);
  th2->Draw("same");

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
  	TLegend *leg = new TLegend(0.6,0.63,0.90,0.88,NULL,"brNDC");
  	setupLegend(leg);
  	leg->SetTextSize(0.03);
  	leg->AddEntry(th1,opis1.c_str(),"l");
  	leg->AddEntry(th2,opis2.c_str(),"l");
  	leg->SetHeader(Form("#int N = %.2f",th1int));
  	leg->Draw();
  }

  th1->GetYaxis()->SetTitle("Number of events");
  th2->GetYaxis()->SetTitle("Number of events");
  th1->GetYaxis()->SetTitleOffset(2.1);

  float max = th1->GetMaximum();
  th1->SetMaximum(1.3*max);

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

  c->Print(TString::Format("fig_png/%s.png",nazwa.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::getQCDbackgroundAsymm(std::string varName, std::string selName,std::string asymmSel, std::string subSelName){
// This function prepare qcd background for asymmetry analysis.
// selName = ""/"qcdsel; asymmSel="MT"/"Eta"; subSelName="Plus"/"Minus"

  std::cout<<"Calling method: "<<__func__<<std::endl;
  
  int ifscale=0;
  if(selName.find("SS")==std::string::npos) ifscale=1;
  if(selName.find("qcdsel")!=std::string::npos) selName = "";

  std::string hName = "h1D" + varName;

  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+selName+asymmSel+subSelName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+selName+asymmSel+subSelName).c_str());
  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+selName+asymmSel+subSelName).c_str());
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"qcdselSS"+selName+asymmSel+subSelName).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+"qcdselSS"+selName+asymmSel+subSelName).c_str());

  if(!hSoup) return 0;

  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hSoup->Clone((hName+"hDYJetsLowM"+"qcdselSS"+selName+asymmSel+subSelName).c_str()); hDYJetsLowM->Reset();
  }
  if(!hDYJets){
    hDYJets = (TH1F*)hSoup->Clone((hName+"hDYJets"+"qcdselSS"+selName+asymmSel+subSelName).c_str()); hDYJets->Reset();
  }
  if(!hTTbar){
    hTTbar = (TH1F*)hSoup->Clone((hName+"hTTbar"+"qcdselSS").c_str()); hTTbar->Reset();
  }

  //////////////////////////////////////////////////////////////////////
  float lumi = getLumi();
  std::string sampleName = "DYJetsLowM";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJetsLowM->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "WJets";
  scale = getSampleNormalisation(sampleName)*lumi;
  	hWJets->Scale(getWNormalisation("wselSS",subSelName).first * scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar->Scale(scale);

  hDYJets->Add(hDYJetsLowM,1);

  hSoup->SetName(("h1DQCDEstimate"+varName).c_str());
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);
  hSoup->Add(hWJets,-1);


  ///Clean up the QCD shape, and remove fluctuations around 0 counts.
//  for(unsigned int iBinX=0;iBinX<=hSoup->GetNbinsX();++iBinX){
//    if(hSoup->GetBinContent(iBinX)<3.0) hSoup->SetBinContent(iBinX,0);
//  }
  
  scale =1;
  if(ifscale==1){
	  scale = getQCDOStoSS(selName,subSelName).first;
	  hSoup->Scale(scale);
  }

	stringstream ss;
	ss << scale;
	string str = ss.str();

  TCanvas* c = new TCanvas("WEstim","QCD_OStoSS",460,500);
  std::string nazwa = "QCDBack_weighted_by_asymm"+varName+"_" + selName+"_"+asymmSel+"_"+subSelName+"_scale="+str;
  hSoup->Draw();
  c->Print(TString::Format("fig_png/%s.png",nazwa.c_str()).Data());

  return hSoup;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
