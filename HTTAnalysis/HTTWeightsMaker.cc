#include <sstream>

#include "HTTWeightsMaker.h"
#include "HTTWeightHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTWeightsMaker::HTTWeightsMaker(const std::string & aName):Analyzer(aName){

  ///Load ROOT file with PU histograms.
  std::string filePath = "Data_Pileup_2015D_Feb02.root";
  puDataFile_ = new TFile(filePath.c_str());

  filePath = "MC_Spring15_PU25_Startup.root";
  puMCFile_ = new TFile(filePath.c_str());


  ntupleFile_ = 0; 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTWeightsMaker::~HTTWeightsMaker(){

  if(myHistos_) delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTWeightsMaker::clone() const{

  HTTWeightsMaker* clone = new HTTWeightsMaker(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::initialize(TFileDirectory& aDir,
			     pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  myHistos_ = new HTTWeightHistograms(&aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::getPreselectionEff(const EventProxyHTT & myEventProxy){

    TH1F *hStatsFromFile = (TH1F*)myEventProxy.getTTree()->GetCurrentFile()->Get("m2n/hStats");

    std::string hName = "h1DStats"+getSampleName(myEventProxy);
    TH1F *hStats = myHistos_->get1DHistogram(hName.c_str(),true);
    hStats->SetBinContent(2,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(1)));   
    hStats->SetBinContent(3,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(3)));
    delete hStatsFromFile;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTWeightsMaker::getSampleName(const EventProxyHTT & myEventProxy){

  if(myEventProxy.wevent->sample()==0) return "Data";
  if(myEventProxy.wevent->sample()==1){
    int decayModeBoson = myEventProxy.wevent->decayModeBoson();
    if(decayModeBoson==7) return "DYJetsMuMu";
    else if(decayModeBoson==6) return "DYJetsEE";
    else if(decayModeBoson==0) return "DYJetsMuTau";
    else return "DYJetsOther";    
  }
  if(myEventProxy.wevent->sample()==11) return "DYJetsLowM";
  if(myEventProxy.wevent->sample()==2){
    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    if(fileName.find("WJetsHT100to200")!=std::string::npos) return "WJetsHT100to200";
    if(fileName.find("WJetsHT200to400")!=std::string::npos) return "WJetsHT200to400";
    if(fileName.find("WJetsHT400to600")!=std::string::npos) return "WJetsHT400to600";
    if(fileName.find("WJetsHT600toInf")!=std::string::npos) return "WJetsHT600toInf";
    return "WJets";
  }
  if(myEventProxy.wevent->sample()==3) return "TTbar";
  if(myEventProxy.wevent->sample()==5) return "H";
  if(myEventProxy.wevent->sample()==6) return "A";

  return "Unknown";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTWeightsMaker::getPUWeight(const EventProxyHTT & myEventProxy){

  ///Load histogram only once,later fetch it from vector<TH1F*>
  ///At the same time divide the histogram to get the weight.
  ///First load Data PU
  if(!hPUVec_.size())  hPUVec_.resize(64);

  if(!hPUVec_[myEventProxy.wevent->sample()]){
    std::string hName = "pileup";
    TH1F *hPUData = (TH1F*)puDataFile_->Get(hName.c_str());
    TH1F *hPUSample = (TH1F*)puMCFile_->Get(hName.c_str());
    ///Normalise both histograms.
    hPUData->Scale(1.0/hPUData->Integral(0,hPUData->GetNbinsX()+1));
    hPUSample->Scale(1.0/hPUSample->Integral(0,hPUSample->GetNbinsX()+1));
    ///
    hPUData->SetDirectory(0);
    hPUSample->SetDirectory(0);
    hPUData->Divide(hPUSample);
    hPUData->SetName(("h1DPUWeight"+getSampleName(myEventProxy)).c_str());
    ///To get uniform treatment put weight=1.0 for under/overlow bins of
    ///data PU, as nPU for data has a dummy value.
    if(getSampleName(myEventProxy)=="Data"){
      hPUData->SetBinContent(0,1.0);
      hPUData->SetBinContent(hPUData->GetNbinsX()+1,1.0);
    }
    hPUVec_[myEventProxy.wevent->sample()] =  hPUData;
  }

  int iBinPU = hPUVec_[myEventProxy.wevent->sample()]->FindBin(myEventProxy.wevent->npu());
  return  hPUVec_[myEventProxy.wevent->sample()]->GetBinContent(iBinPU);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTWeightsMaker::getGenWeight(const EventProxyHTT & myEventProxy){

  if(myEventProxy.wevent->sample()==0) return 1.0;
  ///generator weight broken in miniAODv2
  /*
  if(myEventProxy.wevent->sample()==1) return myEventProxy.wevent->genevtweight()/23443.423;  
  if(myEventProxy.wevent->sample()==2) return myEventProxy.wevent->genevtweight()/225892.45;  
  if(myEventProxy.wevent->sample()==3) return myEventProxy.wevent->genevtweight()/6383;
  */
  /*
  if(myEventProxy.wevent->sample()==2){
    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    if(fileName.find("WJetsHT250to200")!=std::string::npos) return 1345/(3*20508.9+1345+359.7+48.9+18.77);
    if(fileName.find("WJetsHT200to400")!=std::string::npos) return 359.7/(3*20508.9+1345+359.7+48.9+18.77);
    if(fileName.find("WJetsHT400to600")!=std::string::npos) return 48.9/(3*20508.9+1345+359.7+48.9+18.77);
    if(fileName.find("WJetsHT600toInf")!=std::string::npos) return 18.7/(3*20508.9+1345+359.7+48.9+18.77);
    else return 3*20508.9/(3*20508.9+1345+359.7+48.9+18.77);
  }
  */
  return 1;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::fillControlHistos(float eventWeight, std::string & hNameSuffix){
 
  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aMuon.mt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaMuon"+hNameSuffix,aMuon.eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix,aMuon.iso(),eventWeight);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<std::string> HTTWeightsMaker::getTauDecayName(int decModeMinus, int decModePlus){

  std::vector<std::string> types;

  if(decModeMinus==tauDecay1ChargedPion0PiZero && decModePlus==tauDecay1ChargedPion0PiZero) types.push_back("PiPi0Pi0");

  if(isOneProng(decModeMinus) && isOneProng(decModePlus) ) types.push_back("1Prong1Prong");

  if( (decModeMinus==tauDecay1ChargedPion0PiZero && isLepton(decModePlus) ) ||
      (isLepton(decModeMinus) && decModePlus==tauDecay1ChargedPion0PiZero)) types.push_back("Lepton1Prong0Pi0");
    
  if( (isOneProng(decModeMinus) && isLepton(decModePlus) ) ||
      ( isLepton(decModeMinus) && isOneProng(decModePlus) ) ) types.push_back("Lepton1Prong");

  if(decModeMinus==tauDecay1ChargedPion1PiZero && decModePlus==tauDecay1ChargedPion1PiZero ) types.push_back("PiPlusPiMinus2Pi0");


  if( isOneProng(decModeMinus) && decModeMinus!=tauDecay1ChargedPion0PiZero && 
      isOneProng(decModePlus) && decModePlus!=tauDecay1ChargedPion0PiZero )   types.push_back("1Prong1ProngXPi0");

  if(isLepton(decModePlus) && isLepton(decModeMinus)) types.push_back("LeptonLepton");

  return types;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTWeightsMaker::isOneProng(int decMode){
  if(decMode==tauDecay1ChargedPion0PiZero ||
     decMode==tauDecay1ChargedPion1PiZero ||
     decMode==tauDecay1ChargedPion2PiZero ||
     decMode==tauDecay1ChargedPion3PiZero ) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTWeightsMaker::isLepton(int decMode){
  if(decMode==tauDecaysElectron || decMode==tauDecayMuon) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<Wjet> HTTWeightsMaker::getSeparatedJets(const EventProxyHTT & myEventProxy, const Wtau & aTau,
						const Wmu & aMuon, float deltaR){

  std::vector<Wjet> separatedJets;
  
  for(auto aJet: *myEventProxy.wjet){
    float dRTau = sqrt(pow(aJet.eta() - aTau.eta(),2) + pow(aJet.phi() - aTau.phi(),2));
    float dRMu = sqrt(pow(aJet.eta() - aMuon.eta(),2) + pow(aMuon.phi() - aTau.phi(),2));
    if(dRTau>deltaR && dRMu>deltaR) separatedJets.push_back(aJet);
  }

  return separatedJets;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::setAnalysisObjects(const EventProxyHTT & myEventProxy){

  aEvent = *myEventProxy.wevent;  
  aPair = (*myEventProxy.wpair)[0];
  aTau = (*myEventProxy.wtau)[0];
  if(myEventProxy.wtauGen && myEventProxy.wtauGen->size()){
    aGenNegativeTau = (*myEventProxy.wtauGen)[0];
    aGenPositiveTau = (*myEventProxy.wtauGen)[1];
  }
  aMuon = (*myEventProxy.wmu)[0];
  aMET = (*myEventProxy.wmet)[0];
  aSeparatedJets = getSeparatedJets(myEventProxy, aTau, aMuon, 0.5);
  aJet = aSeparatedJets.size() ? aSeparatedJets[0] : Wjet();
  nJets30 = count_if(aSeparatedJets.begin(), aSeparatedJets.end(),[](const Wjet & aJet){return aJet.pt()>30;});

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<bool, bool> HTTWeightsMaker::checkTauDecayMode(const EventProxyHTT & myEventProxy){

  bool goodGenDecayMode = false;
  bool goodRecoDecayMode = false;
  std::vector<std::string> decayNamesGen = getTauDecayName(myEventProxy.wevent->decModeMinus(), myEventProxy.wevent->decModePlus());
  std::vector<std::string> decayNamesReco = getTauDecayName(aTau.decayMode(), tauDecayMuon);
  for(auto it: decayNamesGen) if(it.find("Lepton1Prong0Pi0")!=std::string::npos) goodGenDecayMode = true;
  for(auto it: decayNamesReco) if(it.find("Lepton1Prong0Pi0")!=std::string::npos) goodRecoDecayMode = true;

  return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::addBranch(TTree *tree){

  tree->Branch("muonPt",&muonPt);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTWeightsMaker::analyze(const EventProxyBase& iEvent){

  //clearTTreeVariables();

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

  std::string sampleName = getSampleName(myEventProxy);
  std::string hNameSuffix = sampleName;
  float puWeight = getPUWeight(myEventProxy);
  float genWeight = getGenWeight(myEventProxy);
  float eventWeight = puWeight*genWeight;

  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,0,eventWeight);
  getPreselectionEff(myEventProxy);
  /////////////////////////////////////////////////////////////////////////////
  if(!myEventProxy.wpair->size() || !myEventProxy.wtau->size() || !myEventProxy.wmu->size()) return true;

  setAnalysisObjects(myEventProxy);

  std::pair<bool, bool> goodDecayModes = checkTauDecayMode(myEventProxy);
  bool goodGenDecayMode = goodDecayModes.first;
  bool goodRecoDecayMode = goodDecayModes.second;

//h1DPhi_nVectorsWJetsGenNoOfflineSel

  ///This stands for core selection, that is common to all regions.
  bool tauKinematics = aTau.pt()>30 && fabs(aTau.eta())<2.3;
  bool tauID = 1; // for 16.02. data
//  bool tauID = aTau.tauID(byMediumCombinedIsolationDeltaBetaCorr3Hits);
  bool muonKinematics = aMuon.pt()>19 && fabs(aMuon.eta())<2.1;
  bool trigger = aPair.trigger(HLT_IsoMu17_eta2p1);
  if(sampleName=="Data") trigger = aPair.trigger(HLT_IsoMu18);
  bool extraRequirements = aTau.decayMode()!=5 && aTau.decayMode()!=6 && nJets30==0;

  if(!myEventProxy.wpair->size()) return true;
  if(!tauKinematics || !muonKinematics || !trigger) return true;
  //if(!tauKinematics || !tauID || !muonKinematics || !trigger) return true;
  //if(!extraRequirements) return true;

  ///Note: parts of the signal/control region selection are applied in the following code.
  ///FIXME AK: this should be made in a more clear way.
  bool baselineSelection = aPair.diq()==-1 && aMuon.mt()<250 && aMuon.iso()<0.1;
  bool wSelection = aMuon.mt()>60 && aMuon.iso()<0.1;
  bool qcdSelectionSS = aPair.diq()==1;
  bool qcdSelectionOS = aPair.diq()==-1;

  std::string hNameSuffixP = sampleName+"Plus";
  std::string hNameSuffixN = sampleName+"Minus";

  ///Fill variables stored in TTree
  muonPt = aMuon.pt();

  ///Histograms for the baseline selection  
  if(baselineSelection){
    fillControlHistos(eventWeight, hNameSuffix);
	if(aMuon.charge()==1) fillControlHistos(eventWeight, hNameSuffixP);
	if(aMuon.charge()==-1) fillControlHistos(eventWeight, hNameSuffixN);
  }

  // to create global QCD histograms; it is not important what is inside
  std::string hNameSuffixQP= "QCDDiffPlus"; // in this context it is a sum
  std::string hNameSuffixQN= "QCDDiffMinus"; // in this context it is a difference

  //Fill empty histos and "Plus" "Minus" control histos
    if(aMuon.charge()==1 && aMuon.mt()<250 && aMuon.iso()<0.1)  {
      fillControlHistos(eventWeight, hNameSuffixQP);
      hNameSuffixP=sampleName+"DiffPlus";
      fillControlHistos(eventWeight, hNameSuffixP);}
//      myHistos_->fill1DHistogram("h1DMassTransDiff"+hNameSuffixP,aMuon.mt(),eventWeight);}
    if(aMuon.charge()==-1 && aMuon.mt()<250 && aMuon.iso()<0.1) {
      hNameSuffixN=sampleName+"DiffMinus";
      fillControlHistos(eventWeight, hNameSuffixQN);
      fillControlHistos(eventWeight, hNameSuffixN);}
//      myHistos_->fill1DHistogram("h1DMassTransDiff"+hNameSuffixN,aMuon.mt(),eventWeight);}

  ///Histograms for the QCD control region
  if(qcdSelectionSS){
    hNameSuffix = sampleName+"qcdselSS";
    hNameSuffixP = sampleName+"qcdselSSPlus";
    hNameSuffixN = sampleName+"qcdselSSMinus";
    ///SS ans OS isolation histograms are filled only for mT<40 to remove possible contamnation
    //from TT in high mT region.
    if(aMuon.mt()<250) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.iso(),eventWeight);
    if(aMuon.mt()<250 && aMuon.charge()==1) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix+"Plus",aMuon.iso(),eventWeight);
    if(aMuon.mt()<250 && aMuon.charge()==-1) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix+"Minus",aMuon.iso(),eventWeight);
    ///Fill SS histos in signal mu isolation region. Those histograms
    ///provide shapes for QCD estimate in signal region and in various control regions.
    ///If control region has OS we still use SS QCD estimate.
    if(aMuon.mt()<250 && aMuon.iso()<0.1) fillControlHistos(eventWeight, hNameSuffix);
    if(aMuon.mt()<250 && aMuon.iso()<0.1 && aMuon.charge()==1) fillControlHistos(eventWeight, hNameSuffixP);
    if(aMuon.mt()<250 && aMuon.iso()<0.1 && aMuon.charge()==-1) fillControlHistos(eventWeight, hNameSuffixN);
    if(aMuon.mt()>60){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselSS",aMuon.mt(),eventWeight);    
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselOS",aMuon.mt(),eventWeight);   
      if(aMuon.charge()==1)  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselSS"+"Plus",aMuon.mt(),eventWeight); 
      if(aMuon.charge()==-1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselSS"+"Minus",aMuon.mt(),eventWeight);
      if(aMuon.charge()==1)  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselOS"+"Plus",aMuon.mt(),eventWeight); 
      if(aMuon.charge()==-1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselOS"+"Minus",aMuon.mt(),eventWeight); 
    }
  }

  ///Make QCD shape histograms for specific selection.
  ///Using the same SS/OS scaling factor for now.    
  if(qcdSelectionOS){
    hNameSuffix = sampleName+"qcdselOS";
    hNameSuffixP = sampleName+"qcdselOSPlus";
    hNameSuffixN = sampleName+"qcdselOSMinus";    
    if(aMuon.mt()<250) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.iso(),eventWeight);
    if(aMuon.mt()<250 && aMuon.charge()==1) {
	myHistos_->fill1DHistogram("h1DIso"+hNameSuffix+"Plus",aMuon.iso(),eventWeight);
	 fillControlHistos(eventWeight, hNameSuffixP);}
    if(aMuon.mt()<250 && aMuon.charge()==-1){
	myHistos_->fill1DHistogram("h1DIso"+hNameSuffix+"Minus",aMuon.iso(),eventWeight);
	 fillControlHistos(eventWeight, hNameSuffixN); }
  }

  ///Histograms for the WJet control region. 
  if(wSelection){
    hNameSuffix = sampleName+"wsel";
    if(aPair.diq()==-1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OS",aMuon.mt(),eventWeight);
    if(aPair.diq()== 1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SS",aMuon.mt(),eventWeight);
	if(aMuon.charge()==1){
    if(aPair.diq()==-1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OSPlus",aMuon.mt(),eventWeight);
    if(aPair.diq()== 1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SSPlus",aMuon.mt(),eventWeight);
	}
	if(aMuon.charge()==-1){
    if(aPair.diq()==-1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OSMinus",aMuon.mt(),eventWeight);
    if(aPair.diq()== 1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SSMinus",aMuon.mt(),eventWeight);
	}
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
