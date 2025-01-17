#include <sstream>
#include <bitset>

#include "RooRealVar.h"

#include "HTTAnalyzerTT.h"
#include "HTTHistogramsTT.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::getPreselectionEff(const EventProxyHTT & myEventProxy){

  if(true || ntupleFile_!=myEventProxy.getTTree()->GetCurrentFile()){
    ntupleFile_ = myEventProxy.getTTree()->GetCurrentFile();
    if(hStatsFromFile) delete hStatsFromFile;
    hStatsFromFile = (TH1F*)ntupleFile_->Get("hStats");
    
    std::string hName = "h1DStats"+getSampleName(myEventProxy);
    TH1F *hStats = myHistos_->get1DHistogram(hName.c_str(),true);
    
    float genWeight = getGenWeight(myEventProxy);
    
    hStats->SetBinContent(2,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(1))*genWeight);   
    hStats->SetBinContent(3,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(3))*genWeight);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzerTT::getSampleName(const EventProxyHTT & myEventProxy){

  return getSampleNameFromFileName(myEventProxy);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzerTT::getSampleNameFromFileName(const EventProxyHTT & myEventProxy){

  std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

  if(fileName.find("W1JetsToLNu")!=std::string::npos) return "W1Jets";
  else if(fileName.find("W2JetsToLNu")!=std::string::npos) return "W2Jets";
  else if(fileName.find("W3JetsToLNu")!=std::string::npos) return "W3Jets";
  else if(fileName.find("W4JetsToLNu")!=std::string::npos) return "W4Jets";
  else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) return "W0Jets";
  else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) return "WAllJets";
  //TEST else if(fileName.find("WJetsToLNu")!=std::string::npos) return "WJets";
  
  else if(fileName.find("Run201")!=std::string::npos) return "Data";
  else if(fileName.find("SUSYGluGluToHToTauTau")!=std::string::npos) return "A";
  else if(fileName.find("GluGluHToTauTauM120")!=std::string::npos) return "ggH120";
  else if(fileName.find("VBFHToTauTauM120")!=std::string::npos) return "qqH120";  
  else if(fileName.find("GluGluHToTauTauM125")!=std::string::npos) return "ggH125";
  else if(fileName.find("VBFHToTauTauM125")!=std::string::npos) return "qqH125";
  else if(fileName.find("GluGluHToTauTauM130")!=std::string::npos) return "ggH130";
  else if(fileName.find("VBFHToTauTauM130")!=std::string::npos) return "qqH130";
  else if(fileName.find("GluGluHToTauTauM135")!=std::string::npos) return "ggH135";
  else if(fileName.find("VBFHToTauTauM135")!=std::string::npos) return "qqH135";
  else if(fileName.find("STtWantitop")!=std::string::npos) return "Wantitop";
  else if(fileName.find("STtWtop")!=std::string::npos) return "Wtop";
  else if(fileName.find("STtchannel__antitop")!=std::string::npos) return "t-channel_antitop";
  else if(fileName.find("STtchannel__top")!=std::string::npos) return "t-channel_top";   
  else if(fileName.find("ZZTo2L2Q")!=std::string::npos) return "ZZTo2L2Q";
  else if(fileName.find("ZZTo4L")!=std::string::npos) return "ZZTo4L";
  else if(fileName.find("WZTo1L3Nu")!=std::string::npos) return "WZTo1L3Nu";
  else if(fileName.find("WZJToLLLNu")!=std::string::npos) return "WZJToLLLNu";
  else if(fileName.find("WWTo1L1Nu2Q")!=std::string::npos) return "WWTo1L1Nu2Q";
  else if(fileName.find("WZTo1L1Nu2Q")!=std::string::npos) return "WZTo1L1Nu2Q";
  else if(fileName.find("VVTo2L2Nu")!=std::string::npos) return "VVTo2L2Nu";
  else if(fileName.find("WZTo2L2Q")!=std::string::npos) return "WZTo2L2Q";
  else if(fileName.find("QCD")!=std::string::npos) return "QCD_MC";
  else if(fileName.find("DY")!=std::string::npos) return getDYSampleName(myEventProxy);
  else if(fileName.find("TTTune")!=std::string::npos) return "TTbar";
  
  return "Unknown";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzerTT::getDYSampleName(const EventProxyHTT & myEventProxy){

  std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

  std::string jetsName = "";    
  if(fileName.find("DY1JetsToLL")!=std::string::npos) jetsName ="1Jets";
  else if(fileName.find("DY2JetsToLL")!=std::string::npos) jetsName = "2Jets";
  else if(fileName.find("DY3JetsToLL")!=std::string::npos) jetsName = "3Jets";
  else if(fileName.find("DY4JetsToLL")!=std::string::npos) jetsName = "4Jets";
  else if(fileName.find("DYJetsToLLM10to50")!=std::string::npos) jetsName =  "LowM";
  //TEST else if(fileName.find("DYJetsToLLM50")!=std::string::npos) jetsName =  "Jets";  
  else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) jetsName =  "0Jets";
  else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) jetsName =  "AllJets";  
  
  int decayModeBoson = myEventProxy.event->getDecayModeBoson();
  int tauMCMatch1 = 6;
  int tauMCMatch2 = 6;
  if(myEventProxy.pairs->size()){
    HTTPair aPair = (*myEventProxy.pairs)[0];
    HTTParticle aTau1 = aPair.getLeg1();
    HTTParticle aTau2 = aPair.getLeg2();
    tauMCMatch1 = aTau1.getProperty(PropertyEnum::mc_match);
    tauMCMatch2 = aTau1.getProperty(PropertyEnum::mc_match);
  }    
  std::string decayName = "Unknown";
  if(tauMCMatch1==5 && tauMCMatch2==5) decayName = "ZTT";
  else if(tauMCMatch1<6 && tauMCMatch2<6) decayName = "ZL";
  else decayName = "ZJ";

  return "DY"+decayName+jetsName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzerTT::getPUWeight(const EventProxyHTT & myEventProxy){

  if(getSampleName(myEventProxy)=="Data") return 1.0;

  if(!puDataFile_ || !puMCFile_ || 
     puDataFile_->IsZombie() || 
     puMCFile_->IsZombie()){ return 1.0; }

  if(!hPUVec_.size())  hPUVec_.resize(1);

  if(!hPUVec_[0]){
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
    hPUVec_[0] =  hPUData;
  }

  int iBinPU = hPUVec_[0]->FindBin(myEventProxy.event->getNPU());  
  return hPUVec_[0]->GetBinContent(iBinPU);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzerTT::getGenWeight(const EventProxyHTT & myEventProxy){

  ///MC weights cab be quite large, but are always +-const.
  ///to avoid counter overflow we keep only sign.
  return myEventProxy.event->getMCWeight()/std::abs(myEventProxy.event->getMCWeight());
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<std::string> HTTAnalyzerTT::getTauDecayName(int decModeMinus, int decModePlus){

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
bool HTTAnalyzerTT::isOneProng(int decMode){
  if(decMode==tauDecay1ChargedPion0PiZero ||
     decMode==tauDecay1ChargedPion1PiZero ||
     decMode==tauDecay1ChargedPion2PiZero ||
     decMode==tauDecay1ChargedPion3PiZero ) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzerTT::isLepton(int decMode){
  if(decMode==tauDecaysElectron || decMode==tauDecayMuon) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzerTT::getLeptonCorrection(float eta, float pt, hadronicTauDecayModes tauDecayMode){

  if(sampleName.find("Data")!=std::string::npos) return 1.0;
  
  if(tauDecayMode == tauDecayMuon){
    scaleWorkspace->var("m_pt")->setVal(pt);
    scaleWorkspace->var("m_eta")->setVal(eta);
    float muon_id_scalefactor = scaleWorkspace->function("m_id_ratio")->getVal();
    float muon_iso_scalefactor = scaleWorkspace->function("m_iso_ratio")->getVal();
    float muon_trg_efficiency = scaleWorkspace->function("m_trgOR_data")->getVal();//OR of the HLT_IsoMu22 and HLT_IsoTkMu22
    return  muon_id_scalefactor*muon_iso_scalefactor*muon_trg_efficiency;
  }
  else if(tauDecayMode == tauDecaysElectron) return 1.0;  
  else{
    float tau_id_scalefactor = 1.0;
    float tau_trg_efficiency = 1.0;
    scaleWorkspace->var("t_pt")->setVal(pt);
    scaleWorkspace->var("t_eta")->setVal(eta);
    scaleWorkspace->var("t_dm")->setVal(tauDecayMode);
    if(sampleName.find("H")==std::string::npos &&
       sampleName.find("ZTT")==std::string::npos
       ) {
      tau_trg_efficiency = scaleWorkspace->function("t_trgTightIsoSS_data")->getVal();      
    }
    else{
      tau_id_scalefactor = scaleWorkspace->function("t_iso_mva_t_pt40_eta2p1_sf")->getVal();
      //tau_id_scalefactor = 0.9;//according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#MC_corrections
      tau_trg_efficiency = scaleWorkspace->function("t_trgTightIso_data")->getVal();
    }
    return tau_id_scalefactor*tau_trg_efficiency;
  }
  return 1.0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
