#include <sstream>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::HTTAnalyzer(const std::string & aName):Analyzer(aName){

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::~HTTAnalyzer(){

  if(myHistos_) delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::initialize(TFileDirectory& aDir,
				 pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  ///The histograms for this analyzer will be saved into "HTTAnalyzer"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTHistograms(&aDir, selectionFlavours_);  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  
  float puWeight = myEventProxy.puWeight;
  float genWeight = myEventProxy.puWeight;
  float eventWeight = puWeight*genWeight;
  
  float lumi = 1.0;
  float scaleFactor = 1.0;
  
  eventWeight *=lumi*scaleFactor;
  
  std::string sampleName = "MC";
  if(myEventProxy.run>1) sampleName = "Data";
  
  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.npv,eventWeight);

  ///Fill SVfit mass
  myHistos_->fill1DHistogram("h1DSVfit"+sampleName,myEventProxy.svfit,eventWeight);
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

