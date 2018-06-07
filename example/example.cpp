#include "CalFR.cpp"

void example(){
  TFile* ifile = new TFile("fakerate_GJet_Cut1.root"); 

  // measured b jet fraction
  TH1F* hfrac1 = (TH1F*)ifile->Get("fraction_GJetData_TypeI");
  TH1F* hfrac2 = (TH1F*)ifile->Get("fraction_GJetData_TypeII");

  // measured fake rate in two samples in GJet Data
  TH1F* hfr1 = (TH1F*)ifile->Get("fakerate_GJetData_TypeI");
  TH1F* hfr2 = (TH1F*)ifile->Get("fakerate_GJetData_TypeII");

  TH1F* hfrb = (TH1F*)hfr1->Clone("hfrb");
  hfrb->SetTitle("fake rate of b jets in GJet Data");
  TH1F* hfrl = (TH1F*)hfr1->Clone("hfrl");
  hfrl->SetTitle("fake rate of light jets in GJet Data");

  FrHistoCal(hfrac1, hfrac2, hfr1, hfr2, hfrb, hfrl);

  TFile *ofile = new TFile("result.root", "RECREATE");
  hfrb->SetDirectory(0);
  hfrb->Write();
  hfrl->SetDirectory(0);
  hfrl->Write();

  std::cout << " -------  Finished :-) ------ " << std::endl;
}
