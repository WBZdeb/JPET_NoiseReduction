#include <TROOT.h>
#include <TChain.h>
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif
#include "ROOT/RDataFrame.hxx"
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif

#include <iostream>
#include <string>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <cassert>
#include <memory>
using namespace std;

void load_data(const std::string inFile ="")
{
  std::unique_ptr<TCanvas> canv(new TCanvas("canv", "canv", 1920, 1080));

  string treeName = "FlatTree";
  //TFile file(inFile.c_str(), "READ");
  //auto tree = (TTree*) (file.Get(treeName.c_str()));
  //assert(tree);
  //tree->Scan();
  TChain chain(treeName.c_str());
  chain.Add(inFile.c_str());
  ROOT::RDataFrame df(chain);
  auto columns = df.GetColumnNames();
  for (auto c:columns) {
    std::cout << c << std::endl;
    
  }
  auto df_selected = df.Filter("containsPrompt");
  std::cout << "total nb of events:" << *(df.Count()) << std::endl;
  std::cout << "nb of events after selection:" << *(df_selected.Count()) << std::endl;
  //Define("Hour", getUTCTimeInHours, {"GpsTime"});
  //Define("radius", "TMath::Sqrt(x*x + y*y)");
  auto hist = df.Histo1D({"energy", ";energy [?]; nevents", 1000, 0, 2000}, "energy");
  auto histSelected = df_selected.Histo1D({"energySelected", ";energy [?]; nevents", 1000, 0, 2000}, "energy");
  //hist->SetCanExtend(TH1::kAllAxes);

  hist->Draw();
  histSelected->Draw("same");
  canv->SaveAs("bla.png");
}

int main()
{
  load_data("YOUR_FILE"); 
  return 0;
}
