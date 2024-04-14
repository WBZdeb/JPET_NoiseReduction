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

void DTW_prep(){
	std::unique_ptr<TCanvas> canv(new TCanvas("canv", "canv", 1920, 1080));
	std::string inFile = "flatTree_40800.unk.evt.root";
	std::string treeName = "FlatTree";

	TChain chain(treeName.c_str());
	chain.Add(inFile.c_str());
	ROOT::RDataFrame df(chain);

	auto columns = df.GetColumnNames();
	for (auto c:columns) {
		cout << c << endl;
		
	}
	//plot number of events with for hit numbers
	auto hist = df.Histo1D({"HitCounts", "HitCounts", 1000, 0, 2000}, "numberOfHits");
	
	//2-hit events, no OPs, no scatters, no prompts
	//auto df_fltr = df.Filter("numberOfHits == 2 && !isOPs && !isScattered && !containsPrompt");
	
	//Accidentals
	//std::cout << "total nb of events:" << *(df_fltr.Count()) << std::endl;	// error: no member named 'Count'	--> what to do here?
	//std::cout << "total nb of accidentals:" << *(df_fltr.Filter("isAcc").Count()) << std::endl;
	
	
	//Energy plot
	//std::unique_ptr<TCanvas> canv2(new TCanvas("canv2", "canv2", 1920, 1080));
	//auto hist2 = df_fltr.Histo1D({"Energy", "Energy", 1000, 0, 2000}, "energy");
}




