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
	gROOT->Reset();

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
	//plot number of events for given hit numbers
	auto hist = df.Histo1D({"HitCounts", "HitCounts", 6, 0, 6}, "numberOfHits");
	hist->GetXaxis()->SetTitle("Hit number");
	hist->GetYaxis()->SetTitle("Event count");
	hist->Draw();
	canv->SaveAs("histEventCount.png");
	
	//2-hit events, no OPs, no scatters, no prompts
	auto df_fltr = df.Filter("numberOfHits == 2 && !isOPs && !isScattered && !containsPrompt");
	
	//Accidentals
	std::cout << "total nb of events:" << *(df_fltr.Count())/2 << std::endl;
	std::cout << "total nb of accidentals:" << *(df_fltr.Filter("isAcc").Count())/2 << std::endl;
	
	
	//Energy plot
	std::unique_ptr<TCanvas> canv2(new TCanvas("canv2", "canv2", 1920, 1080));
	auto hist2 = df_fltr.Histo1D({"Energy", "Energy", 1000, 0, 1000}, "energy");
	hist2->GetXaxis()->SetTitle("Energy [keV]");
	hist2->GetYaxis()->SetTitle("Event count");
	hist2->Draw();
	canv2->SaveAs("histEnergy.png");
	
	//test event numbers
	/*
	auto disp = df_fltr.Range(0, 50).Display({"numberOfHits", "eventNumber", "time"});
	disp->Print();
	
	auto eventNum = *(df_fltr.Take<int>("eventNumber"));
	for (int ind = 0; ind < 50; ind++){
		std::cout << eventNum[ind] << std::endl;
	}*/
	
	//Plot time difference between first and second hit
	auto time_Vec = *(df_fltr.Take<vector<float>>("time"));
	
	auto timeDiff = new TCanvas("timeDiff", "timeDiff", 0, 0, 1920, 1080);
	TH1F *histTD = new TH1F("t1 - t2","t1 - t2",150,0,1500);
	
	for (int i = 0; i < time_Vec.size(); i++){
		float diff = time_Vec[i][1] - time_Vec[i][0];
		if (diff < 0){ diff = -diff; }
		histTD->Fill(diff);
	}
	histTD->GetXaxis()->SetTitle("T1 - T2");
	histTD->GetYaxis()->SetTitle("Event count");
	histTD->Draw();
	timeDiff->SaveAs("histT1-T2.png");

}




