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
#include <fstream>
#include <string>
#include <vector>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <cassert>
#include <memory>

#include "calc_lifetime.h"
using namespace std;

void DTW_datainfo(){
	gROOT->Reset();
	
	//open and read file with flatTree paths
	ifstream inFile("flatTrees.txt");
	
	if(!inFile.is_open() ){
		cerr << "Error opening the file! Are you sure file 'flatTrees.txt' exists?" << endl;
	}	
	
	std::string filePath;
	std::string treeName = "FlatTree";
	
	TChain chain(treeName.c_str());
	
	while(getline(inFile, filePath)){
		chain.Add(filePath.c_str());
	}
	
	ROOT::RDataFrame df(chain);
	
	//close file
	inFile.close();
	
	
	//========================
	//	Info about data
	//========================
	
	// Ratio of Acc/All
	std::cout << "isAcc/All before filters" << std::endl;
	std::cout << "Acc: " << (*(df.Filter("isAcc").Take<int>("eventNumber"))).size() << std::endl;
	std::cout << "All: " << (*(df.Take<int>("eventNumber"))).size() << std::endl;
	
	std::cout << 100.0*(*(df.Filter("isAcc").Take<int>("eventNumber"))).size() / (*(df.Take<int>("eventNumber"))).size() << std::endl;
	std::cout << std::endl;
	
	
	// Ratio of Acc/All first cut
	auto df_f1 = df.Filter("(numberOfHits == 3) && (containsPrompt)");
	
	std::cout << "isAcc/All first cut" << std::endl;
	std::cout << "Acc: " << (*(df_f1.Filter("isAcc").Take<int>("eventNumber"))).size() << std::endl;
	std::cout << "All: " << (*(df_f1.Take<int>("eventNumber"))).size() << std::endl;
	
	std::cout << 100.0*(*(df_f1.Filter("isAcc").Take<int>("eventNumber"))).size() / (*(df_f1.Take<int>("eventNumber"))).size() << std::endl;
	std::cout << std::endl;
	
	
	// Ratio of Acc/All second cut
	auto df_f2 = df.Filter("(numberOfHits == 3) && (isPickOff) && (containsPrompt)");
	
	std::cout << "isAcc/All Dsecond cut" << std::endl;
	std::cout << "Acc: " << (*(df_f2.Filter("isAcc").Take<int>("eventNumber"))).size() << std::endl;
	std::cout << "All: " << (*(df_f2.Take<int>("eventNumber"))).size() << std::endl;
	
	std::cout << 100.0*(*(df_f2.Filter("isAcc").Take<int>("eventNumber"))).size() / (*(df_f2.Take<int>("eventNumber"))).size() << std::endl;
	std::cout << std::endl;
	
	// Ratio of Acc/All after all filters
	auto df_fltr = df.Filter("(numberOfHits == 3) && (isPickOff) && (!isScattered) && (!isSecondary) && (containsPrompt)");
	
	std::cout << "isAcc/All after all filters" << std::endl;
	std::cout << "Acc: " << (*(df_fltr.Filter("isAcc").Take<int>("eventNumber"))).size() << std::endl;
	std::cout << "All: " << (*(df_fltr.Take<int>("eventNumber"))).size() << std::endl;
	
	std::cout << 100.0*(*(df_fltr.Filter("isAcc").Take<int>("eventNumber"))).size() / (*(df_fltr.Take<int>("eventNumber"))).size() << std::endl;
	std::cout << std::endl;
	
	//Random coincidences
	std::cout << "\n Type 1 DTW:" << std::endl;
	DTW_type1(df, 2);
	std::cout << "\n Type 2 DTW:" << std::endl;
	DTW_type2(df, 2);
	std::cout << "\n Type 3 DTW:" << std::endl;
	DTW_type3(df, 2);
	std::cout << "\n Type 4 DTW:" << std::endl;
	DTW_type4(df, 2);

	
	//Registered source activity
	auto window_num = *(df.Take<int>("timeWindowNumber"));
	auto window_count = window_num.back();
	double runtime = 0.00005 * window_count;	//[s]
	std::cout << "Registered source activity: "	<< 1022681.0/runtime << " Bq" << std::endl;
	
	
	//========================
	//	Test for p-ps
	//========================
	
	std::cout << "\n Test p-ps" << std::endl;
	std::cout << "Flagi:	2 hits, !containsPrompt, !isScattered, !isOPs" << "\nDTW: " << std::endl;
	
	auto df_pps = df.Filter("(numberOfHits == 2) && (!isScattered) && (!containsPrompt) && (!isOPs)");
	DTW_type1(df_pps, 2);
	std::cout << "isAcc: " << (*(df_pps.Filter("isAcc").Take<int>("eventNumber"))).size() << std::endl;
	
	std::cout << "\n Po fladze !isPickOff" << std::endl;
	std::cout << "All: " << (*(df_pps.Filter("(!isPickOff)").Take<int>("eventNumber"))).size() << std::endl;
	//DTW_type1(df_pps.Filter("!isPickOff"), 2);
	std::cout << "isAcc: " << (*(df_pps.Filter("(!isPickOff) && (isAcc)").Take<int>("eventNumber"))).size() << std::endl;
	
	//========================
	//	Histograms
	//========================
	
	ofstream hitsFile("Filter_hist/Hit_count.txt");
	
	if(!hitsFile.is_open() ){
		cerr << "Error opening the file! Are you sure file 'Filter_hist/Hit_count.txt' exists?" << endl;
	}
	
	std::vector<const char*> saveLocation = {"Filter_hist/isPickOff.png", "Filter_hist/Cut_1.png", "Filter_hist/Cut_2.png", "Filter_hist/Cut_3.png", "Filter_hist/Cut_4.png"};
	std::vector<std::string> filters = {"(isPickOff)", "(!isScattered)", "(!isSecondary)", "(containsPrompt)", "(numberOfHits == 3)"};
	auto df_filt = df.Filter("isPickOff");
	
	hitsFile << "Event count:		" << (*(df.Take<int>("eventNumber"))).size() << std::endl;
	
	for(int i = 0; i < filters.size(); i++){
		auto filt = filters[i];
		df_filt = df_filt.Filter(filt);
		hitsFile << "After " << filt << ": " << (*(df_filt.Take<int>("eventNumber"))).size() << std::endl;
		
		//histogram
		std::unique_ptr<TCanvas> histo(new TCanvas("Hit multiplicity", "Hit multiplicity", 1920, 1080));
		auto hist = df_filt.Histo1D({"numberOfHits", "numberOfHits", 7, -0, 7}, "numberOfHits");
		hist->GetXaxis()->SetTitle("numberOfHits");
		hist->Draw();
		histo->SaveAs(saveLocation[i]);
	}
	
	hitsFile.close();
}

int main()
{
  DTW_datainfo();
  return 0;  
}

