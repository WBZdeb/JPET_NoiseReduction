/// guards
#ifndef HELPER_TOOLS_H
#define HELPER_TOOLS_H

#include <TVector3.h>
#include <TROOT.h>
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx" /// for RNode
using RDF = ROOT::RDataFrame;
/// This was added in root v6.16
using RNode = ROOT::RDF::RNode;


TVector3 calculateAnnihilationPoint(const TVector3& hitA, const TVector3& hitB, double tof);
void saveReportToFile(const std::string& outFileName, ROOT::RDF::RCutFlowReport& report);
RNode applyCuts(const std::vector<std::string>& cuts, RNode& df);



#endif
