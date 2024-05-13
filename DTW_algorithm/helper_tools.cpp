#include <TVector3.h>

#include <TROOT.h>
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx" /// for RNode

const double kLightVelocity_cm_ps = 0.0299792458;

using RDF = ROOT::RDataFrame;

TVector3 calculateAnnihilationPoint(const TVector3& hitA, const TVector3& hitB, double tof)
{
  TVector3 middleOfLOR = 0.5 * (hitA + hitB);
  TVector3 versorOnLOR = (hitB - hitA).Unit()  ;

  double shift = 0.5 * tof  * kLightVelocity_cm_ps;
  TVector3 annihilationPoint(middleOfLOR.X() + shift * versorOnLOR.X(),
                             middleOfLOR.Y() + shift * versorOnLOR.Y(),
                             middleOfLOR.Z() + shift * versorOnLOR.Z());
  return annihilationPoint;
}

void saveReportToFile(const std::string& outFileName, ROOT::RDF::RCutFlowReport& report)
{
  std::ofstream file(outFileName.c_str());
  file << std::left << std::setw(30) << "Cut condition" << std::setw(30) << "All" << std::setw(20) << "Passed" << std::setw(15) << "Percentege" << std::endl;
  for (auto  cutInfo : report) {
    file << std::left << std::setw(30) << cutInfo.GetName() << std::setw(30) << cutInfo.GetAll() << std::setw(20) << cutInfo.GetPass() << std::setw(15)
         << cutInfo.GetEff() << " %" << std::endl;
  }
}
