#ifndef analysis_test_analysis_h
#define analysis_test_analysis_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class test_analysis : public AnalyzerBase
{
  INIT_ANALYSIS(test_analysis,"test_analysis")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
  std::vector<RecLeptonFormat> GetBaselineElectrons(const std::vector<RecLeptonFormat>& candidate_electrons);
  std::vector<RecLeptonFormat> GetBaselineMuons(const std::vector<RecLeptonFormat>& candidate_muons);
  std::vector<RecJetFormat> GetBaselineJets(const std::vector<RecJetFormat>& candidate_jets);
  std::vector<RecLeptonFormat> GetSignalElectrons(const std::vector<RecLeptonFormat>& baseline_electrons);
  std::vector<RecLeptonFormat> GetSignalMuons(const std::vector<RecLeptonFormat>& baseline_muons);
  std::vector<RecJetFormat> GetSignalJets(const std::vector<RecJetFormat>& baseline_jets);

  /**
   * @brief Removes overlaping jets and leptons
   * @param jets
   * @param electrons
   * @param muons
   */
  void RemoveOverlap(std::vector<RecJetFormat>& jets, std::vector<RecLeptonFormat>& electrons,
                     std::vector<RecLeptonFormat>& muons);
};
}

#endif