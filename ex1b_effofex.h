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
  virtual void filterBaseLine(const std::vector<RecJetFormat> &to, std::vector<RecJetFormat> &from,
                              double min_pt, double maxAbsEta);
  virtual void filterBaseLine(const std::vector<RecLeptonFormat> &to, std::vector<RecLeptonFormat> &from,
                              double min_pt, double maxAbsEta);
};
}

#endif
