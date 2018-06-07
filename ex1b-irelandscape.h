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

    void CheckOverlap (std::vector<RecLeptonFormat>& baseline_electrons,
                       std::vector<RecLeptonFormat>& baseline_muons,
                       std::vector<RecJetFormat>& baseline_jets);

    static inline bool jl_overlap_criteria_met (const RecJetFormat& j,
                                                const RecLeptonFormat& l)
    {
        double dr = j.dr(l);

        return dr < std::min(0.4, 0.04 + (10.0 / l.pt()));
    }
                                        
};
}

#endif
