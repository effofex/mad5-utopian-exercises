#include "SampleAnalyzer/User/Analyzer/test_analysis.h"
using namespace MA5;
using namespace std;


std::vector<RecLeptonFormat> test_analysis::GetBaselineElectrons(const std::vector<RecLeptonFormat> &candidate_electrons) {
  std::vector<RecLeptonFormat> result;
  for (auto& candidate : candidate_electrons) {
    if(candidate.pt() <= 5.7)
      continue;
    if(candidate.abseta() >= 2.47)
      continue;
    result.push_back(candidate);
  }
  return result;
}

std::vector<RecLeptonFormat> test_analysis::GetBaselineMuons(const std::vector<RecLeptonFormat> &candidate_muons) {
  std::vector<RecLeptonFormat> result;
  for (auto& candidate : candidate_muons) {
    if(candidate.pt() <= 4.0)
      continue;
    if(candidate.abseta() >= 2.7)
      continue;
    result.push_back(candidate);
  }
  return result;
}

std::vector<RecJetFormat> test_analysis::GetBaselineJets(const std::vector<RecJetFormat> &candidate_jets) {
  std::vector<RecJetFormat> result;
  for (auto& candidate : candidate_jets) {
    if(candidate.pt() <= 20.0)
      continue;
    result.push_back(candidate);
  }
  return result;
}

void test_analysis::RemoveOverlap(std::vector<RecJetFormat> &jets,
                                  std::vector<RecLeptonFormat> &electrons,
                                  std::vector<RecLeptonFormat> &muons) {
  // Second step from Table 3.
  for(auto& electron: electrons) {
    for(auto jet_it = jets.begin(); jet_it != jets.end(); jet_it++) {
      if(jet_it->dr(electron) < 0.2) {
        jets.erase(jet_it);
        // we have to move iterator back after removing the element
        jet_it--;
      }
    }
  }

  for(auto& jet: jets) {
    // first compare all electrons with the jet
    for(auto electron_it = electrons.begin(); electron_it != electrons.end(); electron_it++) {
      MAfloat32 threshold = min(0.4, 0.04 + 10.0 / electron_it->pt());
      if(electron_it->dr(jet) < threshold) {
        electrons.erase(electron_it);
        electron_it--;
      }
    }
    // and then all muons
    for(auto muon_it = muons.begin(); muon_it != muons.end(); muon_it++) {
      MAfloat32 threshold = min(0.4, 0.04 + 10.0 / muon_it->pt());
      if(muon_it->dr(jet) < threshold) {
        muons.erase(muon_it);
        muon_it--;
      }
    }
  }
}

std::vector<RecLeptonFormat> test_analysis::GetSignalElectrons(const std::vector<RecLeptonFormat> &baseline_electrons) {
  std::vector<RecLeptonFormat> result;
  for (auto& electron : baseline_electrons) {
    if(electron.pt() <= 25.0)
      continue;
    result.push_back(electron);
  }
  return result;
}

std::vector<RecLeptonFormat> test_analysis::GetSignalMuons(const std::vector<RecLeptonFormat> &baseline_muons) {
  std::vector<RecLeptonFormat> result;
  for (auto& muon : baseline_muons) {
    if(muon.pt() <= 25.0)
      continue;
    result.push_back(muon);
  }
  return result;
}

std::vector<RecJetFormat> test_analysis::GetSignalJets(const std::vector<RecJetFormat> &baseline_jets) {
  std::vector<RecJetFormat> result;
  for (auto& jet : baseline_jets) {
    if(jet.pt() <= 25.0)
      continue;
    if(jet.abseta() >= 2.5)
      continue;
    result.push_back(jet);
  }
  return result;
}
// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool test_analysis::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void test_analysis::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool test_analysis::Execute(SampleFormat& sample, const EventFormat& event)
{
  cout << "Initial electrons count: " << event.rec()->electrons().size()
       << ", initial muons count: " << event.rec()->muons().size()
       << ", initial jets count: " << event.rec()->jets().size() << endl;

  auto baseline_electrons = GetBaselineElectrons(event.rec()->electrons());
  auto baseline_muons = GetBaselineMuons(event.rec()->muons());
  auto baseline_jets = GetBaselineJets(event.rec()->jets());

  cout << "Baseline electrons count: " << baseline_electrons.size()
       << ", baseline muons count: " << baseline_muons.size()
       << ", baseline jets count: " << baseline_jets.size() << endl;

  RemoveOverlap(baseline_jets, baseline_electrons, baseline_muons);

  cout << "After overlap removal:" << endl;
  cout << "Baseline electrons count: " << baseline_electrons.size()
       << ", baseline muons count: " << baseline_muons.size()
       << ", baseline jets count: " << baseline_jets.size() << endl;

  auto signal_electrons = GetSignalElectrons(baseline_electrons);
  auto signal_muons = GetSignalMuons(baseline_muons);
  auto signal_jets = GetSignalJets(baseline_jets);

  cout << "signal electrons count: " << signal_electrons.size()
       << ", signal muons count: " << signal_muons.size()
       << ", signal jets count: " << signal_jets.size() << endl;


  
  return true;
}

