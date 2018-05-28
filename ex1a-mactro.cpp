#include "SampleAnalyzer/User/Analyzer/test_analysis.h"
using namespace MA5;
using namespace std;

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
  cout << "Electrons count: " << event.rec()->electrons().size() << endl;
  cout << "Photons count: " << event.rec()->photons().size() << endl;
  cout << "Muons count: " << event.rec()->muons().size() << endl;
  
  return true;
}

