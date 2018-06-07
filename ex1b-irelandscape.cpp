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

void test_analysis::CheckOverlap (std::vector<RecLeptonFormat>& electrons,
                                  std::vector<RecLeptonFormat>& muons,
                                  std::vector<RecJetFormat>& jets)
{
  // (e,j) check
  for (std::vector<RecLeptonFormat>::const_iterator it_electron = electrons.begin();
       it_electron != electrons.end();
       )
  {
      bool electron_removed = false;

      for (std::vector<RecJetFormat>::const_iterator it_jet = jets.begin();
           it_jet != jets.end();
           )
      {
          bool jet_removed = false;
          double dr = it_electron->dr(*it_jet);

          // DELTA R < 0.2
          if (dr < 0.2)
          {
              cout << "(e,j) overlap detected: DR = " << dr << endl;
              if (it_jet->btag())
              {
                  cout << "Jet b-tagged, keeping jet" << endl;
                  it_electron = electrons.erase(it_electron);
                  electron_removed = true;
                  break;
              }
              else
              {
                  cout << "Jet is not b-tagged, keeping electron" << endl;
                  it_jet = jets.erase(it_jet);
                  jet_removed = true;
              }
          }

          if (!jet_removed)
              ++it_jet;
      }

      if (!electron_removed)
          ++it_electron;
  }
  
  // (j,l) check
  for (std::vector<RecJetFormat>::const_iterator it_jet = jets.begin();
       it_jet != jets.end();
       ++it_jet)
  {
      // First, iterate electrons
      for (std::vector<RecLeptonFormat>::const_iterator it_electron = electrons.begin();
           it_electron != electrons.end();
           )
      {
          if (test_analysis::jl_overlap_criteria_met(*it_jet, *it_electron))
          {
              cout << "(j,e) overlap detected. Keeping jet" << endl;
              it_electron = electrons.erase(it_electron);
          }
          else
              ++it_electron;
      }
      
      // Next, iterate muons
      for (std::vector<RecLeptonFormat>::const_iterator it_muon = muons.begin();
           it_muon != muons.end();
           )
      {
          if (test_analysis::jl_overlap_criteria_met(*it_jet, *it_muon))
          {
              cout << "(j,u) overlap detected. Keeping jet" << endl;
              it_muon = muons.erase(it_muon);
          }
          else
              ++it_muon;
      }
  }

}


// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool test_analysis::Execute(SampleFormat& sample, const EventFormat& event)
{
  if (event.rec() == NULL)
      return true;

  cout << "Initial Electrons: " << event.rec()->electrons().size() << ", "
       << "Initial Muons: " << event.rec()->muons().size() << ", "
       << "Initial Jets: " << event.rec()->jets().size() << endl;

  // Extract baseline electrons (pt > 5GeV and |n| < 2.47
  std::vector<RecLeptonFormat> electrons;
  for (std::vector<RecLeptonFormat>::const_iterator it_electron = event.rec()->electrons().begin();
       it_electron != event.rec()->electrons().end();
       ++it_electron)
  {
      MAfloat32 pt = it_electron->pt();
      MAfloat32 abseta = it_electron->abseta();

      if (pt > 5 && abseta < 2.47)
      {
          cout << "Adding baseline electron (pt = " << pt << ", |n| = " << abseta << ")" << endl;
          electrons.push_back(*it_electron);
      }
      else
      {
          cout << "Rejecting electron (pt = " << pt << ", |n| = " << abseta << ")" << endl;
      }
  }

  cout << "Baseline electrons: " << electrons.size() << endl;
  
  // Extract baseline muons (pt > 4GeV and |n| <= 2.7
  std::vector<RecLeptonFormat> muons;
  for (std::vector<RecLeptonFormat>::const_iterator it_muon = event.rec()->muons().begin();
       it_muon != event.rec()->muons().end();
       ++it_muon)
  {
      MAfloat32 pt = it_muon->pt();
      MAfloat32 abseta = it_muon->abseta();

      if (pt > 4 && abseta <= 2.47)
      {
          cout << "Adding baseline muon (pt = " << pt << ", |n| = " << abseta << ")" << endl;
          muons.push_back(*it_muon);
      }
      else
      {
          cout << "Rejecting muon (pt = " << pt << ", |n| = " << abseta << ")" << endl;
      }
  }

  cout << "Baseline muons: " << muons.size() << endl;
  
  // Extract baseline jets (pt > 20GeV)
  std::vector<RecJetFormat> jets;
  for (std::vector<RecJetFormat>::const_iterator it_jet = event.rec()->jets().begin();
       it_jet != event.rec()->jets().end();
       ++it_jet)
  {
      MAfloat32 pt = it_jet->pt();

      if (pt > 20)
      {
          cout << "Adding baseline jet (pt = " << pt << ")" << endl;
          jets.push_back(*it_jet);
      }
      else
      {
          cout << "Rejecting jet (pt = " << pt << ")" << endl;
      }
  }

  cout << "Baseline jets: " << jets.size() << endl;

  CheckOverlap(electrons,
               muons,
               jets);

  cout << "***" << endl;
  cout << "Baseline electrons (after overlap removal): " << electrons.size() << endl;
  cout << "Baseline muons (after overlap removal): " << muons.size() << endl;
  cout << "Baseline jets (after overlap removal): " << jets.size() << endl;
  cout << "***" << endl;
  
  // Extract signal electrons (pt >= 25GeV)
  for (std::vector<RecLeptonFormat>::const_iterator it_electron = electrons.begin();
       it_electron != electrons.end();
       )
  {
      MAfloat32 pt = it_electron->pt();

      if (pt < 25)
      {
          cout << "Rejecting non-signal electron with pt " << pt << endl;
          it_electron = electrons.erase(it_electron);
      }
      else
      {
          ++it_electron;
      }
  }
  
  // Extract signal muons (pt >= 25GeV)
  for (std::vector<RecLeptonFormat>::const_iterator it_muon = muons.begin();
       it_muon != muons.end();
       )
  {
      MAfloat32 pt = it_muon->pt();

      if (pt < 25)
      {
          cout << "Rejecting non-signal muon with pt " << pt << endl;
          it_muon = muons.erase(it_muon);
      }
      else
      {
          ++it_muon;
      }
  }
  
  // Extract signal jets (pt > 25GeV and |n| < 2.5)
  for (std::vector<RecJetFormat>::const_iterator it_jet = jets.begin();
       it_jet != jets.end();
       )
  {
      MAfloat32 pt = it_jet->pt();
      MAfloat32 abseta = it_jet->abseta();

      if (pt <= 25 || abseta >= 2.5)
      {
          cout << "Rejecting non-signal jet with pt=" << pt << ", |n|=" << abseta << endl;
          it_jet = jets.erase(it_jet);
      }
      else
      {
          ++it_jet;
      }
  }

  cout << "***" << endl;
  cout << "Signal electrons: " << electrons.size() << endl;
  cout << "Signal muons: " << muons.size() << endl;
  cout << "Signal jets: " << jets.size() << endl;
  cout << "***" << endl;

  return true;
}

