// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Projections/GeneratedPercentileProjection.hh"
#include "Rivet/Projections/UserCentEstimate.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#define _USE_MATH_DEFINES

namespace Rivet {
  /// @brief Centrality projection for PHENIX AuAu.
class ALICE_2021_AverageEventParticles : public Analysis {

  public:

    //RHIC_2019_CentralityCalibration() : Analysis("ALICE_2021_AverageEventParticles") { };
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2021_AverageEventParticles);
  /// Book histograms and initialise projections before the run
  void init() {
    // One projection for the actual observable, and one for the
    // generated impact parameter.

    //string cuts = getOption<string>("cuts","");
    //MSG_INFO("RHIC Experiment: " << experiment);


    const ALICE::PrimaryParticles app(Cuts::abseta < 0.8 && Cuts::pT > 0.150*GeV && Cuts::abscharge > 0);
    declare(app, "app");

    // The calibration histogram:
    book(_calib, "AVERAGE_N", 1, 0., 1.);


  }

  /// Perform the per-event analysis
  void analyze(const Event& event) {

    // The alternative centrality based on generated impact
    // parameter, assumes that the generator does not describe the
    // full final state, and should therefore be filled even if the
    // event is not triggered.

      _calib->fill(0.5, apply<ALICE::PrimaryParticles>(event, "app").particles().size());

      Particles particles = apply<ALICE::PrimaryParticles>(event, "app").particlesByPt();
      cout << "New event" << endl;
      for(auto p : particles)
      {
          cout << p.pT()/GeV << endl;
      }

  }

  /// Finalize
  void finalize() {

  }

  /// The calibration histograms.
  Profile1DPtr _calib;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2021_AverageEventParticles);
}
