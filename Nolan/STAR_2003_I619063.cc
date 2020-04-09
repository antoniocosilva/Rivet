// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Centrality/RHICCentrality.hh" //external header for Centrality calculation
namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2003_I619063 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2003_I619063);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      //Centrality
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");  
        
      // the basic final-state projection: 
      // all final-state particles within 
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9 && Cuts::abscharge > 0);
      declare(fs,"fs");
      

      // Book histograms
      // specify custom binning
      book(chSpectrum["chSpectrum0_5"], 1, 1, 1);
      book(chSpectrum["chSpectrum5_10"], 1, 1, 2);
      book(chSpectrum["chSpectrum10_20"], 1, 1, 3);
      book(chSpectrum["chSpectrum20_30"], 1, 1, 4);
      book(chSpectrum["chSpectrum30_40"], 1, 1, 5);
      book(chSpectrum["chSpectrum40_60"], 1, 1, 6);
      book(sow["sow0_5"], "sow0_5");
      book(sow["sow5_10"], "sow5_10");
      book(sow["sow0_5"], "sow10_20");
      book(sow["sow5_10"], "sow20_30");
      book(sow["sow0_5"], "sow30_4");
      book(sow["sow5_10"], "sow40_60");
      book(sow["sow5_10"], "sow60_80");
      book(Rcp["Rcp0_5_over_5_10"], "Rcp0_5_over_5_10");
      
      string refnameCentRatio05_4060 = mkAxisCode(4,1,1);
      book(Rcp["Rcp0_5_over_40_60"], refnameCentRatio05_4060);
      
      string refnameCentRatio05_6080 = mkAxisCode(4,1,2);
      const Scatter2D& refdataCentRatio05_6080 = refData(refnameCentRatio05_6080);
      book(chSpectrum["chSpectrum60_80"], refnameCentRatio05_6080 + "_6080", refdataCentRatio05_6080);
      book(Rcp["Rcp0_5_over_60_80"], refnameCentRatio05_6080);
      
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      //book(_h["AAAA"], 1, 1, 1);
      //book(_p["BBBB"], 2, 1, 1);
      //book(_c["CCCC"], 3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      
      //Get cenrality  
      const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
      const double c = cent();  

      FinalState fs = applyProjection<FinalState>(event,"fs");
      
      Particles particles = fs.particles();
      
      for(Particle p : particles)
      {
          if(c < 5.)
          {
              chSpectrum["chSpectrum0_5"]->fill(p.pT()/GeV);
              sow["sow0_5"]->fill();
          }
          else if(c >= 5 && c < 10)
          {
              chSpectrum["chSpectrum5_10"]->fill(p.pT()/GeV);
              sow["sow5_10"]->fill();
          }
          else if(c >= 10 && c < 20)
          {
              chSpectrum["chSpectrum10_20"]->fill(p.pT()/GeV);
              sow["sow10_20"]->fill();
          }
          else if(c >= 20 && c < 30)
          {
              chSpectrum["chSpectrum20_30"]->fill(p.pT()/GeV);
              sow["sow20_30"]->fill();
          }
          else if(c >= 30 && c < 40)
          {
              chSpectrum["chSpectrum30_40"]->fill(p.pT()/GeV);
              sow["sow30_40"]->fill();
          }
          else if(c >= 40 && c < 60)
          {
              chSpectrum["chSpectrum40_60"]->fill(p.pT()/GeV);
              sow["sow40_60"]->fill();
          }
          else if(c >= 60 && c < 80)
          {
              chSpectrum["chSpectrum60_80"]->fill(p.pT()/GeV);
              sow["sow60_80"]->fill();
          }
      }
        

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      chSpectrum["chSpectrum0_5"]->scaleW(1./sow["sow0_5"]->sumW());
      chSpectrum["chSpectrum5_10"]->scaleW(1./sow["sow5_10"]->sumW());
      chSpectrum["chSpectrum10_20"]->scaleW(1./sow["sow10_20"]->sumW());
      chSpectrum["chSpectrum20_30"]->scaleW(1./sow["sow20_30"]->sumW());
      chSpectrum["chSpectrum30_40"]->scaleW(1./sow["sow30_40"]->sumW());
      chSpectrum["chSpectrum40_60"]->scaleW(1./sow["sow40_60"]->sumW());
      chSpectrum["chSpectrum60_80"]->scaleW(1./sow["sow60_80"]->sumW());
        
      //normalize(_h["YYYY"]); // normalize to unity
      //scale(_h["ZZZZ"], crossSection()/picobarn/sumOfWeights()); // norm to cross section
      
      //
      divide(chSpectrum["chSpectrum0_5"], chSpectrum["chSpectrum40_60"], Rcp["Rcp0_5_over_40_60"]);
      divide(chSpectrum["chSpectrum0_5"], chSpectrum["chSpectrum60_80"], Rcp["Rcp0_5_over_60_80"]);

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> chSpectrum;
    map<string, CounterPtr> sow;
    map<string, Scatter2DPtr> Rcp;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_2003_I619063);


}
