// -*- C++ -*-
#ifndef RIVET_AverageN_HH
#define RIVET_AverageN_HH

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"


namespace Rivet
{

    /// @brief Project Average Event Particles
    class AverageN : public Projection
    {
        public:

            AverageN(const ALICE::PrimaryParticles app, YODA::Profile1DPtr profileN);
            AverageN(const FinalState fs, YODA::Profile1DPtr profileN);

            /// Clone on the heap.
            DEFAULT_RIVET_PROJ_CLONE(AverageN);

            double averageN() const { return _averageN; }

            double getLeadDeltaPhi(Particle& particle);

            Particle getLeadParticle() const { return _leadingParticle; }

            double operator()() const {
              return averageN();
            }




        protected:

            /// Apply the projection to the event.
            void project(const Event& e);

            CmpState compare(const Projection& p) const;

            string _jetProjName;
            string _particlesProjName;
            double _averageN;
            Particle _leadingParticle;

    };
}

#endif
