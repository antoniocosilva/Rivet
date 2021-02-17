// -*- C++ -*-
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/AverageN.hh"

namespace Rivet {

    //Constructors
    AverageN::AverageN(const ALICE::PrimaryParticles app, YODA::Profile1DPtr profileN) :
    _particlesProjName("particlesALICEPrimary")
    {
        setName("AverageN");
        declare(app,_particlesProjName);

        _averageN = profileN->bin(0).mean();
    }

    AverageN::AverageN(const FinalState fs, YODA::Profile1DPtr profileN) :
    _particlesProjName("particlesFinalState")
    {
        setName("AverageN");
        declare(fs,_particlesProjName);

        _averageN = profileN->bin(0).mean();

    }

    double AverageN::getLeadDeltaPhi(Particle& particle)
    {
        double deltaPhi = abs(_leadingParticle.phi() - particle.phi());

        if(deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

        return deltaPhi;
    }

    void AverageN::project(const Event& e)
    {

        Particles particles;

        if(_particlesProjName.compare("particlesALICEPrimary") != 0) particles = apply<ALICE::PrimaryParticles>(e, _particlesProjName).particlesByPt();
        else if(_particlesProjName.compare("FinalState") != 0) particles = apply<FinalState>(e, _particlesProjName).particlesByPt();

        _leadingParticle = particles[0];

    }

    CmpState AverageN::compare(const Projection& p) const
    {
        const AverageN& other = dynamic_cast<const AverageN&>(p);
        if(_particlesProjName.compare(other._particlesProjName) != 0) return CmpState::NEQ;

        return CmpState::EQ;
    }



}
