 #!/bin/bash
ANALYSIS_DIR=$PWD
PYTHIA_EX_DIR=/nics/a/proj/UTK0019/Pythia8/pythia8243/examples
#cd $ANALYSIS_DIR/Centrality
#rm *.so
#export RIVET_ANALYSIS_PATH=$PWD
#rivet-buildplugin RivetRHIC_2019_CentralityCalibration.so RHIC_2019_CentralityCalibration.cc
#rivet --pwd -a RHIC_2019_CentralityCalibration:exp=STAR -o calibration_AuAu_200GeV_STAR.yoda --ignore-beams $PYTHIA_EX_DIR/hepmc_AuAu_MB_1.hepmc $PYTHIA_EX_DIR/hepmc_AuAu_MB_2.hepmc $PYTHIA_EX_DIR/hepmc_AuAu_MB_3.hepmc $PYTHIA_EX_DIR
#cd ..
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetSTAR_2003_I619063.so STAR_2003_I619063.cc
rivet --pwd -p $ANALYSIS_DIR/Centrality/calibration_AuAu_200GeV_STAR.yoda -a STAR_2003_I619063:cent=GEN -o Angantyr.yoda $PYTHIA_EX_DIR/hepmc_AuAu_200GeV_500Events.hepmc
#rivet-mkhtml --pwd Rivet.yoda
