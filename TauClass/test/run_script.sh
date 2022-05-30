#!/bin/sh

mkdir files
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cd ${CMSSW_BASE}/src
cmsenv
cd -

for i in i in {1..299..$2}
    ./NohupParser.sh $1 ${i} $(expr ${i} + $2)
done

hadd -f final.root files/histos_*.root
run_plot final.root