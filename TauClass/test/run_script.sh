#!/bin/sh

mkdir files
mkdir jobs
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cd ${CMSSW_BASE}/src
cmsenv
cd -

start=1
end=$(expr $2 - 1)
for counter in {1..299}
do
    echo $counter. "\$start = " $start " and \$end = " $end
    ./NohupParser.sh $1 ${start} ${end}
    if [[ $start -eq 1 ]];
    then
        start=0
    fi
    start=$(($start+$2))
    end=$(($end+$2))
    if [[ $end -ge 299 ]];
    then 
        end=299
    fi
    if [[ $start -gt 299 ]];
    then
        break
    fi
done

