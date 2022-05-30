#!/bin/bash

mkdir files

for filename in /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/lowpt/*.root do
    ./run_gen $filename ${filename//[^0-9]/}
done

hadd -f final.root files/histos_*.root
./run_plot final.root
    