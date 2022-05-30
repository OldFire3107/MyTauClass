#!/bin/sh 
# $1 - input files path
# $2 - first index
# $3 - last index

cat > jobs/runfrom${2}to${3}.sh <<EOF1
for i in {$2..$3}; do
    echo \${i}
    run_gen $1/flatTuple_\${i}.root \${i}
done
EOF1

chmod u+x jobs/runfrom${2}to${3}.sh
nohup ./jobs/runfrom${2}to${3}.sh > jobs/nohup_runfrom${2}to${3}.out &
