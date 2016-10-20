#!/usr/bin/env bash
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters, needs 2 argument 1) name of mip server number, 2) num of threads to use"
    exit
fi

mipster mipIllumExtractByArmAndFilterMultiple       --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile extractRun_run1 --allowableErrors 6
mipster mipBarcodeCorrectionMultiple                --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipBarcodeCorrecting_run1 --allowableErrors 6  
mipster mipClusteringMultiple                       --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipClustering_run1
mipster mipPopulationClusteringMultiple             --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipPopClustering_run1 
nohup mipster mav			                        --masterDir $(realpath ./)  --numThreads $2 --port 1000$1 --name mip$1 --verbose &

