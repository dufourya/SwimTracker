addpath(genpath('/mnt/research/dufourlab/GitLab/SwimTracker'))   
mcc -R -nodisplay -R -singleCompThread -m runTrackSingle.m -a /mnt/research/dufourlab/GitLab/SwimTracker/lib/u-track/kalman*.m -a /mnt/research/dufourlab/GitLab/SwimTracker/lib/u-track/costMat*.m
exit
