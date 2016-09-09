#!/bin/tcsh
################################################################################################
# This gets the latest L1 track fitting code from the tracklet group.
# https://github.com/EmanuelPerez/cmssw/tree/TTI_62X_TrackTriggerObjects/SimDataFormats/SLHC/
#
# Run this periodically to get the latest code and compare it with that in our
# SVN repository, updating the latter if necessary.
################################################################################################

# Name of directory where Louise's analysis code is.
setenv fitGitDir "https://raw.githubusercontent.com/EmanuelPerez/cmssw/TTI_62X_TrackTriggerObjects/SimDataFormats/SLHC/interface/"

# Copy code from the git respository to the local linux directory here.
echo 'Code being copied from git to files named L1*.DontRun in ../interface/' 
echo 'You should compare these with corresponding files in svn and check that you understand the differences'
curl -k ${fitGitDir}/L1TTrack.hh  >! ../interface/L1TTrack.hh.DontRun
