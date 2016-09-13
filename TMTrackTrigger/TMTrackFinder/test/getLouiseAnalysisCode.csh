#!/bin/tcsh
######################################################################################
# This gets the official (CMS L1 Track group approved) code from Louise Skinnari for 
# analysis L1 track objects (TTTrack objects).
# https://github.com/skinnari/cmssw/tree/TTI_62X_TrackTriggerObjects/SLHCUpgradeSimulations/L1TrackTrigger/test/
#
# Run this periodically to get her latest code and compare it with that in our
# SVN repository, updating the latter if necessary.
#
# Her code is documented in 
# https://twiki.cern.ch/twiki/bin/view/CMS/TrackerUpgradeL1TriggerPerformance
######################################################################################

# Name of directory where Louise's analysis code is.
setenv louiseGitDir "https://raw.githubusercontent.com/skinnari/cmssw/L1TK_62X_SLHC28/SLHCUpgradeSimulations/L1TrackTrigger/test/"

echo 'Code being copied from git to files named L1*.DontRun' 
echo 'You should compare these with corresponding files in svn and check that you understand the differences'
echo 'The _cfg.py file will be significantly different, as it must run our TMT tracking instead of tracklet tracking'

# curl -k ${louiseGitDir}/L1TrackNtupleMaker_cfg.py >! L1TrackNtupleMaker_cfg.py.DontRun 
curl -k ${louiseGitDir}/L1TrackNtupleMaker.cc     >! L1TrackNtupleMaker.cc.DontRun
curl -k ${louiseGitDir}/L1TrackNtuplePlot.C       >! L1TrackNtuplePlot.C.DontRun
# curl -k ${louiseGitDir}/BuildFile.xml             >! BuildFile.xml.DontRun

# Since Louise keeps her code in SLHCUpgradeSimulations/L1TrackTrigger/test/, whereas we copy it to
# TMTrackTrigger/TMTrackFinder/test/ , change any references to the former directory name.

sed -i "s|SLHCUpgradeSimulationsL1TrackTriggerTests|TMTrackTriggerTMTrackFinderTests|g" BuildFile.xml.DontRun
