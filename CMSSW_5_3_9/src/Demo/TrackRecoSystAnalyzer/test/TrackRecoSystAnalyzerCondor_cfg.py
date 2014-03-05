import FWCore.ParameterSet.Config as cms

process = cms.Process("DemoTwo")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

import sys
import os
import glob
import fnmatch
processNumber = int (sys.argv[2])
processPostfix = str (processNumber)  

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
##For five hits
#'file:/data/users/jbrinson/FiveHits/WtoLNuFiveHits_Reco_' + processPostfix + '.root') # this directory contains electrons and muons

##for all hits
'file:/mnt/hadoop/se/store/user/jbrinson/WToTauNu_PU/WToLNu_TuneZ2star_8TeV_pythia6_tauola_RECO_' + processPostfix + '.root') # this directory contains electrons and muons

)

process.demo2 = cms.EDAnalyzer("TrackRecoSystAnalyzer",
                              muons        = cms.untracked.InputTag("muons"),
                              genParts        = cms.untracked.InputTag("genParticles"),
                              
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('analyzerOutputAllHits/hist_' + processPostfix + '.root')
                                   
                                   )

process.p = cms.Path(process.demo2)
