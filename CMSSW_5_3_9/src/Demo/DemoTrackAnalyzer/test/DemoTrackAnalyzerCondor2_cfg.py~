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

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:/data/users/jbrinson/FiveHits_3/WtoLNuFiveHits_Five_Reco' + processPostfix + '.root')
'file:/data/users/jbrinson/FiveHits_2/WtoLNuFiveHits_Reco_' + processPostfix + '.root')

)

process.demo = cms.EDAnalyzer("DemoTrackAnalyzer",
                              tracks       = cms.untracked.InputTag("generalTracks"),
                              muons        = cms.untracked.InputTag("globalMuons"),
                              smuons        = cms.untracked.InputTag("standAloneMuons"),
                              genParts     = cms.untracked.InputTag("genParticles"),
                              simvtxTag    = cms.untracked.InputTag("g4SimHits"),
                              isPiPartGun  = cms.untracked.bool(False),
                              
    )

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string('analyzerOutput_2/hist_' + processPostfix + '.root')
                                   fileName = cms.string('analyzerOutput/hist_' + processPostfix + '.root')
                                   )

process.p = cms.Path(process.demo)
