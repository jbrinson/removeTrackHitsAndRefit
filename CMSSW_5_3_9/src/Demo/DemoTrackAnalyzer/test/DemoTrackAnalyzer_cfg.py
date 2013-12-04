import FWCore.ParameterSet.Config as cms

process = cms.Process("DemoTwo")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

import sys
import os
import glob
import fnmatch

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
	#fileNames = cms.untracked.vstring(dir+'WToLNu_TuneZ2star_8TeV_pythia6_tauola_RECO_'+processPostfix+'_'+ " ".join(glob.glob('*.root')))
#                            fileNames = cms.untracked.vstring(testString.join(glob.glob('*.root')))
#                            fileNames = cms.untracked.vstring()
#    'file:/afs/cern.ch/work/j/jbrinson/public/disappTrks/AnTest/removeTrackHits/CMSSW_5_3_9/src/pionPartGunFiveHits_Reco.root')
#    'file:/afs/cern.ch/work/j/jbrinson/public/disappTrks/AnTest/removeTrackHits/CMSSW_5_3_9/src/singlePionPartGunFiveHits_Reco_10K.root')
#'file:/afs/cern.ch/work/j/jbrinson/public/disappTrks/AnTest/studyMisOutHits/CMSSW_5_3_3/src/pionPartGun/batch/output/SinglePiPt20_cfi_py_Ideal_RECO_7.root ')
#'file:/afs/cern.ch/work/j/jbrinson/public/disappTrks/AnTest/studyMisOutHits/CMSSW_5_3_3/src/pionPartGun/batch/output/reallySinglePiPt20_cfi_py_Ideal_RECO_8.root ')
#'file:/afs/cern.ch/work/j/jbrinson/public/disappTrks/AnTest/removeTrackHits/CMSSW_5_3_9/src/WtoLNuFiveHits_1K_Reco.root')
#'file:/home/jbrinson/DisTrack/mainAN/AnTemp/removeHits/CMSSW_5_3_9/src/WtoLNuFiveHits_Reco_' + processPostfix + '.root')
#'file:singlePionPartGunFiveHits_Reco_10K.root')
#'file:SinglePiFiveHits_Pt50_Reco.root')
#'file:SinglePiFiveHits_Pt40_Reco.root')
#'file:SinglePiFiveHits_Pt30_Reco.root')
#'file:SinglePiFiveHits_Pt10_Reco.root')
#'file:SinglePiSixHits_Pt10_Reco.root')
#'file:SinglePiSixHits_Pt30_Reco.root')
#'file:SinglePiSixHits_Pt50_Reco.root')
#'file:SinglePiSixHits_Pt10_Reco.root')
'file:SinglePiSixHits_Pt40_Reco.root')
#    'root://xrootd.unl.edu//store/user/jbrinson/WToTauNu_PU_2/WToLNu_TuneZ2star_8TeV_pythia6_tauola_RECO_9_1_TLw.root')
#    'root://xrootd.unl.edu//store/user/jbrinson/WToTauNu_PU_2/WToLNu_TuneZ2star_8TeV_pythia6_tauola_RECO_9' + glob.glob('*.root') )
#    )
)
#dir = '/data/users/jbrinson/FiveHits/'
#dir = '/mnt/hadoop/se/store/user/jbrinson/WToTauNu_PU/'

#for file in os.listdir(dir):
#    print file
    #if file.find(".root") != -1:  # Skip over files that do not contain .root.
 #   print file
 #if file.find('WToLNu_TuneZ2star_8TeV_pythia6_tauola_RECO_'+processPostFix+'_') != -1:  # Skip over files that do not contain .root.
    
        #                process.source.fileNames.extend(cms.untracked.vstring('file:' + dir + file))  
 #               process.source.fileNames.extend(cms.untracked.vstring(dir + file))  


process.demo = cms.EDAnalyzer("DemoTrackAnalyzer",
                              tracks       = cms.untracked.InputTag("generalTracks"),
                              taus         = cms.untracked.InputTag("hpsPFTauProducer"),
                              muons        = cms.untracked.InputTag("globalMuons"),
                              smuons        = cms.untracked.InputTag("standAloneMuons"),
                              genParts     = cms.untracked.InputTag("genParticles"),
                              simvtxTag    = cms.untracked.InputTag("g4SimHits"),
                              isPiPartGun  = cms.untracked.bool(True),
    )

process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string('histPiPt50.root')
                                   #fileName = cms.string('histPiPt40.root')
                                   #fileName = cms.string('histPiPt30.root')
#                                   fileName = cms.string('histPiPt10.root')
#                                   fileName = cms.string('histPiPt30SixHit.root')
#          fileName = cms.string('histPiPt50SixHit.root')
#          fileName = cms.string('histPiPt10SixHit.root')
         fileName = cms.string('histPiPt40SixHit.root')
                                   #fileName = cms.string('histPiPt20.root')
                                   )

process.p = cms.Path(process.demo)
