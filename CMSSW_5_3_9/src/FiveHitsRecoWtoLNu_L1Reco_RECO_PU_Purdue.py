# Auto generated configuration file
# using: 
# Revision: 1.381.2.18 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: FiveHitsRecoWtoLNu --conditions START53_V27::All -s L1Reco,RECO --eventcontent RECOSIM --datatier GEN-SIM-RECO -n 1000 --fileout pionPartGunFiveHits_Reco.root --filein file:/afs/cern.ch/work/j/jbrinson/public/disappTrks/AnTest/studyMisOutHits/CMSSW_5_3_3/src/pionPartGun/batch/output/SinglePiPt20_cfi_py_Ideal_RECO_7.root --pileup 2012_Summer_50ns_PoissonOOTPU --no_exec
import FWCore.ParameterSet.Config as cms
import os

import sys
processNumber = int (sys.argv[2])
processPostfix = str (processNumber)

process = cms.Process('RECO')
#print os.environ['HOME']
os.system('printenv')
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/user/jbrinson/WToTauNu_PU_2/WToLNu_TuneZ2star_8TeV_pythia6_tauola_RECO_'+processPostfix+'.root')
)
process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.18 $'),
    annotation = cms.untracked.string('FiveHitsRecoWtoLNu nevts:1000'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/data/users/jbrinson/FiveHits_3/WtoLNuFiveHits_Five_Reco'+processPostfix+'.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V27::All', '')

# Path and EndPath definitions
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)

