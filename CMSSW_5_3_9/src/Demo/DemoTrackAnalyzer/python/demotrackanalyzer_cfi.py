import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('DemoTrackAnalyzer'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
