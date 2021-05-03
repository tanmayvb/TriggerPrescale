import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


# source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00B6F8B6-90F1-E311-B72C-0025905A6092.root'
'/store/data/Run2017B/JetHT/MINIAOD/UL2017_MiniAODv2-v1/260000/00EA7817-B2DB-E44A-9FD7-3FFB3778C024.root'
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/50000/00280B8B-36E0-E711-A549-02163E01421E.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FEA2ED14-5CDF-E711-ACA6-02163E012AF0.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FAD0B094-CBDF-E711-8E6A-A0369FC5FBA4.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/50000/8261436C-85DF-E711-8BCA-02163E01A1EB.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/50000/8257B6C4-04E0-E711-82B7-02163E0137A0.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/50000/82335A18-C8DF-E711-8880-A4BF0112BD40.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FEA2ED14-5CDF-E711-ACA6-02163E012AF0.root'
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FEE33912-EFDF-E711-91A9-02163E01351F.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/750/00000/28938773-BD72-E511-A479-02163E01432A.root'
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/301A497D-70B0-E511-9630-002590D0AFA8.root',
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/7210C351-67B0-E511-A34C-7845C4FC37AF.root',
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/745E2A4F-67B0-E511-9DA3-0090FAA57620.root',
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/7E46D250-67B0-E511-BB96-0025905C3E66.root'
 )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v6')
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v32')
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.MessageLogger = cms.Service("MessageLogger",  
 cout = cms.untracked.PSet(  
  default = cms.untracked.PSet( ## kill all messages in the log  
 
   limit = cms.untracked.int32(0)  
  ),  
  FwkJob = cms.untracked.PSet( ## but FwkJob category - those unlimitted  
 
   limit = cms.untracked.int32(-1)  
  )  
 ),  
 categories = cms.untracked.vstring('FwkJob'), 
 destinations = cms.untracked.vstring('cout')  
)  

process.TFileService=cms.Service("TFileService",
    fileName=cms.string("Test_QCD_data.root")
)
print "test1"
process.analyzeBasicPat = cms.EDAnalyzer("Trigger",
        bits = cms.InputTag("TriggerResults","","HLT"),
        objects = cms.InputTag("slimmedPatTrigger"),
        prescales = cms.InputTag("patTrigger"),
	LHEEventProductInputTag   = cms.InputTag('externalLHEProducer'),
        LHERunInfoProductInputTag = cms.InputTag('externalLHEProducer'),
#      objects = cms.InputTag("selectedPatTrigger"),
	RootFileName = cms.untracked.string('pythia8_test_13tev.root'),

      
 )


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
process.p = cms.Path(process.analyzeBasicPat)
print "test2"
#process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)



