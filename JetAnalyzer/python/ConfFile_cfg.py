import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("Events")

# QG likelihood
process.load("JetNtupleProducerTool.JetAnalyzer.QGLikelihood_cfi")
process.load("RecoJets.JetProducers.QGTagger_cfi")
process.QGTagger.srcJets = cms.InputTag("slimmedJets")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")

# File service
process.load("CommonTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName=cms.string("JetNtuple_RunIISummer16_13TeV_MC.root")

# Load up the filelist
filePath=os.environ["CMSSW_BASE"]+"/src/JetNtupleProducerTool/JetAnalyzer/python/"
fileList = FileUtils.loadListFromFile(filePath+"RunIISummer16MiniAODv2_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_filelist.txt")

process.source = cms.Source("PoolSource",
	## Process whole data set
	fileNames = cms.untracked.vstring(*fileList)
	## Process just one file
	# fileNames = cms.untracked.vstring(fileList[0])
)

process.Events = cms.EDAnalyzer("JetAnalyzer",
	saveJetTree = cms.bool(False),
	## jet, PF and generator level collections ##
	jets = cms.InputTag("slimmedJets"),
	pfCands = cms.InputTag("packedPFCandidates"),
	genJets = cms.InputTag("slimmedGenJets"),
	genEventInfo = cms.InputTag("generator"),
	## extra track collections
	extraTracks = cms.InputTag("lostTracks"),
	## gen particles collection
	genParticles = cms.InputTag("prunedGenParticles"),
	## good primary vertices ##
	vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
	confGoodVtxNdof = cms.double(4),
	confGoodVtxZ = cms.double(24),
	confGoodVtxRho = cms.double(2),
	## pileup and rhos ##
	pileupInfo = cms.InputTag("slimmedAddPileupInfo"),
	pfRhoAll = cms.InputTag("fixedGridRhoFastjetAll"),
	pfRhoCentral = cms.InputTag("fixedGridRhoFastjetCentral"),
	pfRhoCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
	pfRhoCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
)

# Choose how many events to process (-1 = all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Report execution progress
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.p = cms.Path(process.QGTagger + process.Events)
