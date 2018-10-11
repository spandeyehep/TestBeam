import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:


options.register('dataFile',
                 #'/home/tquast/tbJune2018_H2/EUDAQ_unpacked/run000223.root',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_june/unpacked/run000313.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 313,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('timingFile',
                 '/afs/cern.ch/work/s/spandey/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/output_files/timingFiles/timing000313.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('outputFile',
                 '/afs/cern.ch/work/s/spandey/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/output_files/pedestals/pedestals_000313.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

# options.register('unpackedFile',
#                  '/afs/cern.ch/work/s/spandey/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/output_files/CMSSW_unpacked/unpacked_000313.root',
#                  VarParsing.VarParsing.multiplicity.singleton,
#                  VarParsing.VarParsing.varType.string,
#                  'Output file where pedestal histograms are stored')


options.register('beamEnergy',
                30,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Beam energy.'
                )

options.register('beamParticlePDGID',
                11,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Beam particles PDG ID.'
                )

options.register('runType',
                 "Beam",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Run type: Pedestal, Beam, Simulation.'
                )

options.register('setupConfiguration',
                18,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'setupConfiguration (1: July - 4: 20 Layers in October in H6A".'
                )

options.register('pedestalHighGainFile',
                 '/afs/cern.ch/work/s/spandey/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/output_files/pedestals/pedestalHG_000313.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('pedestalLowGainFile',
                 '/afs/cern.ch/work/s/spandey/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/output_files/pedestals/pedestalLG_000313.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('noisyChannelsFile',
                 '/afs/cern.ch/work/s/spandey/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/output_files/pedestals/noisyChannels_000313.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('electronicMap',
                 "map_CERN_Hexaboard_June_28Sensors_28EELayers_V1.txt",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')

options.register('NHexaBoards',
                28,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('hgcalLayout',
                 'layerGeom_june2018_h2_28layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')

options.register('adcCalibrations',
                 'hgcal_calibration_June2018_v2.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal ADC to MIP calibration file in HGCal/CondObjects/data/')

options.register('SubtractPedestal',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Subtract the pedestals.'
                )

options.register('MaskNoisyChannels',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )


options.register('layerPositionFile',
                 '/afs/cern.ch/user/s/spandey/work/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/HGCal/CondObjects/data/layer_distances_CERN_Hexaboard_June2018_28Layers_v2.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('isSimulation',
                0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 ''
                )


options.register('outputNtuple',
                 #'/afs/cern.ch/work/s/spandey/public/hgcal/2018_TB/OctoberTB/CMSSW_9_3_0/src/output_files/ntuples/ntuple_000223.root',
                 'test_ntuple_000313.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output Ntuple file Name')


options.register('reportEvery',
                100,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '.'
                )

options.maxEvents = 10000

options.parseArguments()
print options

electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap
hgcalLayout="HGCal/CondObjects/data/%s" % options.hgcalLayout
adcCalibrations="HGCal/CondObjects/data/%s" % options.adcCalibrations
layerPositionFile=options.layerPositionFile

################################
process = cms.Process("CMSSWUnpacker")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

################################
#process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputNtuple))

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')


process.source = cms.Source("HGCalTBEUDAQDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:%s"%options.dataFile),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NSkipEvents=cms.untracked.uint32(0),
                            runNumber=cms.untracked.int32(options.runNumber),
                            beamEnergy=cms.untracked.double(options.beamEnergy),
                            beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),
                            runType=cms.untracked.string(options.runType),
                            setupConfiguration=cms.untracked.uint32(options.setupConfiguration)
)

process.timingfilewriter = cms.EDAnalyzer("HGCalTBTimingFileWriter",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        TimingFilePath=cms.untracked.string(options.timingFile)
)

process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(True),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(electronicMap),
                                         HighGainPedestalFileName=cms.untracked.string(options.pedestalHighGainFile),
                                         LowGainPedestalFileName=cms.untracked.string(options.pedestalLowGainFile),
                                         WriteNoisyChannelsFile=cms.untracked.bool(True),
                                         NoisyChannelsFileName=cms.untracked.string(options.noisyChannelsFile),
                                         NTSForPedestalComputation=cms.untracked.int32(0)
)

process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        GlobalTimestampCollectionName=cms.string("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicMap=cms.untracked.string(electronicMap),
                                        SubtractPedestal=cms.untracked.bool(bool(options.SubtractPedestal)),
                                        MaskNoisyChannels=cms.untracked.bool(bool(options.MaskNoisyChannels)),
                                        HighGainPedestalFileName=cms.untracked.string(options.pedestalHighGainFile),
                                        LowGainPedestalFileName=cms.untracked.string(options.pedestalLowGainFile),
                                        ChannelsToMaskFileName=cms.untracked.string(options.noisyChannelsFile)
)


process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag("rawhitproducer", "HGCALTBRAWHITS"),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        calibrationPerChannel=cms.untracked.bool(True),
                                        ExpectedMaxTimeSample=cms.untracked.int32(3),
                                        MaxADCCut=cms.untracked.double(20)
)


if options.isSimulation==0:
    #rundata_tag = cms.InputTag("wirechamberproducer", "FullRunData" )
    rundata_tag = cms.InputTag("source", "RunData" )
    rechit_tag = cms.InputTag("rechitproducer","HGCALTBRECHITS" )
    #dwc_tag = cms.InputTag("wirechamberproducer","DelayWireChambers" )
    #dwc_track_tag = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" )
else:
    rundata_tag = cms.InputTag("source", "FullRunData" )
    rechit_tag = cms.InputTag("source","HGCALTBRECHITS" )
    #dwc_tag = cms.InputTag("source","DelayWireChambers" )
    #dwc_track_tag = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" )



process.eventdisplay = cms.EDAnalyzer("EventDisplay",
                                RUNDATA = rundata_tag, 
                                HGCALTBRECHITS = rechit_tag,
                                electronicsMap = cms.untracked.string(electronicMap),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                eventsToPlot=cms.vint32(range(1, 4))
                              )


process.rechitntupler = cms.EDAnalyzer("RecHitNtupler",
                                       InputCollection=rechit_tag,
                                       RUNDATA = rundata_tag,
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       layerPositionFile = cms.untracked.string(layerPositionFile),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(True),
                                       MipThreshold=cms.untracked.double(5.0),
                                       NoiseThreshold=cms.untracked.double(0.0)
)


####################################
# # Load the standard sequences
# process.load('HGCal.StandardSequences.LocalReco_cff')
# process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################




# process.p = cms.Path(process.rawhitproducer)

# process.end = cms.EndPath(process.output)


process.p = cms.Path( process.timingfilewriter * process.pedestalplotter*process.rawhitproducer*process.rechitproducer*process.eventdisplay*process.rechitntupler)

# process.output = cms.OutputModule("PoolOutputModule",
#                                   fileName = cms.untracked.string("test.root"),
#                                   outputCommands = cms.untracked.vstring('drop *',
#                                                                          'keep *_*_HGCALTBRAWHITS_*',
#                                                                          'keep *_*_HGCALGLOBALTIMESTAMPS_*',
#                                                                          'keep *_*_RunData_*')
# )

# process.output = cms.OutputModule("PoolOutputModule",
#                                   fileName = cms.untracked.string("test.root")
# )


# process.output = cms.OutputModule("PoolOutputModule",
#                                   #fileName = cms.untracked.string(options.processedFile),
#                                   fileName = cms.untracked.string("test_RECO.root"),
#                                   outputCommands = cms.untracked.vstring('drop *',
#                                                                          'keep *_*_HGCALTBRECHITS_*',
#                                                                          'keep *_*_HGCALGLOBALTIMESTAMPS_*',
#                                                                          'keep *_*_RunData_*')
# )

#process.end = cms.EndPath(process.output)

