import os
import logging
import miqScoreShotgunPublicSupport
import defaults.standard as default
import miqScoreNGSReadCountPublic

validApplicationModes = ["PE", "SE", "LONG"]

def getApplicationMode():
    applicationModeParamter = miqScoreShotgunPublicSupport.parameters.environmentParameterParser.EnvParameters()
    applicationModeParamter.addParameter("mode", str, default="PE", externalValidation=True)
    applicationMode = applicationModeParamter.mode.value
    applicationMode = applicationMode.upper()
    if not applicationMode in validApplicationModes:
        raise ValueError("ERROR: Application mode must be one of the following: %s. %s was given." %(validApplicationModes, applicationMode))
    return applicationMode


def getApplicationParametersSE():
    parameters = miqScoreShotgunPublicSupport.parameters.environmentParameterParser.EnvParameters()
    parameters.addParameter("sampleName", str, required=True, externalValidation=True)
    parameters.addParameter("maxReadCount", int, default=default.maxReadCount, lowerBound=0)
    parameters.addParameter("workingFolder", str, default=default.workingFolder, expectedDirectory=True)
    parameters.addParameter("reads", str, default = default.forwardReads, expectedFile=True)
    parameters.addParameter("sequenceFolder", str, default.sequenceFolder, expectedDirectory=True)
    parameters.addParameter("outputFolder", str, default=default.outputFolder, createdDirectory=True)
    parameters.addParameter("referenceGenome", str, default=default.referenceGenome, expectedFile=True)
    parameters.addParameter("fileNamingStandard", str, default="zymo", externalValidation=True)
    parameters.addParameter("referenceDataFile", str, default = default.referenceDataFile, expectedFile=True)
    parameters.addParameter("bacteriaOnly", bool, required=False, default=False)
    if parameters.bacteriaOnly.value:
        print("ANALYZING BACTERIA ONLY")
        defaultBadExample = default.badMiqExampleBacteriaOnly
        defaultGoodExample = default.goodMiqExampleBacteriaOnly
    else:
        defaultBadExample = default.badMiqExample
        defaultGoodExample = default.goodMiqExample
    parameters.addParameter("goodMiqExample", str, default=defaultGoodExample, expectedFile=True)
    parameters.addParameter("badMiqExample", str, default=defaultBadExample, expectedFile=True)
    if not validSampleName(parameters.sampleName.value):
        logger.error("Invalid sample name given: %s" %parameters.sampleName.value)
        raise ValueError("Invalid sample name given: %s" %parameters.sampleName.value)
    parameters.checkCreatedFileStructures()
    return parameters


def getApplicationParametersPE():
    parameters = miqScoreShotgunPublicSupport.parameters.environmentParameterParser.EnvParameters()
    parameters.addParameter("sampleName", str, required=True, externalValidation=True)
    parameters.addParameter("maxReadCount", int, default=default.maxReadCount, lowerBound=0)
    parameters.addParameter("workingFolder", str, default=default.workingFolder, expectedDirectory=True)
    parameters.addParameter("forwardReads", str, default = default.forwardReads, expectedFile=True)
    parameters.addParameter("reverseReads", str, default=default.reverseReads, expectedFile=True)
    parameters.addParameter("sequenceFolder", str, default.sequenceFolder, expectedDirectory=True)
    parameters.addParameter("outputFolder", str, default=default.outputFolder, createdDirectory=True)
    parameters.addParameter("referenceGenome", str, default=default.referenceGenome, expectedFile=True)
    parameters.addParameter("fileNamingStandard", str, default="zymo", externalValidation=True)
    parameters.addParameter("referenceDataFile", str, default=default.referenceDataFile, expectedFile=True)
    parameters.addParameter("bacteriaOnly", bool, required=False, default=False)
    if parameters.bacteriaOnly.value:
        print("ANALYZING BACTERIA ONLY")
        defaultBadExample = default.badMiqExampleBacteriaOnly
        defaultGoodExample = default.goodMiqExampleBacteriaOnly
    else:
        defaultBadExample = default.badMiqExample
        defaultGoodExample = default.goodMiqExample
    parameters.addParameter("goodMiqExample", str, default=defaultGoodExample, expectedFile=True)
    parameters.addParameter("badMiqExample", str, default=defaultBadExample, expectedFile=True)
    if not validSampleName(parameters.sampleName.value):
        logger.error("Invalid sample name given: %s" %parameters.sampleName.value)
        raise ValueError("Invalid sample name given: %s" %parameters.sampleName.value)
    parameters.checkCreatedFileStructures()
    return parameters


def getApplicationParametersLong():
    parameters = miqScoreShotgunPublicSupport.parameters.environmentParameterParser.EnvParameters()
    parameters.addParameter("sampleName", str, required=True, externalValidation=True)
    parameters.addParameter("maxReadCount", int, default=default.maxReadCount, lowerBound=0)
    parameters.addParameter("workingFolder", str, default=default.workingFolder, expectedDirectory=True)
    parameters.addParameter("reads", str, default = default.forwardReads, expectedFile=True)
    parameters.addParameter("sequenceFolder", str, default.sequenceFolder, expectedDirectory=True)
    parameters.addParameter("outputFolder", str, default=default.outputFolder, createdDirectory=True)
    parameters.addParameter("referenceGenome", str, default=default.referenceGenome, expectedFile=True)
    parameters.addParameter("fileNamingStandard", str, default="zymo", externalValidation=True)
    parameters.addParameter("referenceDataFile", str, default=default.referenceDataFileHMW, expectedFile=True)
    parameters.addParameter("bacteriaOnly", bool, required=False, default=False)
    if parameters.bacteriaOnly.value:
        print("ANALYZING BACTERIA ONLY")
        defaultBadExample = default.badMiqExampleBacteriaOnlyHMW
        defaultGoodExample = default.goodMiqExampleBacteriaOnlyHMW
    else:
        defaultBadExample = default.badMiqExampleHMW
        defaultGoodExample = default.goodMiqExampleHMW
    parameters.addParameter("goodMiqExample", str, default=defaultGoodExample, expectedFile=True)
    parameters.addParameter("badMiqExample", str, default=defaultBadExample, expectedFile=True)
    if not validSampleName(parameters.sampleName.value):
        logger.error("Invalid sample name given: %s" %parameters.sampleName.value)
        raise ValueError("Invalid sample name given: %s" %parameters.sampleName.value)
    parameters.checkCreatedFileStructures()
    return parameters


applicationModeParameterTable = {"PE"   : getApplicationParametersPE,
                                 "SE"   : getApplicationParametersSE,
                                 "LONG" : getApplicationParametersLong}


def validateFastqPair(forwardPath:str, reversePath:str):
    readCount = miqScoreShotgunPublicSupport.formatReaders.fastq.fastqHandler.validFastqPair(forwardPath, reversePath)
    if not readCount:
        if readCount == False:
            errorMessage = "Fastq file failed validation checks"
        elif readCount == 0:
            errorMessage = "Fastq files appear to be empty of reads"
        else:
            raise RuntimeError("Fastq validation should only return either False or a non-negative integer. This code should be unreachable and this is a bug.")
        logger.error(errorMessage)
        raise miqScoreShotgunPublicSupport.formatReaders.fastq.fastqHandler.FastqFormatError(errorMessage)
    if parameters.maxReadCount.value > 0:
        if readCount > parameters.maxReadCount.value:
            raise RuntimeError("Max fastq read count exceeded for this sample. Max: %s. Read count: %s" %(parameters.maxReadCount.value, readCount))
    return True


def validateFastqSingle(fastqPath:str):
    readCount = miqScoreShotgunPublicSupport.formatReaders.fastq.fastqHandler.validFastqFile(fastqPath)
    if not readCount:
        if readCount == False:
            errorMessage = "Fastq file failed validation checks"
        elif readCount == 0:
            errorMessage = "Fastq files appear to be empty of reads"
        else:
            raise RuntimeError("Fastq validation should only return either False or a non-negative integer. This code should be unreachable and this is a bug.")
        logger.error(errorMessage)
        raise miqScoreShotgunPublicSupport.formatReaders.fastq.fastqHandler.FastqFormatError(errorMessage)
    if parameters.maxReadCount.value > 0:
        if readCount > parameters.maxReadCount.value:
            raise RuntimeError("Max fastq read count exceeded for this sample. Max: %s. Read count: %s" %(parameters.maxReadCount.value, readCount))
    return True


def validSampleName(name:str):
    validNameCharacters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.-_ "
    invalidStartingCharacters = ".-_ "
    if name[0] in invalidStartingCharacters:
        return False
    for character in name:
        if not character in validNameCharacters:
            return False
    return True


def getLoggingParameters():
    loggingParameters = miqScoreShotgunPublicSupport.parameters.environmentParameterParser.EnvParameters()
    loggingParameters.addParameter("logFile", str, default=default.logFile, createdFile=True)
    loggingParameters.addParameter("logLevel", str, default=default.loggingLevel, logLevel=True)
    loggingParameters.addParameter("streamOff", bool, default=False)
    loggingParameters.addParameter("streamLoglevel", str, default=default.loggingLevel, logLevel=True)
    loggingParameters.addParameter("fileLogLevel", str, default=default.loggingLevel, logLevel=True)
    logFilePath = os.path.split(loggingParameters.logFile.value)[0]
    if not os.path.isdir(logFilePath):
        os.makedirs(logFilePath)
    loggingParameters.checkCreatedFileStructures()
    return loggingParameters


def loadDefaultPackage():
    defaultParameters = miqScoreShotgunPublicSupport.parameters.environmentParameterParser.EnvParameters()
    defaultParameters.addParameter("defaultPackageName", str, default="standard", externalValidation=True)
    return miqScoreShotgunPublicSupport.parameters.defaultParser.loadDefaultModule(defaultParameters.defaultPackageName.value)


def setLogging():
    loggingParameters = getLoggingParameters()
    formatter = logging.Formatter(loggingFormat)
    logStreamHandle = logging.StreamHandler()
    logStreamHandle.setFormatter(formatter)
    if not loggingParameters.streamLogLevel.usingDefaultValue:
        logStreamHandle.setLevel(loggingParameters.streamLogLevel.value)
    else:
        logStreamHandle.setLevel(loggingParameters.logLevel.value)
    logFileHandle = logging.FileHandler(loggingParameters.logFile.value)
    logFileHandle.setFormatter(formatter)
    if not loggingParameters.fileLogLevel.usingDefaultValue:
        logFileHandle.setLevel(loggingParameters.fileLogLevel.value)
    else:
        logFileHandle.setLevel(loggingParameters.logLevel.value)
    logger.addHandler(logFileHandle)
    if not loggingParameters.streamOff:
        logger.addHandler(logStreamHandle)


def getReadLengthsFromFastq(path:str):
    return miqScoreShotgunPublicSupport.formatReaders.fastq.fastqHandler.estimateReadLength(path)


readFatePrintNames = {"filteredReads": "Failed Quality Filter",
                      "unmergedReads": "Failed To Merge",
                      "chimericReads": "Chimeric",
                      "Reference": "Aligned To Reference"}


def analyzeStandardResult(resultTable:dict):
    if not parameters.bacteriaOnly.value:
        analysisMethod = "Genomic"
    else:
        analysisMethod = "GenomicBacteriaOnly"
    referenceData = miqScoreNGSReadCountPublic.referenceHandler.StandardReference(parameters.referenceDataFile.value)
    calculator = miqScoreNGSReadCountPublic.MiqScoreCalculator(referenceData, analysisMethod=analysisMethod, floor=0)
    miqScoreResult = calculator.calculateMiq(resultTable, parameters.sampleName.value)
    miqScoreResult.makeReadFateChart(readFatePrintNames=readFatePrintNames)
    miqScoreResult.makeRadarPlots()
    goodComposition, badComposition = miqScoreNGSReadCountPublic.loadReferenceCompositionFromExampleMiq(parameters.goodMiqExample.value, parameters.badMiqExample.value)
    miqScoreResult.makeCompositionBarPlot(goodComposition, badComposition)
    return miqScoreResult


def saveResult(result:miqScoreNGSReadCountPublic.MiqScoreData):
    outputFilePath = os.path.join(parameters.outputFolder.value, "%s.json" %parameters.sampleName.value)
    print("Output results to %s" %outputFilePath)
    outputFile = open(outputFilePath, 'w')
    outputFile.write(result.jsonOutput())
    outputFile.close()
    return outputFilePath


def generateReportReplacementTable(sampleMiq, goodExampleMiq, badExampleMiq, readFatePrintNames:dict=None):
    readFateTable = miqScoreShotgunPublicSupport.reporting.generateReadFateChartBody(miqScoreShotgunPublicSupport.reporting.generateAbsoluteReadFateCounts(sampleMiq), readFatePrintNames)
    replacementTable = {"SAMPLENAME": sampleMiq.sampleID,
                        "MIQSCORE": str(round(sampleMiq.miqScore)),
                        "READFATETABLE": readFateTable,
                        "READFATECHART": sampleMiq.plots["readFates"],
                        "COMPOSITIONBARPLOT": sampleMiq.plots["compositionPlot"],
                        "GOODRADARPLOTLYSIS": goodExampleMiq.plots["radarPlots"]["Lysis Difficulty"],
                        "SAMPLERADARPLOTLYSIS": sampleMiq.plots["radarPlots"]["Lysis Difficulty"],
                        "BADRADARPLOTLYSIS": badExampleMiq.plots["radarPlots"]["Lysis Difficulty"],
                        "GOODRADARPLOTGC": goodExampleMiq.plots["radarPlots"]["GC Content"],
                        "SAMPLERADARPLOTGC": sampleMiq.plots["radarPlots"]["GC Content"],
                        "BADRADARPLOTGC": badExampleMiq.plots["radarPlots"]["GC Content"]}
    return replacementTable


def generateReport(result:miqScoreNGSReadCountPublic.MiqScoreData):
    referenceData = miqScoreNGSReadCountPublic.referenceHandler.StandardReference(parameters.referenceDataFile.value)
    templateFilePath = os.path.join(os.path.split(os.path.abspath(__file__))[0], "reference", "shotgunReportTemplate.html")
    templateFile = open(templateFilePath, 'r')
    template = templateFile.read()
    templateFile.close()
    if not parameters.bacteriaOnly.value:
        analysisMethod = "Genomic"
    else:
        analysisMethod = "GenomicBacteriaOnly"
    goodMiq, badMiq = miqScoreNGSReadCountPublic.loadExampleData(parameters.goodMiqExample.value, parameters.badMiqExample.value, referenceData, analysisMethod)
    replacementTable = generateReportReplacementTable(result, goodMiq, badMiq, readFatePrintNames=readFatePrintNames)
    report = miqScoreNGSReadCountPublic.reportGeneration.generateReport(template, replacementTable)
    reportFilePath = os.path.join(parameters.outputFolder.value, "%s.html" % parameters.sampleName.value)
    print("Output report to %s" % reportFilePath)
    outputFile = open(reportFilePath, 'w')
    outputFile.write(report)
    outputFile.close()
    return reportFilePath


if __name__ == "__main__":
    default = loadDefaultPackage()
    loggingFormat = "%(levelname)s:%(name)s:%(message)s"
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)  # Do not change this line unless you know exactly what you are doing any why you are doing it. This will mess up logging in a way that can be hard to trace back.
    setLogging()
    applicationMode = getApplicationMode()
    if applicationMode == "PE":
        parameters = getApplicationParametersPE()
    elif applicationMode == "SE":
        parameters = getApplicationParametersSE()
    elif applicationMode == "LONG":
        parameters = getApplicationParametersLong()
    else:
        raise RuntimeError("This should never be reached as the application mode checker should have caught an invalid mode. Please investigate how this happened. Application mode: %s" %applicationMode)
    logger.debug("Starting analysis")
    bamFilePath = os.path.join(parameters.outputFolder.value, "%s.bam" %parameters.sampleName.value)
    if applicationMode == "PE":
        miqScoreShotgunPublicSupport.alignmentAnalysis.bwaHandler.bwaAlignPE(parameters.forwardReads.value, parameters.reverseReads.value, parameters.workingFolder.value, bamFilePath, parameters.referenceGenome.value)
        readTable = miqScoreShotgunPublicSupport.alignmentAnalysis.alignmentAnalysisPE.bamFileProcessor(bamFilePath)
    elif applicationMode == "SE":
        miqScoreShotgunPublicSupport.alignmentAnalysis.bwaHandler.bwaAlignSE(parameters.reads.value, parameters.workingFolder.value, bamFilePath, parameters.referenceGenome.value)
        readTable = miqScoreShotgunPublicSupport.alignmentAnalysis.alignmentAnalysisSE.bamFileProcessor(bamFilePath)
    elif applicationMode == "LONG":
        miqScoreShotgunPublicSupport.alignmentAnalysis.minimap2.minimapAlign(parameters.reads.value, parameters.workingFolder.value, bamFilePath, parameters.referenceGenome.value)
        readTable = miqScoreShotgunPublicSupport.alignmentAnalysis.alignmentAnalysisSE.bamFileProcessor(bamFilePath)
    else:
        raise RuntimeError("This should never be reached as the application mode checker should have caught an invalid mode. Please investigate how this happened. Application mode: %s" %applicationMode)
    standardAnalysisResults = analyzeStandardResult(readTable)
    saveResult(standardAnalysisResults)
    generateReport(standardAnalysisResults)
    exit(0)
