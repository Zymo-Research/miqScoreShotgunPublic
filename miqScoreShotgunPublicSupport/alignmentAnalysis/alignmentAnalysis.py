import os
import pysam
import typing

'''
Read species calls:
Chimera_like: Both paired ends disagree, no good single-species explanation
Poor_quality: Depending on if call was being made on single or both ends, quality of read(s) for making call was insufficient
Unmapped_reads: No usable alignment data
Ambiguous_reads: Both reads support two or more species without a big enough delta in confidence to be decisive
If a species is given, that was a successful call
Under certain parameter settings, a tuple of multiple species may be given indicating reads that support multiple species
'''


class ReadAlignmentData(object):
    __slots__ = ["qname",
                 "mapq",
                 "secondaryAlignment",
                 "isRead2",
                 "species"]

    def __init__(self, contig:str, mapq:int, qname:str, secondaryAlignment:bool, isRead2:bool):
        self.species = extractSpecies(contig)
        self.mapq = mapq
        self.qname = qname
        self.secondaryAlignment = secondaryAlignment
        self.isRead2 = isRead2

    def __str__(self):
        secondary = ""
        read = "R1"
        if self.isRead2:
            read = "R2"
        if self.secondaryAlignment:
            secondary = ", Secondary"
        return "%s, %s, %s%s, %s " %(self.species, self.mapq, read, secondary, self.qname)


class ReadPairData(object):
    __slots__ = ["forwardRead",
                 "reverseRead",
                 "speciesConflict",
                 "secondaryAlignment",
                 "mapped",
                 "bothMapped",
                 "calledSpecies",
                 "poorQualityDrop"]

    def __init__(self, read1:ReadAlignmentData, read2:ReadAlignmentData):
        forwardRead = None
        reverseRead = None
        for read in [read1, read2]:
            if read1.isRead2:
                reverseRead = read1
            else:
                forwardRead = read1
            if read2.isRead2:
                reverseRead = read2
            else:
                forwardRead = read2
        if not forwardRead and reverseRead:
            raise MispairedReadError("Two reads of the same direction were given: %s and %s" %(read1, read2))
        self.forwardRead = forwardRead
        self.reverseRead = reverseRead
        self.speciesConflict = self.checkSpeciesConflict()
        self.secondaryAlignment = self.checkSecondaryAlignment()

    def checkSpeciesConflict(self, readPairDisputeDecisionThresholdQDelta:int=40, minimumMapqForSingleReadConfidence:int=40, minimumMapqForPairedReadConfidence:int=30):
        if bool(self.forwardRead.species) != bool(self.reverseRead.species): #Functionally a logical XOR, handles cases where one read was unmappable
            self.mapped = True
            self.bothMapped = False
            if self.forwardRead.species:
                if self.forwardRead.mapq < minimumMapqForSingleReadConfidence:
                    self.mapped = False
                    self.bothMapped = False
                    self.calledSpecies = None
                    self.poorQualityDrop = True
                    return False
                self.calledSpecies = (self.forwardRead.species)
            else:
                if self.reverseRead.mapq < minimumMapqForSingleReadConfidence:
                    self.mapped = False
                    self.bothMapped = False
                    self.calledSpecies = None
                    self.poorQualityDrop = True
                    return False
                self.calledSpecies = (self.reverseRead.species)
            return False
        if not self.forwardRead.species and not self.reverseRead.species:
            self.mapped = False
            self.bothMapped = False
        else:
            self.mapped = True
            self.bothMapped = True
        if self.forwardRead.species == self.reverseRead.species:
            if self.forwardRead.mapq < minimumMapqForPairedReadConfidence or self.reverseRead.mapq < minimumMapqForPairedReadConfidence:
                self.calledSpecies = None
                self.poorQualityDrop = True
            else:
                self.poorQualityDrop = False
                self.calledSpecies = (self.forwardRead.species)
            return False
        else:
            mapqs = (self.forwardRead.mapq, self.reverseRead.mapq)
            mapqDelta = abs(mapqs[0] - mapqs[1])
            if max(mapqs) < minimumMapqForSingleReadConfidence:
                self.calledSpecies = None
                self.poorQualityDrop = True
            if mapqDelta >= readPairDisputeDecisionThresholdQDelta:
                if mapqs[0] > mapqs[1]:
                    self.calledSpecies = self.forwardRead.species
                else:
                    self.calledSpecies = self.reverseRead.species
                self.poorQualityDrop = False
                return False
            self.calledSpecies = (self.forwardRead.species, self.reverseRead.species)
            self.poorQualityDrop = False
            return True

    def checkSecondaryAlignment(self):
        return self.forwardRead.secondaryAlignment or self.reverseRead.secondaryAlignment


class MispairedReadError(Exception):
    pass


def extractSpecies(referenceName:str):
    if referenceName is None:
        return None
    return "_".join(referenceName.split("_")[:2])


def readParallelProcessor(read:pysam.AlignedRead):
    return ReadAlignmentData(read.reference_name, read.mapping_quality, read.query_name, read.is_secondary, read.is_read2)


def listBamFiles(folder:str):
    folderFiles = os.listdir(folder)
    bamFiles = [os.path.join(folder, file) for file in folderFiles if file.endswith(".bam")]
    return bamFiles


def generateAnalyzedReadList(bamFilePath:str):
    import datetime
    startTime = datetime.datetime.now()
    bamFile = pysam.AlignmentFile(bamFilePath, "rb")
    analyzedReads = []
    readCount = 0
    for read in bamFile:
        analyzedReads.append(ReadAlignmentData(read.reference_name, read.mapping_quality, read.query_name, read.is_secondary, read.is_read2))
        readCount += 1
        if readCount % 500000 == 0:
            analysisTime = datetime.datetime.now() - startTime
            print("Analyzed %s reads in %s" %(readCount, analysisTime), flush=True)
    bamFile.close()
    analysisTime = datetime.datetime.now() - startTime
    print("Analyzed %s reads in %s" % (readCount, analysisTime))
    return analyzedReads


def readSorter(readList:typing.List[ReadAlignmentData]):
    sortedReads = {}
    readList.insert(0, None)
    read = readList.pop()
    while read is not None:
        if not read.qname in sortedReads:
            sortedReads[read.qname] = []
        sortedReads[read.qname].append(read)
        read = readList.pop()
    return sortedReads


def sortedReadsDictToList(readDict:dict):
    readList = []
    qnames = list(readDict.keys())
    qnames.insert(0, None)
    qname = qnames.pop()
    while qname:
        readList.append(readDict[qname])
        del readDict[qname]
        qname = qnames.pop()
    return readList


def getSpeciesCallCounts(readSetList:typing.List[typing.List[ReadAlignmentData]]):
    import collections
    speciesCalls = [callSpeciesFromReadSet(readSet) for readSet in readSetList]
    speciesCounts = collections.Counter(speciesCalls)
    return speciesCounts


def callSpeciesFromReadSet(readList:typing.List[ReadAlignmentData], reportSpeciesConflictListAsChimeraLike:bool=True):
    if len(readList) == 2:
        readPair = ReadPairData(*readList)
        if readPair.speciesConflict:
            if reportSpeciesConflictListAsChimeraLike:
                return "Chimera_like"
            else:
                return readPair.calledSpecies
        if not readPair.mapped:
            return "Unaligned_reads"
        if readPair.poorQualityDrop:
            return "Poor_quality"
        if not readPair.calledSpecies:
            return "Unaligned_reads"
        return readPair.calledSpecies
    else:
        return multimapDisputeResolution(readList)


def splitForwardAndReverse(readList:typing.List[ReadAlignmentData], removeUnaligned:bool=True):
    forward = []
    reverse = []
    if removeUnaligned:
        forward = [read for read in readList if not read.isRead2 and read.species]
        reverse = [read for read in readList if read.isRead2 and read.species]
    else:
        forward = [read for read in readList if not read.isRead2]
        reverse = [read for read in readList if read.isRead2]
    return forward, reverse


def getConfidentReads(readList:typing.List[ReadAlignmentData], multimapDisputeResolutionQDelta:int=20):
    if len(readList) == 1:
        if readList[0].species:
            return set([readList[0]])
        else:  #Handles a case where one of the sets came back as unaligned
            return set()
    else:
        mapqs = [read.mapq for read in readList]
        topMapq = max(mapqs)
        confidentReads = [read for read in readList if topMapq - read.mapq <= multimapDisputeResolutionQDelta]
        return set(confidentReads)


def getConfidentSpecies(readList:typing.List[ReadAlignmentData], multimapDisputeResolutionQDelta:int=20):
    if len(readList) == 1:
        if readList[0].species:
            return set([readList[0]])
        else:  #Handles a case where one of the sets came back as unaligned
            return set()
    else:
        mapqs = [read.mapq for read in readList]
        topMapq = max(mapqs)
        confidentSpecies = [read.species for read in readList if topMapq - read.mapq <= multimapDisputeResolutionQDelta]
        return set(confidentSpecies)


def getMaxMapqFromReadList(readList:typing.List[ReadAlignmentData]):
    return max([read.mapq for read in readList])


def multimapDisputeResolution(readList:typing.List[ReadAlignmentData], minimumMapqForSingleReadConfidence:int=40, minimumMapqForPairedReadConfidence:int=30):
    forwardReads, reverseReads = splitForwardAndReverse(readList, removeUnaligned=True)
    if not reverseReads:
        forwardMapqMax = getMaxMapqFromReadList(forwardReads)
        if not forwardMapqMax >= minimumMapqForSingleReadConfidence:
            return "Poor_quality"
        confidentForwards = getConfidentReads(forwardReads)
        if len(confidentForwards) == 1:
            confidentForwards = list(confidentForwards)
            return confidentForwards[0].species
        else:
            return "Ambiguous_reads"
    if not forwardReads:
        reverseMapqMax = getMaxMapqFromReadList(reverseReads)
        if not reverseMapqMax >= minimumMapqForSingleReadConfidence:
            return "Poor_quality"
        confidentReverses = getConfidentReads(reverseReads)
        if len(confidentReverses) == 1:
            confidentReverses = list(confidentReverses)
            return confidentReverses[0].species
        else:
            return "Ambiguous_reads"
    forwardMapqMax = getMaxMapqFromReadList(forwardReads)
    reverseMapqMax = getMaxMapqFromReadList(reverseReads)
    if not forwardMapqMax >= minimumMapqForPairedReadConfidence or reverseMapqMax >= minimumMapqForPairedReadConfidence:
        return "Poor_quality"
    forwardConfidentSpecies = getConfidentSpecies(forwardReads)
    reverseConfidentSpecies = getConfidentSpecies(reverseReads)
    speciesIntersection = forwardConfidentSpecies.intersection(reverseConfidentSpecies)
    if len(speciesIntersection) == 0:
        return "Chimera_like"
    elif len(speciesIntersection) == 1:
        speciesIntersection = list(speciesIntersection)
        return speciesIntersection[0]
    else:
        return "Ambiguous_reads"


def bamFileProcessor(bamFile:str):
    print("Starting analysis of %s" %bamFile, flush=True)
    readList = generateAnalyzedReadList(bamFile)
    sortedReads = readSorter(readList)
    del readList
    readSets = sortedReadsDictToList(sortedReads)
    del sortedReads
    speciesCallCounts = getSpeciesCallCounts(readSets)
    return speciesCallCounts