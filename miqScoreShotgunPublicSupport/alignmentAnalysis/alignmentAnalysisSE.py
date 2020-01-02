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


def callSpeciesFromReadSet(readList:typing.List[ReadAlignmentData], minimumMapqForSingleReadConfidence:int=40):
    if len(readList) == 1:
       if readList[0].species:
           if readList[0].mapq >= minimumMapqForSingleReadConfidence:
               return readList[0].species
           else:
               return "Poor_quality"
       else:
           return "Unaligned_reads"
    else:
       return multimapDisputeResolution(readList, minimumMapqForSingleReadConfidence)


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


def multimapDisputeResolution(readList:typing.List[ReadAlignmentData], minimumMapqForSingleReadConfidence:int=40):
    mapqMax = getMaxMapqFromReadList(readList)
    if not mapqMax >= minimumMapqForSingleReadConfidence:
        return "Poor_quality"
    confidentReads = getConfidentReads(readList)
    if len(confidentReads) == 1:
        confidentReads = list(confidentReads)
        return confidentReads[0].species
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