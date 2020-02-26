import os
import datetime
import __main__
timestamp = str(datetime.datetime.now().timestamp()).replace(".", "")
dataFolder = "/data"
projectFolder = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
inputFolder = os.path.join(dataFolder, "input")
workingFolder = os.path.join(dataFolder, "working")
outputFolder = os.path.join(dataFolder, "output")
sequenceFolder = os.path.join(inputFolder, "sequence")
forwardReads = os.path.join(sequenceFolder, "standard_submitted_R1.fastq")
reverseReads = os.path.join(sequenceFolder, "standard_submitted_R2.fastq")
referenceFolder = os.path.join(os.path.split(__main__.__file__)[0], "reference")
referenceGenome = os.path.join(referenceFolder, "zrCommunityStandard.fa")
#referenceGenome = "/opt/miqScoreShotgun/reference/zrCommunityStandard.fa" #for external debugging
referenceDataFile = os.path.join(referenceFolder, "zrCommunityStandard.json")
referenceDataFileHMW = os.path.join(referenceFolder, "zrCommunityStandardHMW.json")
goodMiqExample = os.path.join(referenceFolder, "goodMiq.json")
badMiqExample = os.path.join(referenceFolder, "badMiq.json")
goodMiqExampleHMW = os.path.join(referenceFolder, "goodMiqHMW.json")
badMiqExampleHMW = os.path.join(referenceFolder, "badMiqHMW.json")
goodMiqExampleBacteriaOnly = os.path.join(referenceFolder, "goodMiqBacteriaOnly.json")
badMiqExampleBacteriaOnly = os.path.join(referenceFolder, "badMiqBacteriaOnly.json")
goodMiqExampleBacteriaOnlyHMW = os.path.join(referenceFolder, "goodMiqBacteriaOnlyHMW.json")
badMiqExampleBacteriaOnlyHMW = os.path.join(referenceFolder, "badMiqBacteriaOnlyHMW.json")
logFile = os.path.join(outputFolder, "dada2.%s.log" %timestamp)
