# This file is part of HLA-allele-analysis.
#
# HLA-allele-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HLA-allele-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with HLA-allele-analysis. If not, see <http://www.gnu.org/licenses/>.

import xml.etree.ElementTree
import sys
import os
from os import listdir
from os.path import isfile, join, split
from HLA_Allele import HLA_Allele
from HLA_Allele_Group import HLA_Allele_Group
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
import traceback
import subprocess

#Class I Genes
classIGenes = ['A','B','C']
justBForTesting = ['B']
justEForTesting = ['E']

genesForAnalysis = classIGenes
#genesForAnalysis = justEForTesting
#genesForAnalysis = justBForTesting

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not os.path.isdir(tempDir):
        os.mkdir(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

# Print a list of alleles in Fasta format.
def printFasta(alleleList, fileName, printAPDSequences):
    outputFile = createOutputFile(fileName)     
    for currentAllele in alleleList:    
        

        # Print all sequence information provided in the XML
        #if(not printAPDSequences):
            #TODO fix featuresInFullSequence.  Might work this way.  Done.
            #currentFastaHeader = currentAllele.alleleName + ' (' + currentAllele.featuresInFullSequence + ')'
        #    currentFastaHeader = currentAllele.alleleName + ' (' + currentAllele.featuresInFullSequence + ')'
        #    for featureKey in currentAllele.featuresInAPDSequence.keys():
        #       currentFastaHeader += str(featureKey)
        #    currentSequence = currentAllele.sequence
        # Print just the antigen presenting domain stuff.  
        #else:
            #TODO fix featuresInFullSequence.  Might work this way.  Done.
        #    currentFastaHeader = currentAllele.alleleName + ' (' + currentAllele.featuresInAPDSequence + ')'
        #    currentSequence = currentAllele.APDSequence

        outputFile.write('>' + str(currentAllele.getFastaHeader(printAPDSequences)) + '\n')
        outputFile.write(str(currentAllele.getFastaSequence(printAPDSequences)) + '\n')

    outputFile.close()

# Load the XML file, and parse to extract the HLA alleles.
def createAlleleList(inputFileName, outputDirectory):
    #document root is an "alleles" node.  It has a series of "allele" nodes underneath it.  Nothing else, according to the schema.
    print('Loading input XML file...')
    documentRoot = xml.etree.ElementTree.parse(inputFileName).getroot()
    print(str(len(documentRoot)) + ' alleles found.')

    alleleFullList = []
    # Loop through allele nodes and create allele objects
    for index, alleleNode in enumerate(documentRoot):

        #if (index % 1000 == 0):
        #    This goes quickly, don't need to print a progress bar.
        #    print ('Allele#:' + str(index)) 
        #    pass
        currentAllele = HLA_Allele()
        currentAllele.outputDirectory = outputDirectory
        currentAllele.parseAlleleNomenclature(alleleNode)
        #print ('Assigning output directory:' + outputDirectory)
        alleleFullList.append(currentAllele)
    
    return alleleFullList

# This is a short method used for sorting alleles into allele groups.
def getGroupIndex(alleleGroups, Gene, AlleleGroup):
    for i in range(0,len(alleleGroups)):
        currentAlleleGroup = alleleGroups[i]
        if (currentAlleleGroup.AlleleGroup == AlleleGroup and currentAlleleGroup.Gene == Gene):
            return i
    return -1
def getGeneIndex(alleleGenes, Gene):
    for i in range(0,len(alleleGenes)):
        currentAlleleGene = alleleGenes[i]
        if (currentAlleleGene.Gene == Gene):
            return i
    return -1

# Split the big list of alleles into allele groups.
# In this instance, an Allele Group is a combination of the Gene and the first field of HLA nomenclature.
def getAlleleGroups(alleleFullList):
    print ('Splitting Alleles into Groups.')
    alleleGroups = []   
    for allele in alleleFullList:
        groupIndex = getGroupIndex(alleleGroups, allele.geneName, allele.alleleGroup)

        # No, the group doesn't already exist
        if (groupIndex == -1):

            currentAlleleGroup = HLA_Allele_Group()
            currentAlleleGroup.Alleles = []
            #print('Group should be empty:' + str(len(currentAlleleGroup.Alleles)))
            currentAlleleGroup.Gene = allele.geneName
            currentAlleleGroup.AlleleGroup = allele.alleleGroup
            currentAlleleGroup.Alleles.append(allele)
            currentAlleleGroup.FileName = 'HLA_' + allele.geneName + '_' + allele.alleleGroup + '.fasta'
            alleleGroups.append(currentAlleleGroup)

        # Yes, this group exists already.  Add this to the list.
        else:
            alleleGroups[groupIndex].Alleles.append(allele)

    return alleleGroups

def getAlleleGenes(alleleFullList):
    print ('Splitting Alleles into Genes.')
    alleleGenes = []   
    for allele in alleleFullList:
        geneIndex = getGeneIndex(alleleGenes, allele.geneName)

        # No, the group doesn't already exist
        if (geneIndex == -1):

            currentAlleleGroup = HLA_Allele_Group()
            currentAlleleGroup.Alleles = []
            #print('Group should be empty:' + str(len(currentAlleleGroup.Alleles)))
            currentAlleleGroup.Gene = allele.geneName
            currentAlleleGroup.Alleles.append(allele)
            currentAlleleGroup.FileName = 'HLA_' + allele.geneName + '.fasta'
            alleleGenes.append(currentAlleleGroup)

        # Yes, this group exists already.  Add this to the list.
        else:
            alleleGenes[geneIndex].Alleles.append(allele)

    return alleleGenes

def reverseComplementList(alleleFullList):
    reversedAlleles = []   
    for allele in alleleFullList:
        reversedAlleles.append(allele.reverseComplement())
    return reversedAlleles

def forwardComplementList(alleleFullList):
    forwardAlleles = []   
    for allele in alleleFullList:
        forwardAlleles.append(allele.forwardComplement())
    return forwardAlleles

def reverseList(alleleFullList):
    reversedAlleles = []   
    for allele in alleleFullList:
        reversedAlleles.append(allele.reverse())
    return reversedAlleles


def performClustalWAlignmentsForGroupwiseReference(outputDirectory, alleleFullList, runAlignments):
    print ('Performing clustalW Alignments and finding consensus sequences')

    # Create the output directories for clustalw.
    clustalwOutputDirectory = join(outputDirectory, 'ClustalwAlignmentsAPD')
    if not os.path.isdir(clustalwOutputDirectory):
        os.mkdir(clustalwOutputDirectory)

    clustalwConsensusOutputDirectory = clustalwOutputDirectory.replace('Alignments','Consensus')
    if not os.path.isdir(clustalwConsensusOutputDirectory):
        os.mkdir(clustalwConsensusOutputDirectory)
    #clustalwAlignmentScriptFile = createOutputFile(clustalwAlignmentScriptFileName)


    alleleGroups = getAlleleGroups(alleleFullList)
    for index, alleleGroup in enumerate(alleleGroups):
        
        print('(' + str(index + 1) + '/' + str(len(alleleGroups)) + '): HLA-' + alleleGroup.Gene + '*' + alleleGroup.AlleleGroup)

        if (True):
        #if (alleleGroup.Gene in genesForAnalysis):
            outputGroupFileName = join(outputDirectory, 
                join('AlleleGroupsAPD',alleleGroup.FileName))

            clustalwAlignmentOutputFileName = outputGroupFileName.replace(
                '.fasta','.aln').replace('/AlleleGroupsAPD/','/ClustalwAlignmentsAPD/')

            clustalwConsensusOutputFileName = outputGroupFileName.replace('/AlleleGroupsAPD/','/ClustalwConsensusAPD/')
     
            # If the alignment does not already exist
            if not (os.path.isfile(clustalwAlignmentOutputFileName)):

                # if there is more than one allele in the group
                if (len(alleleGroup.Alleles) > 1):
                    print (str(len(alleleGroup.Alleles)) + ' Alleles Found.')                

                    clustalwCommandline = ClustalwCommandline("clustalw", infile=outputGroupFileName, outfile=clustalwAlignmentOutputFileName)
                    print ('ClustalW Alignment Commandline:\n' + str(clustalwCommandline))

                    if (runAlignments):

                       # print ('Performing Clustalw Alignment...')
                        #clustalwAlignmentScriptFile.write(str(clustalwCommandline) + '\n') 
                        #Perform the alignment
                        clustalwCommandline()

                        if (os.path.isfile(clustalwAlignmentOutputFileName)):  
                            # If consensus does not exist yet
                            if not (os.path.isfile(clustalwConsensusOutputFileName)):  
                                #Perform the consensus
                                alignmentType = 'clustal'

                                align = AlignIO.read(clustalwAlignmentOutputFileName, alignmentType)
                            
                                #print ('Consensus FileName = ' + clustalwConsensusOutputFileName)
                            
                                summary_align = AlignInfo.SummaryInfo(align)
                                
                                dumb_consensus = summary_align.dumb_consensus()
                                #print('LengthDumbConsensus:' + str(len(dumb_consensus)))
                                gap_consensus = summary_align.gap_consensus()
                                #print('LengthGapConsensus:' + str(len(gap_consensus)))
                                #print ('Consensus=' + str(gap_consensus))

                                # Print Consensus to fasta.

                                # I can cheat and just create an HLA_Allele object, and print that.
                                currentAllele = HLA_Allele()
                                # I think I'll use the dumb consensus.  The only difference is that a gap consensus allows gaps.
                                currentAllele.APDSequence = str(dumb_consensus)
                                currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                                currentAllele.outputDirectory = outputDirectory
                                #print ('Consensus2=' + currentAllele.APDSequence)
                                printFasta([currentAllele], clustalwConsensusOutputFileName, True)

                                pass
                            else:
                                print ('Consensus file ' + clustalwConsensusOutputFileName + ' already exists.  Moving on...')
                        else:
                            print ('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                            #raise Exception('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                            pass

                    else:
                        print ('Not running Alignments because you told me not to.')

                # There is only one allele in this group.
                else:
                    print ('Only one allele found')
                    currentGene = alleleGroup.Alleles[0].geneName

                    # Only class 1.
                    #if (currentGene in genesForAnalysis):
                    currentAllele = HLA_Allele()
                    currentAllele.sequence = alleleGroup.Alleles[0].sequence
                    currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                    printFasta([currentAllele], clustalwConsensusOutputFileName, False)
          

            else:
                print ('Alignment file ' + clustalwAlignmentOutputFileName + ' already exists.  Moving on...')

        else:
            print ('Skipping alignment, because this gene isnt included in genesForAnalysis')



# Output Allele information in Fasta format.
def printAlleleGroupsAndInfo(alleleFullList, outputDirectory):
    print ('Creating a fasta reference for all HLA Alleles:' + join(outputDirectory, 'HLA_Alleles.fasta') )

    # Allele Information Output File
    alleleInfoOutputFilename = join(outputDirectory, 'AlleleGroupInfo.txt')
    alleleInfoOutputFile = createOutputFile(alleleInfoOutputFilename)      

    # Main full allele output file
    outputFullFileName = join(outputDirectory, 'HLA_Alleles.fasta')
    printFasta(alleleFullList, outputFullFileName, False)

    # Antigen Presenting Domain allele output file
    outputAPDFileName = join(outputDirectory, 'HLA_Alleles_APD.fasta')
    printFasta(alleleFullList, outputAPDFileName, True)

    # Reverses and Complements
    outputAPDRevComFileName = join(outputDirectory, 'HLA_Alleles_APD_RevCom.fasta')
    printFasta(reverseComplementList(alleleFullList), outputAPDRevComFileName, True)
    outputAPDForComFileName = join(outputDirectory, 'HLA_Alleles_APD_ForCom.fasta')
    printFasta(forwardComplementList(alleleFullList), outputAPDForComFileName, True)
    outputAPDRevFileName = join(outputDirectory, 'HLA_Alleles_APD_Rev.fasta')
    printFasta(reverseList(alleleFullList), outputAPDRevFileName, True)

    alleleInfoOutputFile.write('Allele Extraction\n')
    alleleInfoOutputFile.write('Input File:' + inputFileName + '\n')
    alleleInfoOutputFile.write('Output Directory:' + outputDirectory + '\n')

    # Print outputfiles and info for each allele group.
    print ('Generating output files for each HLA Allele Group')
    alleleGroups = getAlleleGroups(alleleFullList)
    for index, alleleGroup in enumerate(alleleGroups):
        
        print('(' + str(index + 1) + '/' + str(len(alleleGroups)) + '): HLA-' + alleleGroup.Gene + '*' + alleleGroup.AlleleGroup)
        alleleInfoOutputFile.write('HLA-' + alleleGroup.Gene + '*' 
            + alleleGroup.AlleleGroup + ' contains ' + str(len(alleleGroup.Alleles))
            + ' alleles.\n')

        outputGroupFileName = join(outputDirectory, 
            join('AlleleGroupsAPD',alleleGroup.FileName))

        alleleInfoOutputFile.write('Sorted Group fasta: ' + outputGroupFileName + '\n')
        #alleleInfoOutputFile.write('Alignment: ' + clustalwAlignmentOutputFileName + '\n')
        #alleleInfoOutputFile.write('Consensus: ' + clustalwConsensusOutputFileName + '\n')
        
        # if there is more than one allele in the group
        if (len(alleleGroup.Alleles) > 1):
            print (str(len(alleleGroup.Alleles)) + ' Alleles Found.')

            # Print allele group to a fasta file
            printFasta(alleleGroup.Alleles, outputGroupFileName, True)

   

        # There is only one allele in this group.
        else:
            print ('Only one allele found')

            printFasta([alleleGroup.Alleles[0]], outputGroupFileName, True)

    alleleInfoOutputFile.close()

# Combine the consensus sequences into an HLA Groupwise Reference
def combineGroupConsensusIntoReference(outputDirectory):
    print("Combining Group Consensus\' into a Groupwise Reference")
    #consensusList = ['ClustalwConsensus', 'MuscleConsensus']
    consensusList = ['ClustalwConsensusAPD']

    for consensusName in consensusList:
        consensusOutputSubDir = join(outputDirectory, consensusName)

        referenceOutputFileName = join(
            outputDirectory
            , consensusName + '.HLA.Groupwise.Reference.fasta')

        try:
            sequenceList = []
            fileNames = [f for f in listdir(consensusOutputSubDir) if isfile(join(consensusOutputSubDir, f))]
            fileNames = sorted(fileNames)
            for fileName in fileNames:
                if('.fasta' in fileName): 
                    records = SeqIO.parse( join(consensusOutputSubDir, fileName) , "fasta")
                    for index, record in enumerate(records):
                        currentAllele = HLA_Allele()
                        currentAllele.sequence = str(record.seq)
                        currentAllele.alleleName = str(record.id)
                        currentAllele.outputDirectory = outputDirectory
                        sequenceList.append(currentAllele)


            #I should check if this printFasta boolean parameter is right.  These should be just the APD Sequences right?  Should
            #No maybe not, because I'm feeding them the consensus sequence directly into the .sequence allele.
            printFasta(sequenceList, referenceOutputFileName, False)

            
        except Exception:
            print ('Unexpected problem when creating HLA Reference for ' + consensusName)
            print(sys.exc_info()[0])
            print(sys.exc_info()[1])
            print(sys.exc_info()[2])


def generateIntron2Consensus(alleleFullList, outputDirectory):
    for featureName in ['Intron 2']:
        shortFeatureName = featureName.replace(' ', '')
        
        #Im deciding to quit here.  Late enough.  I want to fix this method tomorrow.
        
        print ('Creating a ' + featureName + ' Reference:' + join( join(outputDirectory,shortFeatureName + 'References'), 'HLA_Intron2.fasta') )
        
        intron2Alleles = []
    
        for allele in alleleFullList:
            #TODO fix featuresInFullSequence.  Might work this way.
            if('Intron 2' in allele.featuresInFullSequence):
         
                currentIntron2Allele = allele.copy()
                #TODO I don't know if i'm still gonna use in2Sequence.
                currentIntron2Allele.sequence = allele.in2Sequence
                intron2Alleles.append(currentIntron2Allele)
    
        # Intron 2 output file, for analyizing *just* the intron 2
        outputIn2FileName = join(join(outputDirectory,'Intron2References'), 'HLA_Intron2.fasta')
        printFasta(intron2Alleles, outputIn2FileName, False)
    
    
        # Print outputfiles and info for each allele group.
        print ('Generating output files for each HLA Allele Group')
        alleleGroups = getAlleleGroups(intron2Alleles)
        alleleGenes = getAlleleGenes(intron2Alleles)
        combinedAlleleGroups = alleleGroups + alleleGenes
        for index, alleleGroup in enumerate(combinedAlleleGroups):
            
            print('(' + str(index + 1) + '/' + str(len(combinedAlleleGroups)) + '): ' + alleleGroup.FileName)
    
            outputGroupFileName = join(outputDirectory, 
                join('Intron2References',alleleGroup.FileName))
    
            clustalwAlignmentOutputFileName = outputGroupFileName.replace('.fasta','.aln')
            clustalwConsensusOutputFileName = outputGroupFileName.replace('.fasta','.consensus.fasta')
            # if there is more than one allele in the group
            if (len(alleleGroup.Alleles) > 1):
                print (str(len(alleleGroup.Alleles)) + ' Alleles Found.  Writing to file: ' + outputGroupFileName)
    
                # Print allele group to a fasta file
                # So this should actually be a false, I don't want to use the APD sequence here.
                printFasta(alleleGroup.Alleles, outputGroupFileName, False)
    
                if (not os.path.isfile(clustalwAlignmentOutputFileName)): 
                    clustalwCommandline = ClustalwCommandline("clustalw", infile=outputGroupFileName, outfile=clustalwAlignmentOutputFileName)
                    print ('Performing  ClustalW Alignment : \n' + str(clustalwCommandline))
    
                    #Perform the alignment
                    clustalwCommandline()
        
                    # sanity check to make sure it exists.
                    if (os.path.isfile(clustalwAlignmentOutputFileName)):  
                        # If consensus does not exist yet
                        if not (os.path.isfile(clustalwConsensusOutputFileName)):  
                            #Find consensus
                            alignmentType = 'clustal'    
                            align = AlignIO.read(clustalwAlignmentOutputFileName, alignmentType)
                        
                            print ('Consensus FileName = ' + clustalwConsensusOutputFileName)
                        
                            #print('Alignment:' + str(align))
                            summary_align = AlignInfo.SummaryInfo(align)
    
                            dumb_consensus = summary_align.dumb_consensus()
                            #print('LengthDumbConsensus:' + str(len(dumb_consensus)))
                            gap_consensus = summary_align.gap_consensus()
                            #print('LengthGapConsensus:' + str(len(gap_consensus)))
        
                            # Print Consensus to fasta.    
                            # I can cheat and just create an HLA_Allele object, and print that.
                            currentAllele = HLA_Allele()
                            currentAllele.APDSequence = str(dumb_consensus)
                            currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                            currentAllele.outputDirectory = outputDirectory
                            #print ('Consensus2=' + currentAllele.APDSequence)
                            printFasta([currentAllele], clustalwConsensusOutputFileName, True)
        
                            pass
                        else:
                            print ('Consensus file ' + clustalwConsensusOutputFileName + ' already exists.  Moving on...')
                    else:
                        print ('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                        #raise Exception('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                        pass
    
                else:
                    print('This alignment file ' + clustalwAlignmentOutputFileName + ' already exists.  Moving on...')
                #else:
                #    print ('Not running Alignments because you told me not to.')   
    
            # There is only one allele in this group.
            else:
                print ('Only one allele found. Writing to file: ' + outputGroupFileName)
    
                #writing it out twice, that's kind of silly but whatever.
                printFasta([alleleGroup.Alleles[0]], outputGroupFileName, True)
                printFasta([alleleGroup.Alleles[0]], clustalwConsensusOutputFileName, True)
    
        #alleleInfoOutputFile.close()

def createFeatureReferences(alleleFullList, outputDirectory):
    featuresList = ['5\' UTR', '3\' UTR',
                    'Exon 1', 'Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6', 'Exon 7', 'Exon 8', 
                    'Intron 1', 'Intron 2', 'Intron 3', 'Intron 4', 'Intron 5', 'Intron 6', 'Intron 7'
                    ]
    for featureName in featuresList:
        shortFeatureName = str(featureName.replace('Simulated ','SIM-').replace('\' ','').replace('Intron ','IN').replace('Exon ','EX').replace(' ','_')) 
       
        
        featureAlleles = []
        featureReferenceOutputDirectory = join(outputDirectory,shortFeatureName + '_References')
        
        print ('Creating a ' + featureName + ' Reference:' + join( featureReferenceOutputDirectory, 'HLA_' + shortFeatureName + '.fasta') )
        
    
        for allele in alleleFullList:
            if(featureName in allele.featuresInFullSequence):
         
                currentFeatureAllele = allele.copy()
                currentFeatureAllele.sequence = allele.featuresInFullSequence[featureName]
                featureAlleles.append(currentFeatureAllele)
    
        # Intron 2 output file, for analyizing *just* the intron 2
        outputFeatureFileName = join(featureReferenceOutputDirectory, 'HLA_' + shortFeatureName + '.fasta')
        #outputIn2FileName = join(join(outputDirectory,'Intron2References'), 'HLA_Intron2.fasta')
        printFasta(featureAlleles, outputFeatureFileName, False)
    
    
        # Print outputfiles and info for each allele group.
        print ('Generating output files for each HLA Allele Group')
        alleleGroups = getAlleleGroups(featureAlleles)
        alleleGenes = getAlleleGenes(featureAlleles)
        combinedAlleleGroups = alleleGroups + alleleGenes
        for index, alleleGroup in enumerate(combinedAlleleGroups):
            
            print('(' + str(index + 1) + '/' + str(len(combinedAlleleGroups)) + '): ' + alleleGroup.FileName)
    
            outputGroupFileName = join(featureReferenceOutputDirectory,alleleGroup.FileName)
            
            printFasta(alleleGroup.Alleles, outputGroupFileName, False)




if __name__=='__main__':
    try:
        #The [0] arg is apparently the name of the script.  Interesting, I guess that makes sense.
        #First arg is the input file.  Second arg is the output directory.
        inputFileName = sys.argv[1]
        outputDirectory = sys.argv[2]
        print('*** Generating a Fasta reference file from a IMGT HLA XML. ***')
        print('Input:' + inputFileName + '\nOutput:' + outputDirectory)
        print('Just a second...')


        if not os.path.isdir(outputDirectory):
            os.mkdir(outputDirectory)

        # alleleList is all alleles read from the HLA XML file.
        alleleList = createAlleleList(inputFileName, outputDirectory)
        
        # This is a generic method to just make output files for every gene feature.
        # This shoudl overlap with the generateIntron2Consensus method, so I think I'll need to trim that method 
        # so as not to redo efforts. 
        createFeatureReferences(alleleList, outputDirectory)

        # I need Intron 2 sequences for my APD reference, 
        # and I need to generate a consensus of like alleles for the ones that are missing.
        # This takes a long time.
        
        #generateIntron2Consensus(alleleList, outputDirectory)
        
        # I already ran this line of code, 
        # but I'm doing it again quick so it reloades the alleles with the generated intron 2 sequences.
        # There is certainly a better way to do this.
        
        alleleList = createAlleleList(inputFileName, outputDirectory)
        
        # Output many allele references
        
        printAlleleGroupsAndInfo(alleleList, outputDirectory)
    
        # Align and find consensus for allele groups.
        # This takes a long time.
        
        #performClustalWAlignmentsForGroupwiseReference(outputDirectory, alleleList, True)
        
        # Make a groupwise reference consensus.
        
        #combineGroupConsensusIntoReference(outputDirectory)

        print('Done.  Ben did a great job.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print('Unexpected problem during execution:')
        print(sys.exc_info()[1])
        raise
