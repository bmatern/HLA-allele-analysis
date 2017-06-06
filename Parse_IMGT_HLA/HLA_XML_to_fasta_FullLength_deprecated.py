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

lengthMinCutoff=5000

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not os.path.isdir(tempDir):
        os.mkdir(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

# Print a list of alleles in Fasta format.
def printFasta(alleleList, fileName, printFullSequences):
    outputFile = createOutputFile(fileName)     
    for currentAllele in alleleList:
        if(len(currentAllele.sequence) > lengthMinCutoff):
            outputFile.write('>' + currentAllele.alleleName + '\n')
            if(printFullSequences):
                outputFile.write( currentAllele.sequence + '\n')
            else:
                outputFile.write( currentAllele.shortSequence + '\n')
    outputFile.close()

# Load the XML file, and parse to extract the HLA alleles.
def createAlleleList(inputFileName):
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
        currentAllele.parseAlleleNomenclature(alleleNode)
        alleleFullList.append(currentAllele)
    
    return alleleFullList

# This is a short method used for sorting alleles into allele groups.
def getGroupIndex(alleleGroups, Gene, AlleleGroup):
    for i in range(0,len(alleleGroups)):
        currentAlleleGroup = alleleGroups[i]
        if (currentAlleleGroup.AlleleGroup == AlleleGroup and currentAlleleGroup.Gene == Gene):
            return i
    return -1

# Split the big list of alleles into allele groups.
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

# Output Allele information in Fasta format.
def printAlleleInfo(alleleFullList, outputDirectory):
    print ('Creating a fasta reference for all HLA Alleles:' + join(outputDirectory, 'HLA_Alleles.fasta') )

    # Allele Information Output File
    alleleInfoOutputFilename = join(outputDirectory, 'AlleleInfo.txt')
    alleleInfoOutputFile = createOutputFile(alleleInfoOutputFilename)      

    # Main full allele output file
    outputFullFileName = join(outputDirectory, 'HLA_Alleles.fasta')
    printFasta(alleleFullList, outputFullFileName, True)

    # Short full allele output file
    #outputShortFileName = join(outputDirectory, 'HLA_Alleles_short.fasta')
    #printFasta(alleleShortList, outputShortFileName)

    # Print out an Alignment Script for ClustalW.
    #clustalwAlignmentScriptFileName = join(outputDirectory, 'ClustalwAlignmentScript.sh')
    #print ('Generating ClustalW Alignment Script:' + clustalwAlignmentScriptFileName)
    clustalwOutputDirectory = join(outputDirectory, 'ClustalwAlignments')
    if not os.path.isdir(clustalwOutputDirectory):
        os.mkdir(clustalwOutputDirectory)

    clustalwConsensusOutputDirectory = clustalwOutputDirectory.replace('Alignments','Consensus')
    if not os.path.isdir(clustalwConsensusOutputDirectory):
        os.mkdir(clustalwConsensusOutputDirectory)
    #clustalwAlignmentScriptFile = createOutputFile(clustalwAlignmentScriptFileName)

    # Print out an Alignment Script for Muscle.
    #muscleAlignmentScriptFileName = join(outputDirectory, 'MuscleAlignmentScript.sh')
    #print ('Generating Muscle Alignment Script.' + muscleAlignmentScriptFileName)
    #muscleOutputDirectory = join(outputDirectory, 'MuscleAlignments')
    #if not os.path.isdir(muscleOutputDirectory):
    #    os.mkdir(muscleOutputDirectory)
    #muscleAlignmentScriptFile = createOutputFile(muscleAlignmentScriptFileName)

	#print ('Please execute those scripts to generate alignments.  You may need to install clustalw or muscle.')  

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
            join('AlleleGroups',alleleGroup.FileName))

        clustalwAlignmentOutputFileName = outputGroupFileName.replace(
            '.fasta','.aln').replace('/AlleleGroups/','/ClustalwAlignments/')

        clustalwConsensusOutputFileName = outputGroupFileName.replace('/AlleleGroups/','/ClustalwConsensus/')

        alleleInfoOutputFile.write('Sorted Group fasta: ' + outputGroupFileName + '\n')
        alleleInfoOutputFile.write('Alignment: ' + clustalwAlignmentOutputFileName + '\n')
        alleleInfoOutputFile.write('Consensus: ' + clustalwConsensusOutputFileName + '\n')
        
        # if there is more than one allele in the group
        if (len(alleleGroup.Alleles) > 1):
            print (str(len(alleleGroup.Alleles)) + ' Alleles Found.')

            # Print allele group to a fasta file
            printFasta(alleleGroup.Alleles, outputGroupFileName, False)

            # For now, i want to work with HLA A,B,C
            #if (true):
            if (alleleGroup.Gene in ['A','B','C']):

 #               clustalwAlignmentOutputFileName = outputGroupFileName.replace(
 #                   '.fasta','.aln').replace('/AlleleGroups/','/ClustalwAlignments/')

                # If the alignment does not already exist
                if not (os.path.isfile(clustalwAlignmentOutputFileName)):                

                    clustalwCommandline = ClustalwCommandline("clustalw", infile=outputGroupFileName, outfile=clustalwAlignmentOutputFileName)
                    print ('ClustalW Alignment Commandline:' + str(clustalwCommandline))
                    print ('Performing Clustalw Alignment...')
                    #clustalwAlignmentScriptFile.write(str(clustalwCommandline) + '\n') 
                    #Perform the alignment
                    clustalwCommandline()
                else:
                    print ('Alignment file ' + clustalwAlignmentOutputFileName + ' already exists.  Moving on...')

                # Check to make sure the alignment does exist now.    
                if (os.path.isfile(clustalwAlignmentOutputFileName)):  
                    # If consensus does not exist yet
                    if not (os.path.isfile(clustalwConsensusOutputFileName)):  
                        #Perform the consensus
                        alignmentType = 'clustal'

                        align = AlignIO.read(clustalwAlignmentOutputFileName, alignmentType)
                    
                        #fileName = clustalwAlignmentOutputFileName.split[1]

                        #consensusOutputFileName = join(consensusOutputDirectory
                        #    , fileName.replace('.aln','.fasta') )
      
                        print ('Consensus FileName = ' + clustalwConsensusOutputFileName)
                    
                        #print('Alignment:' + str(align))
                        summary_align = AlignInfo.SummaryInfo(align)
                        #print('Alignment Summary:' + str(summary_align))
                        
                        print ('Finding Dumb Consensus...\n')
                        dumb_consensus = summary_align.dumb_consensus()
                        print ('Finding Gap Consensus...\n')
                        gap_consensus = summary_align.gap_consensus()
                        #print ('Consensus=' + str(consensus))

                        # Print Consensus to fasta.

                        # I can cheat and just create an HLA_Allele object, and print thFat.
                        currentAllele = HLA_Allele()
                        currentAllele.sequence = str(gap_consensus)
                        currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                        printFasta([currentAllele], clustalwConsensusOutputFileName, True)

                        pass
                    else:
                        print ('Consensus file ' + clustalwConsensusOutputFileName + ' already exists.  Moving on...')
                else:
                    print ('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                    #raise Exception('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                    pass

                # Write this specific alignment command to the script
                #muscleAlignmentFileName = outputGroupFileName.replace(
                #    '.fasta','.aln').replace('/AlleleGroups/','/MuscleAlignments/')
                #muscleCommandline = MuscleCommandline(input=outputGroupFileName, out=muscleAlignmentFileName)
                #muscleAlignmentScriptFile.write(str(muscleCommandline) + '\n')  
            else:
                print ('Skipping alignment, because I only care about Class I')

        # There is only one allele in this group.
        else:
            print ('Only one allele found')
            currentGene = alleleGroup.Alleles[0].geneName

            # Only class 1.
            if (currentGene in ['A','B','C']):
                currentAllele = HLA_Allele()
                #currentAllele.sequence = str(consensus)
                currentAllele.sequence = alleleGroup.Alleles[0].sequence
                #currentAllele.alleleName = os.path.basename(consensusOutputFileName).replace('.fasta','')
                currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                printFasta([currentAllele], clustalwConsensusOutputFileName, True)
                #raise Exception('Only one allele found')

            
            
        

    alleleInfoOutputFile.close()
    #clustalwAlignmentScriptFile.close()
    #muscleAlignmentScriptFile.close()

    #They should be executable
    #os.chmod(clustalwAlignmentScriptFileName, 0777)

    #OK lets go ahead and launch the alignment script.
    #print ('Starting the Clustalw alignment script: ' + clustalwAlignmentScriptFileName)
    #subprocess.call(clustalwAlignmentScriptFileName, shell=True)

"""
# Lets try to do this without this method
# Look for alignment files and print them.
def constructConsensusFromAlignments(outputDirectory):
    print('Constructing a Consensus from the Alignments')

    #alignmentList = ['ClustalwAlignments', 'MuscleAlignments']
    alignmentList = ['ClustalwAlignments']

    for alignmentName in alignmentList:
        alignmentOutputSubDir = join(outputDirectory, alignmentName)
        consensusOutputDirectory = alignmentOutputSubDir.replace('Alignments','Consensus')

        if not os.path.isdir(consensusOutputDirectory):
            os.mkdir(consensusOutputDirectory)

        fileNames = [f for f in listdir(alignmentOutputSubDir) if isfile(join(alignmentOutputSubDir, f))]
        fileNames = sorted(fileNames)
        for fileName in fileNames:
            if('.aln' in fileName): 
                combinedFileName = join(alignmentOutputSubDir, fileName)
                print ('Alignment found:' + combinedFileName + '.  Finding Consensus...'  )   
                try:
                    if (alignmentName == 'ClustalwAlignments'):
                        alignmentType = 'clustal'
                    elif(alignmentName == 'MuscleAlignments'):
                        alignmentType = 'muscle'
                    else:
                        print ('Unknown Alignment Type.')
                        raise

                    align = AlignIO.read(combinedFileName, alignmentType)

                    consensusOutputFileName = join(consensusOutputDirectory
                        , fileName.replace('.aln','.fasta') )
  
                    #print ('Consensus FileName = ' + consensusOutputFileName)
                
                    #print('Alignment:' + str(align))
                    summary_align = AlignInfo.SummaryInfo(align)
                    #print('Alignment Summary:' + str(summary_align))
                    
                    #print ('Finding Consensus...\n')
                    consensus = summary_align.dumb_consensus()
                    #consensus = summary_align.gap_consensus()
                    #print ('Consensus=' + str(consensus))

                    # Print Consensus to fasta.
                    # I can cheat and just create an HLA_Allele object, and print thFat.
                    currentAllele = HLA_Allele()
                    currentAllele.sequence = str(consensus)
                    currentAllele.alleleName = os.path.basename(consensusOutputFileName).replace('.fasta','')
                    printFasta([currentAllele], consensusOutputFileName)
                
                except ValueError:
                    print('Looks like ' + fileName 
                        + ' doesn\'t have a proper alignment.'
                        + '  MSA can fail if there is only one sequence to align.'
                        + '  Probably should fix that, since this group is not included in the reference.' )

                except Exception:
                    print('Cannot open alignment.  Debug this, not sure what is wrong.') 
                    traceback.print_exc()  
                    print sys.exc_info()[0]
                    print sys.exc_info()[1]
                    print sys.exc_info()[2]

"""

# Combine the consensus sequences into an HLA Groupwise Reference
def combineGroupConsensusIntoReference(outputDirectory):
    print("Combining Group Consensus\' into a Groupwise Reference")
    #consensusList = ['ClustalwConsensus', 'MuscleConsensus']
    consensusList = ['ClustalwConsensus']

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
                        sequenceList.append(currentAllele)

            printFasta(sequenceList, referenceOutputFileName, True)
        except Exception:
            print ('Unexpected problem when creating HLA Reference for ' + consensusName)
            print sys.exc_info()[0]
            print sys.exc_info()[1]
            print sys.exc_info()[2]
 

  


if __name__=='__main__':
    try:
        #The [0] arg is apparently the name of the script.  Interesting, I guess that makes sense.
        #First arg is the input file.  Second arg is the output directory.
        inputFileName = sys.argv[1]
        outputDirectory = sys.argv[2]
        print('*** Generating a Fasta reference file from a IMGT HLA XML. ***')
        print('Input:' + inputFileName + '\nOutput:' + outputDirectory)
        print('Just a second...')

        # alleleList is all alleles read from the HLA XML file.
        alleleList = createAlleleList(inputFileName)
        printAlleleInfo(alleleList, outputDirectory)
        # Probably right here I should execute the consensus scripts.  
        # If I do that, I should add code to check if they already exist beforehand
        # constructConsensusFromAlignments(outputDirectory)
        

        #combineGroupConsensusIntoReference(outputDirectory)

        print('Done.  Yay.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise
