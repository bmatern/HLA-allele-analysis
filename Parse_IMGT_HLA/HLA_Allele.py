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

import os
import time
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def seqForwardComplement(s):
    origSeq = Seq(s, IUPAC.ambiguous_dna)
    return str(origSeq.complement())

def seqReverseComplement(s):    
    origSeq = Seq(s, IUPAC.ambiguous_dna)
    return str(origSeq.reverse_complement())

def seqReverse(s):
    #Extended Slice Syntax.  Not sure what that means but step size is -1
    return s[::-1]

# Check if a string is a number.
def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    

def readFirstSeqFromFasta(fileName):
    try:

        #time.sleep(1)
        #print('Im about to open a sequence file.')
        #print ('input file is ' + fileName)
        #time.sleep(1)

        #currentWorkingDirectory = os.path.dirname(os.path.realpath(__file__))  
        #print('current directory:'  + currentWorkingDirectory)     
        #time.sleep(1) 

        #fullInputPath = os.path.join(self.outputDirectory, fileName)
       # print('FullPath:' + fullInputPath)

        #time.sleep(1)
        fastaHandler = open(fileName, "rU")
        parsedFile = SeqIO.parse(fastaHandler, "fasta")
        firstRecord = parsedFile.next()
        #print ('SequenceName' + firstRecord.id)

        #print ('Sequence:{' + str(firstRecord.seq) + '}')
        #time.sleep(5)
        
        return str(firstRecord.seq)
        #for record in :SeqIO.parse(fileName, "fasta")
#...     print("%s %i" % (record.id, len(record)))
    except Exception as e:

        print('EXCEPTION READING SEQUENCE FROM FASTA:' + str(e))


class HLA_Allele:
    def __init__(self):
    
        self.alleleName = ''
        self.allelePrefix = ''
        self.geneName = ''
        self.nomenclatureGroupCount = 0
        self.alleleGroup = ''
        self.specificProtein = ''
        self.synonymousSubstitution = ''
        self.noncodingChanges = ''
        self.expressionSuffix = ''
        self.sequence = ''
        #TODO: Make APDSequence a method or whatever.  Done.
        #APDSequence = ''
        #TODO: Get rid of in2 sequnece.  It should just be in the list of feature sequences. Done.
        #TODO: Features can be a dictionary, with lookup and sequnce.  Done  
        #in2Sequence = ''
        #TODO: I think I should make these feature lists a quick dictionary.  Done. 
        self.featuresInFullSequence = {}
        self.featuresInAPDSequence = {}
        self.geneFilter = ''
        self.outputDirectory = ''
    
    def getFastaHeader(self, printAPDSequences):
        currentFastaHeader = self.alleleName + ' ('
        
        if(printAPDSequences):
            currentFeatures = self.featuresInAPDSequence.copy()
            
        else:
            currentFeatures = self.featuresInFullSequence.copy()
            
        for featureKey in sorted(currentFeatures.keys()):
            shortFeatureName = str(featureKey.replace('Simulated ','SIM-').replace('\' ','').replace('Intron','IN').replace('Exon','EX').replace(' ','_')) 
            currentFastaHeader += (shortFeatureName + ', ')
            
        currentFastaHeader = currentFastaHeader[0:len(currentFastaHeader)-2] + ')'


        return currentFastaHeader
    
    def getFastaSequence(self, printAPDSequences):
        if(printAPDSequences):
            return self.APDSequence()            
        else:
            return self.sequence  

    
    def APDSequence(self):
             
        currentAPDSequence = ''
        
        if('Exon 2' in self.featuresInAPDSequence and
           'Intron 2' in self.featuresInAPDSequence and
           'Exon 3' in self.featuresInAPDSequence):
            currentAPDSequence = (
                str(self.featuresInAPDSequence['Exon 2']) +
                str(self.featuresInAPDSequence['Intron 2']) +
                str(self.featuresInAPDSequence['Exon 3']))
            
        elif('Exon 2' in self.featuresInAPDSequence and
           'Simulated Intron 2' in self.featuresInAPDSequence and
           'Exon 3' in self.featuresInAPDSequence):
            currentAPDSequence = (
                str(self.featuresInAPDSequence['Exon 2']) +
                str(self.featuresInAPDSequence['Simulated Intron 2']) +
                str(self.featuresInAPDSequence['Exon 3']))
        
        else:
            print('I cannot construct a reliable APD sequence for this allele.' + self.alleleName)
            for featureKey in self.featuresInAPDSequence.keys():
                currentAPDSequence += str(self.featuresInAPDSequence[featureKey])
        
        return currentAPDSequence

    def copy(self):
        copyAllele = HLA_Allele()
        copyAllele.allelePrefix = self.allelePrefix
        copyAllele.geneName = self.geneName
        copyAllele.nomenclatureGroupCount = self.nomenclatureGroupCount
        copyAllele.alleleGroup = self.alleleGroup
        copyAllele.specificProtein = self.specificProtein
        copyAllele.synonymousSubstitution = self.synonymousSubstitution
        copyAllele.noncodingChanges = self.noncodingChanges
        copyAllele.expressionSuffix = self.expressionSuffix
        copyAllele.featuresInFullSequence = self.featuresInFullSequence.copy()
        copyAllele.featuresInAPDSequence = self.featuresInAPDSequence.copy()
        copyAllele.geneFilter = self.geneFilter
        copyAllele.alleleName = self.alleleName
        #copyAllele.APDSequence = self.APDSequence
        copyAllele.sequence = self.sequence
        #copyAllele.in2Sequence = self.in2Sequence
        copyAllele.outputDirectory = self.outputDirectory
        return copyAllele

    def reverseComplement(self):     
        newAllele = self.copy()
        newAllele.sequence = seqReverseComplement(self.sequence)
        
        for featureKey in self.featuresInFullSequence.keys():
            newAllele.featuresInFullSequence[featureKey] = seqReverseComplement(self.featuresInFullSequence[featureKey])
        
        for featureKey in self.featuresInAPDSequence.keys():
            newAllele.featuresInAPDSequence[featureKey] = seqReverseComplement(self.featuresInAPDSequence[featureKey])
  
        return newAllele

    def forwardComplement(self):     
        newAllele = self.copy()
        newAllele.sequence = seqForwardComplement(self.sequence)
        
        for featureKey in self.featuresInFullSequence.keys():
            newAllele.featuresInFullSequence[featureKey] = seqForwardComplement(self.featuresInFullSequence[featureKey])
        
        for featureKey in self.featuresInAPDSequence.keys():
            newAllele.featuresInAPDSequence[featureKey] = seqForwardComplement(self.featuresInAPDSequence[featureKey])
  
        return newAllele

    def reverse(self):     
        newAllele = self.copy()
        newAllele.sequence = seqReverse(self.sequence)
        
        for featureKey in self.featuresInFullSequence.keys():
            newAllele.featuresInFullSequence[featureKey] = seqReverse(self.featuresInFullSequence[featureKey])
        
        for featureKey in self.featuresInAPDSequence.keys():
            newAllele.featuresInAPDSequence[featureKey] = seqReverse(self.featuresInAPDSequence[featureKey])
  
        return newAllele


    def parseAlleleNomenclatureFiltered(self, alleleNode, geneFilter):
        self.geneFilter = geneFilter
        self.parseAlleleNomenclature(alleleNode)

    def parseAlleleNomenclature(self, alleleNode):

        self.alleleName = alleleNode.get('name')
     
        # Split the HLA allele apart.  See this as a reference:
        # http://hla.alleles.org/nomenclature/naming.html
        alleleSplit = str.split(self.alleleName, '*')
        
        #HLA
        self.allelePrefix  = str.split(alleleSplit[0],'-')[0]
        #A
        self.geneName      = str.split(alleleSplit[0],'-')[1]

        #01:01:01:01
        nomenclatureFields = str.split(alleleSplit[1],':')

        fieldCount = len(nomenclatureFields)
        self.nomenclatureGroupCount = fieldCount

        if (fieldCount == 0):
            print ('For allele ' + self.alleleName + ' no nomenclature fields were found.  This might be a problem.')
        if (fieldCount > 0):
            self.alleleGroup = nomenclatureFields[0]
        if (fieldCount > 1):
            self.specificProtein = nomenclatureFields[1]
        if (fieldCount > 2):
            self.synonymousSubstitution = nomenclatureFields[2]
        if (fieldCount > 3):
            #Testing for suffix at the end of an allele:
            #Is the last character a digit?
            lastToken = nomenclatureFields[3] 
            if (isNumber(lastToken[ len(lastToken)-1 : len(lastToken) ])):
                self.noncodingChanges = lastToken
            else:
                self.noncodingChanges = lastToken[ 0 : len(lastToken)-1 ]
                self.expressionSuffix = lastToken[ len(lastToken)-1 : len(lastToken) ]               
        if (fieldCount > 4):
            print ('For allele ' + self.alleleName + 
                ' I found more than 4 allele groups.  This seems like a nomenclature problem.')
        
        # each 'allele' node should have a 'sequence' underneath it.  I think just one.
        sequenceList = alleleNode.findall('{http://hla.alleles.org/xml}sequence')
        
        # Is it really just one?
        if (len(sequenceList) != 1):
            print( 'Error: More than one sequence node found, thats strange:' + str(len(sequenceList)))
            raise Exception('More than one sequence node found for an allele.' + str(self.alleleName))

        else:            
            for sequenceNode in sequenceList:
                nucSequenceList = sequenceNode.findall('{http://hla.alleles.org/xml}nucsequence')

                # Is it really just one nucleotide sequence?
                if (len(nucSequenceList) != 1):
                    print ('Error: More than one nucsequence node found: '  + str(len(nucSequenceList)))
                    raise Exception('Error: More than one nucsequence node on ' + str(self.alleleName) + ' found: '  + str(len(nucSequenceList)))
                else:
                    self.sequence = nucSequenceList[0].text
                    
                    # This is simple.  Don't filter on gene name, no point to it really.
                    if (len(self.geneFilter)==0 or self.geneName in self.geneFilter):

                        # Calculate APD sequence, which is the sequence including exon 2, intron 2, exon 3.
                        # The Antigen Presenting Domain 
                        featureList = sequenceNode.findall('{http://hla.alleles.org/xml}feature')

                        #TODO Make this generic for all features.
                        #ex2Start = 0
                        #ex2End = 0
                        #in2Start = 0
                        #in2End = 0
                        #ex3Start = 0
                        #ex3End = 0

                        print('Allele:' + self.alleleName)

                        print('Features found: ' + str(len(featureList)))


                        for featureNode in featureList:                        
                            featureName = featureNode.get('name')
                            featureType = featureNode.get('featuretype')

                            # Record what features are present
                            # This excludes the translated protein feature, we just want UTRs,Introns,Exons
                            if (featureType in ['Intron','Exon','UTR']):
                                #featureShortName = (featureName
                                #    .replace('Exon ', 'X')
                                #    .replace('Intron ', 'I')
                                #   .replace('5\' UTR', '5UTR')
                                #    .replace('3\' UTR', '3UTR'))
                                #self.featuresInFullSequence += featureShortName + ' '

                                coordinate = featureNode.findall('{http://hla.alleles.org/xml}SequenceCoordinates')[0]
                                
                                # Subtract one from the start index because 
                                # IMGT XML uses 1-based indexing
                                # Python string indexing uses 0 based.        
                                featureStart = int(coordinate.get('start')) - 1
                                featureEnd = int(coordinate.get('end'))
                                
                                print('Storing this feature:' + str(featureName) + ' which is located here:(' + str(featureStart) + ':' + str(featureEnd) + ')')
                              
                                
                                featureSequence = self.sequence[featureStart:featureEnd]
                                
                                print('sequence=' + featureSequence)
                                #TODO: this is broken.  Parse the features better.  Done I think.
                                
                                self.featuresInFullSequence[featureName] = featureSequence
                                
                                #if(featureName=='Exon 2'):
                                #    ex2Start = int(coordinate.get('start'))
                                #    ex2End = int(coordinate.get('end'))
                                #elif(featureName=='Intron 2'):
                                #    in2Start = int(coordinate.get('start'))
                                #    in2End = int(coordinate.get('end'))
                                    
                                    #We found Intron 2.  I use this information later to generate an Intron 2 reference.  
                                    #Store it away.
                                #    self.in2Sequence = self.sequence[in2Start-1:in2End]
                                #elif(featureName=='Exon 3'):
                                #    ex3Start = int(coordinate.get('start'))
                                #    ex3End = int(coordinate.get('end'))
                                #else:
                                #    pass

                        # If we would like to add this information to the title
                        # self.alleleName += ' ' + self.featuresInFullSequence
                        
                        #TODO: we can check this without using specific indexes. Done.

                        #if (ex2Start != 0 and ex2End != 0 and 
                        #    in2Start != 0 and in2End != 0 and 
                        #   ex3Start != 0 and ex3End != 0):  
                        if('Exon 2' in self.featuresInFullSequence and
                           'Intron 2' in self.featuresInFullSequence and
                           'Exon 3' in self.featuresInFullSequence):
                     
                            print('Ex2-In2-Ex3 Found!')
                            #TODO: Switch this, im not using APD sequence anymore.  Done.
                            #self.APDSequence = (self.sequence[ex2Start-1:ex2End] + 
                            #    self.sequence[in2Start-1:in2End] + 
                            #    self.sequence[ex3Start-1:ex3End])
                            #self.featuresInAPDSequence = 'X2 I2 X3'
                            
                            self.featuresInAPDSequence['Exon 2']   = self.featuresInFullSequence['Exon 2']
                            self.featuresInAPDSequence['Intron 2'] = self.featuresInFullSequence['Intron 2']
                            self.featuresInAPDSequence['Exon 3']   = self.featuresInFullSequence['Exon 3']
                            
                            print('Lengths Ex2:In2:Ex3:Total ' + str(len(self.featuresInAPDSequence['Exon 2'])) 
                                  + ':' + str(len(self.featuresInAPDSequence['Intron 2'])) 
                                  + ':' + str(len(self.featuresInAPDSequence['Exon 3']))
                                  + ':' + str(len(self.APDSequence()))) 

                        #elif (ex2Start != 0 and ex2End != 0 and 
                        #    ex3Start != 0 and ex3End != 0): 
                        elif('Exon 2' in self.featuresInFullSequence and
                           'Exon 3' in self.featuresInFullSequence):
 
                            print('Ex2 and Ex3 found, but no intron between them.  I will replace In2 with a simulated intron.')
                            intronSequence = self.loadConsensusIntron2Sequence()

                            #print ('***FoundIntronSequence' + str(intronSequence) + ' GeneName was:' + self.geneName)

                            #self.APDSequence = (str(self.sequence[ex2Start-1:ex2End]) + 
                            #    str(intronSequence) + 
                            #    str(self.sequence[ex3Start-1:ex3End]))
                            
                            self.featuresInAPDSequence['Exon 2']   = self.featuresInFullSequence['Exon 2']
                            self.featuresInAPDSequence['Simulated Intron 2'] = intronSequence
                            self.featuresInAPDSequence['Exon 3']   = self.featuresInFullSequence['Exon 3']
                            
                            #self.featuresInAPDSequence = 'X2 SIM-I2 X3'
                            
                            print('Lengths Ex2:Sim-In2:Ex3:Total ' + str(len(self.featuresInAPDSequence['Exon 2'])) 
                                  + ':' + str(len(self.featuresInAPDSequence['Simulated Intron 2'])) 
                                  + ':' + str(len(self.featuresInAPDSequence['Exon 3']))
                                  + ':' + str(len(self.APDSequence()))) 


                        else:
                            # TODO  Add some more logic in here.  What can be done?
                            #print('Problem with these coordinates:' + str(ex2Start) + ',' + str(ex2End) + ',' + 
                            #    str(in2Start) + ',' + str(in2End) + ',' + str(ex3Start) + ',' + str(ex3End))
                            #raise Exception('Cannot find Exons 2 and 3 for allele ' + self.alleleName) 
                            print('Cannot find Ex2-In2-Ex3 for allele ' + self.alleleName + 'Guess I will include whatever I have.') 
                            #self.APDSequence = self.sequence
                            self.featuresInAPDSequence = self.featuresInFullSequence
                            print('Lengths Total :' + str(len(self.APDSequence()))) 


                    # Not filtering on gene name anymore.
                    else:
                        print('Not parsing the feature information, because this gene isn\'t in the gene filter.')

    def loadConsensusIntron2Sequence(self):
        #print('Loading an Intron 2 consensus.')
        
        intronSequence = ''
        
        groupwiseIntron2ReferenceFilename = os.path.join('Intron2References', 'HLA_' + self.geneName + '_' + self.alleleGroup + '.consensus.fasta')
        genewiseIntron2ReferenceFilename = os.path.join('Intron2References', 'HLA_' + self.geneName + '.consensus.fasta')

        groupwiseIntron2ReferenceFilenameFull = os.path.join(self.outputDirectory, groupwiseIntron2ReferenceFilename)
        genewiseIntron2ReferenceFilenameFull = os.path.join(self.outputDirectory, genewiseIntron2ReferenceFilename)

        #print('Groupwise Intron 2 Reference Filename:' + groupwiseIntron2ReferenceFilenameFull)
        #print('Genewise Intron 2 Reference Filename:' + genewiseIntron2ReferenceFilenameFull)
        
        #If the groupwise file exists, use that.
        if (os.path.isfile(groupwiseIntron2ReferenceFilenameFull)): 
            print('Groupwise intron 2 reference found.  Using: ' + groupwiseIntron2ReferenceFilename)
            intronSequence = readFirstSeqFromFasta(groupwiseIntron2ReferenceFilenameFull)
        #Otherwise, use the genewise reference.  It's a little less specific than the groupwise consensus
        elif(os.path.isfile(genewiseIntron2ReferenceFilenameFull)):
            print('Groupwise intron 2 reference is missing, using the genewise reference: ' + genewiseIntron2ReferenceFilename)
            intronSequence = readFirstSeqFromFasta(genewiseIntron2ReferenceFilenameFull)
        #I don't know what to do here, I'll use a bunch of ambiguous nucleotides.  Who has a better idea?
        else:
            print('Suitable reference was not found.  I will fill with ambiguous nucleotides.')
            intronSequence = ('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + 
                            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + 
                            'NNNNNNNNNNNNNNNNNNNNNNNN')

        return intronSequence 
