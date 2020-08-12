from __future__ import print_function
import sys, Bio, Bio.SeqUtils, traceback, random, copy
from operator import itemgetter

import nrpcalc
from PyVRNA import PyVRNA
import numpy as np

from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt

def readFasta(filename):
    seqList = []
    with open(filename,'r') as f:
        while True:
            id = f.readline().rstrip()
            seq = f.readline().rstrip()
            if len(id) > 0:
                seqList.append( (id, seq) )
            else:
                break
    return seqList

def revcomp(seq):
    xform = {'G' : 'C', 'C' : 'G', 'A' : 'T', 'T' : 'A'}
    return "".join([xform[letter] for letter in seq[::-1]])
    
def globalFunction_MeltingTemperature(seq, spec):
    
    Tm_settings = spec['Tm_settings']
    Tm = mt.Tm_NN(seq, **Tm_settings)
    if Tm < spec['Tm_target'][0]: return False
    if Tm > spec['Tm_target'][1]: return False
    return True

def globalFunction_GCContent(seq, spec):
    
    if spec['GC_target'][0] <= 0.0 and spec['GC_target'][1] >= 100.0: return True
    GC_content = seq.count("C") + seq.count("G")
    if GC_content < spec['GC_target'][0] * len(seq) / 100: return False
    if GC_content > spec['GC_target'][1] * len(seq) / 100: return False
    return True

def localFunction_excludedSequencesEverywhere(seq, spec):
    for excludedSequence in spec['excludedSequencesEverywhere']:
        pos = seq.find(excludedSequence)
        if pos >= 0:   return (False, pos)
    return (True, None)

def localFunction_excludedSequencesWithExclusions(seq, spec):
    for (excludedSequence, pos_list) in spec['excludedSequencesWithExclusions']:
        pos = seq.find(excludedSequence)
        if pos < 0:
            pass
        elif pos in pos_list:
            pass
        else:
            return (False, pos)
    return (True, None)
            
def createGlobalModelFunction(spec):
    
    fcn_list = []
    if 'GC_target' in spec: fcn_list.append(globalFunction_GCContent)
    if 'Tm_target' in spec: fcn_list.append(globalFunction_MeltingTemperature)
    #add more in order of increasing compute time
    
    if len(fcn_list) == 0: return None
    
    def globalModelFunction(seq, spec):
        for fcn in fcn_list:
            Passed = fcn(seq, spec)
            if Passed:
                pass
            else:
                return False
        return True
    return lambda seq: globalModelFunction(seq, spec)

def createLocalModelFunction(spec):
    fcn_list = []
    if 'excludedSequencesEverywhere' in spec: fcn_list.append(localFunction_excludedSequencesEverywhere)
    if 'excludedSequencesWithExclusions' in spec: fcn_list.append(localFunction_excludedSequencesWithExclusions)
    #add more in order of increasing compute time
    
    if len(fcn_list) == 0: return None
    
    def localModelFunction(seq, spec):
        for fcn in fcn_list:
            (Passed, pos) = fcn(seq, spec)
            if Passed:
                pass
            else:
                return (Passed, pos)
        return (True, None)
    return lambda seq: localModelFunction(seq, spec)
            
class OligoPoolCalculator(object):

    def __init__(self, specs, outputFilename, toleranceDeltaT = 0.2):
        
        self.specs = specs
        self.outputFilename = outputFilename
        self.toleranceDeltaT = toleranceDeltaT
        self.excludeSequences = specs['excludeSequences']
        self.backgroundSequenceList = specs['backgroundSequenceList']
        self.Tm_target = specs['Tm_target']
        self.Tm_settings = specs['Tm_settings']
        self.maxOligoLength = specs['maxOligoLength']
        self.primer_Lmax = 1000
        self.barcode_Lmax = 1000
        
        self.numVariants = 0
        self.parts = {}
        self.selectedPrimers = {}
        self.selectedBarcodes = {}
        self.sequenceVariants = []
    
    def run(self):
        
        self.parseSpecifications()
        self.initializeBackground()
        self.designPrimers()
        self.designBarcodes()
        self.designSequences()
        self.validateSequences()
        self.designOligos()
        self.writeOligoPool()
    
    def parseSpecifications(self):
        numVariants = {}
        
        for (pos, part) in enumerate(self.specs['design']):
            type = part['type']
            if not (type in self.parts): self.parts[type] = []
            self.parts[type].append( (pos, part) )
            
            if type == 'restriction_site':
                self.excludeSequences.append(part['sequence'])
                self.backgroundSequenceList.append(part['sequence'])
            
            elif type == 'primer_binding_site':
                self.primer_Lmax = min(part['Lmax'], self.primer_Lmax)
                
            elif type == 'barcode':
                self.barcode_Lmax = min(part['Lmax'], self.barcode_Lmax)
                
                
            elif type == 'constant':
                self.backgroundSequenceList.append(part['sequence'])
            
            elif type == 'variable':
                numVariants[pos] = len(part['variants'])
                for variant in part['variants']:
                    self.backgroundSequenceList.append(variant['sequence'])
        
        if not (max(numVariants.values()) == min(numVariants.values())):
            raise Exception('Error: The number of variants in each variable region must be the same.')
        else:
            self.numVariants = min(numVariants.values())
            print("INFO: Reading %s variable regions, %s restriction sites, %s primer binding sites, %s barcodes, %s constant regions." % 
            (len(self.parts['variable']), len(self.parts['restriction_site']), len(self.parts['primer_binding_site']), len(self.parts['barcode']), len(self.parts['constant'])))
            print("INFO: Number of Sequence Variants: %s " % self.numVariants)
            
    def initializeBackground(self):
    
        path = 'TEMP_PLYVEL_BACKGROUND1'
        self.primer_background = nrpcalc.background(path, self.primer_Lmax,  verbose=True)
        self.primer_background.drop()
        self.primer_background = nrpcalc.background(path, self.primer_Lmax,  verbose=True)
        self.primer_background.multiadd(self.backgroundSequenceList)
        
        path = 'TEMP_PLYVEL_BACKGROUND2'
        self.barcode_background = nrpcalc.background(path, self.barcode_Lmax,  verbose=True)
        self.barcode_background.drop()
        self.barcode_background = nrpcalc.background(path, self.barcode_Lmax,  verbose=True)
        self.barcode_background.multiadd(self.backgroundSequenceList)
        
        
    
    def designPrimers(self):
        
        self.designed_primers = {}
        energy_model = PyVRNA(dangles=0, gquad=True, parameter_file='dna_mathews2004.par')
        print("INFO: Beginning Design of Primer Binding Sequences.")
        
        modelSpec = {'Tm_target' : self.Tm_target, 'Tm_settings' : self.Tm_settings, 'excludedSequencesEverywhere' : self.excludeSequences}
        
        globalModelFunction = createGlobalModelFunction(modelSpec)
        localModelFunction = createLocalModelFunction(modelSpec)
        
        for (pos, spec) in self.parts['primer_binding_site']:
            structureConstraint = "x" * len(spec['sequenceConstraint'])
            numParts = 100
            
            results =  nrpcalc.maker(seq_constr = spec['sequenceConstraint'], struct_constr = structureConstraint, part_type ='DNA',
                          target_size = numParts, Lmax = spec['Lmax'],
                          internal_repeats=False, background=self.primer_background, struct_type='mfe', seed = None, synth_opt=True,
                          local_model_fn=localModelFunction,
                          global_model_fn=globalModelFunction,
                          jump_count=10, fail_count=1000, output_file=None, verbose=True)
    
            candidates = [seq for (num, seq) in results.items()]
            self.designed_primers[pos] = {}
            self.designed_primers[pos]['sequences'] = candidates
            
            #Calculate Melting Temperatures
            self.designed_primers[pos]['Tm'] = [mt.Tm_NN(seq, **self.Tm_settings) for seq in candidates]
            
            #Calculate Minimum Free Folding Free Energies and Structures
            fold_list = [energy_model.RNAfold(seq) for seq in candidates]
            self.designed_primers[pos]['dG_folding'] = [fold.energy for fold in fold_list]
            self.designed_primers[pos]['structure'] = [fold.structure for fold in fold_list]
            
        
        #Identify Sets of Primer Binding Sites with Similar Melting Temperatures
        self.compatiblePrimers = {}
        self.selectedPrimers = {}
        
        candidate_indexes = {}
        numSteps = np.int( (self.Tm_target[1] - self.Tm_target[0]) / self.toleranceDeltaT) + 1
        
        for Tm_target in np.linspace(self.Tm_target[0], self.Tm_target[1], numSteps):
            self.compatiblePrimers[Tm_target] = {}
            for pos in self.designed_primers.keys():
                candidate_indexes[pos] = [i for (i, Tm) in enumerate(self.designed_primers[pos]['Tm']) if abs(Tm - Tm_target) < self.toleranceDeltaT]
            if all( [len(indexes) > 1 for indexes in candidate_indexes.values()] ):
                for pos in self.designed_primers.keys():
                    self.compatiblePrimers[Tm_target][pos] = {'sequences' :  itemgetter(*candidate_indexes[pos])(self.designed_primers[pos]['sequences']),
                                                              'Tm' :         itemgetter(*candidate_indexes[pos])(self.designed_primers[pos]['Tm']),
                                                              'dG_folding' : itemgetter(*candidate_indexes[pos])(self.designed_primers[pos]['dG_folding']),
                                                              'structure'  : itemgetter(*candidate_indexes[pos])(self.designed_primers[pos]['structure'])
                                                             }
        
        print("INFO: Compatible Primer Binding Sites:")
        for (Tm, posObj) in self.compatiblePrimers.items():
            for (pos, info) in posObj.items():
                print('INFO: Target Tm: %s. Position: %s. Primers Available: %s.' % (Tm, pos, len(info['Tm'])))
        
        
        #Select the First Compatible Primer with the Highest Melting Temperature
        highest_Tm = sorted([Tm for (Tm, primers) in self.compatiblePrimers.items() if len(primers) > 0], reverse=True)[0]
        for pos in self.compatiblePrimers[highest_Tm].keys():
            self.selectedPrimers[pos] = {'sequence' : self.compatiblePrimers[highest_Tm][pos]['sequences'][0],
                                         'Tm' : self.compatiblePrimers[highest_Tm][pos]['Tm'][0],
                                         'dG_folding' : self.compatiblePrimers[highest_Tm][pos]['dG_folding'][0],
                                         'structure' : self.compatiblePrimers[highest_Tm][pos]['structure'][0]
                                        }
            self.primer_background.multiadd(self.selectedPrimers[pos]['sequence'])
            self.barcode_background.multiadd(self.selectedPrimers[pos]['sequence'])
                
        print("INFO: Finished Design of Primer Binding Sites using highest Tm = %s." % highest_Tm)
        return
    
    def designBarcodes(self):
                 
        print("INFO: Beginning the Design of Barcode Sequences")
        self.selectedBarcodes = {}
        
        modelSpec = {'excludedSequencesEverywhere' : self.excludeSequences}
        globalModelFunction = createGlobalModelFunction(modelSpec)
        localModelFunction = createLocalModelFunction(modelSpec)
            
        for (pos, spec) in self.parts['barcode']:    
            structureConstraint = "." * len(spec['sequenceConstraint'])
            numParts = self.numVariants
        
            results =  nrpcalc.maker(seq_constr = spec['sequenceConstraint'], struct_constr = structureConstraint, part_type ='DNA',
                              target_size = numParts, Lmax = spec['Lmax'],
                              internal_repeats=False, background=self.barcode_background, struct_type='mfe', seed = None, synth_opt=True,
                              local_model_fn=localModelFunction,
                              global_model_fn=globalModelFunction,
                              jump_count=10, fail_count=1000, output_file=None, verbose=True)
        
            barcodes = [seq for (num, seq) in results.items()]
            if len(barcodes) < numParts: raise Exception('Error: Only %s barcodes were designed, but %s were needed. Adjust the barcode design specifications to resolve.' % (len(barcodes), numParts))
            self.selectedBarcodes[pos] = barcodes
            self.barcode_background.multiadd(barcodes)
        
        print("INFO: Barcode Design completed. %s Barcodes Designed." % numParts)
        return
        
    def designSequences(self):
        
        sequenceVariants = []
        for n in range(self.numVariants):
            seqInfo = {}
            seq = ''
            for (pos, part) in enumerate(self.specs['design']):
                type = part['type']
        
                if type == 'restriction_site':
                    seqInfo[pos] = part
                    seqInfo[pos]['sequence'] = part['sequence']
                    seq += part['sequence']
                
                elif type == 'primer_binding_site':
                    seqInfo[pos] = self.selectedPrimers[pos]
                    if part['direction'] == 'forward':
                        seqInfo[pos]['sequence'] = self.selectedPrimers[pos]['sequence']
                        seq += self.selectedPrimers[pos]['sequence']
                    elif part['direction'] == 'reverse':
                        seqInfo[pos]['sequence'] = revcomp(self.selectedPrimers[pos]['sequence'])
                        seq += revcomp(self.selectedPrimers[pos]['sequence'])
                    
                elif type == 'barcode':
                    seqInfo[pos] = part
                    seqInfo[pos]['sequence'] = self.selectedBarcodes[pos][n]
                    seq += self.selectedBarcodes[pos][n]
                    
                    
                elif type == 'constant':
                    seqInfo[pos] = part
                    seqInfo[pos]['sequence'] = part['sequence']
                    seq +=  part['sequence']
                    
                elif type == 'variable':
                    seqInfo[pos] = part['variants'][n]
                    seqInfo[pos]['sequence'] = part['variants'][n]['sequence']
                    seq += part['variants'][n]['sequence']
            
            seqInfo['sequence'] = seq
            sequenceVariants.append(copy.deepcopy(seqInfo))
        
        self.sequenceVariants = sequenceVariants[:]
        print('self.sequenceVariants[0:10]:', self.sequenceVariants[0:10])
        return
    
    def validateSequences(self):
    
        expectedSites = {}
        for (pos, part) in self.parts['restriction_site']:
            if part['sequence'] in expectedSites:
                expectedSites[part['sequence']] += 1
            else:
                expectedSites[part['sequence']] = 1
            
        for (num,seqInfo) in enumerate(self.sequenceVariants):
            for (motif,expected_counts) in expectedSites.items():
                counts = seqInfo['sequence'].count(motif)
                if counts == expected_counts:
                    pass
                else:
                    raise Exception('ERROR: Expected %s restriction sites (%s), but found %s sites in sequence id %s (%s)' % (expected_counts, motif, counts, num, seqInfo['sequence']))
    
    def addPadding(self):
        
        print("INFO: Adding Padding to Sequences")
        modelSpec = {'GC_target' : [30.0, 70.0], 'excludedSequencesEverywhere' : self.excludeSequences}
        globalModelFunction = createGlobalModelFunction(modelSpec)
        localModelFunction = createLocalModelFunction(modelSpec)
        
        sequenceVariants = self.sequenceVariants[:]
        paddedSequenceVariants = []
        
        for seqInfo in sequenceVariants:
            paddingLength = self.maxOligoLength - len(seqInfo['sequence'])
            
            while True:
                paddingSequence = "".join([random.choice(['A','G','C','T']) for n in range(paddingLength)])
                (localPassed, pos) = localModelFunction(paddingSequence)
                globalPassed = globalModelFunction(paddingSequence)
                if localPassed and globalPassed: break
            lastPos = len(self.specs['design'])
            seqInfo[lastPos+1] = {'type' : 'pad', 'sequence' : paddingSequence}
            seqInfo['oligos'] = [seqInfo['sequence'] + paddingSequence] #Padding 3' end
            paddedSequenceVariants.append(seqInfo)
        self.sequenceVariants = paddedSequenceVariants[:]
        
        print('Padded self.sequenceVariants[0:10]:', self.sequenceVariants[0:10])
        return
    
    def designOligos(self):
        
        #Do we need to break up sequences into multiple oligos?
        #Do we need to pad sequences?
        
        maxSequenceLength = max( [len(seqInfo['sequence']) for seqInfo in self.sequenceVariants] )
        minSequenceLength = min( [len(seqInfo['sequence']) for seqInfo in self.sequenceVariants] )
        
        if maxSequenceLength > self.maxOligoLength:                                     #Yes, we need to break up into multiple oligos with overlapping regions
            #Add functionality here
            pass
        elif self.specs['enablePadding'] and maxSequenceLength > minSequenceLength:     #Yes, we need to pad sequences
            self.addPadding()
        else:                                                                           #No, each sequence variant is a single unpadded oligonucleotide
            for seqInfo in self.sequenceVariants:
                seqInfo['oligos'] = [seqInfo['sequence']]
            
    
    def writeOligoPool(self, replicates = 1):
    
        import csv
        with open(self.outputFilename + '.csv', 'w') as csvfile:
            
            #Define column names here
            fieldnames = ['id','name','oligo','sequence']
            
            for (pos, part) in enumerate(self.specs['design']):
                fieldnames.append(part['name'])
            
            fieldnames.append('replicate')
            fieldnames.extend( ['Tm degC %s' % part['name'] for (pos, part) in self.parts['primer_binding_site']] )            
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
        
            for rep in range(replicates):
                for (seqNum, seqInfo) in enumerate(self.sequenceVariants):
                    
                    for (n, oligo) in enumerate(seqInfo['oligos']):
                        if n == 0: #first oligo for a sequence variant
                            line = {'id': seqNum,
                                    'oligo' : oligo,
                                    'sequence' : seqInfo['sequence'],
                                    'replicate' : rep
                                   }
                            
                            name_list = [seqInfo[pos]['name'] for (pos, part) in self.parts['variable']]
                            line['name'] = ":".join(name_list)
                            
                            for (pos, part) in self.parts['primer_binding_site']:
                                line['Tm degC %s' % part['name']] = seqInfo[pos]['Tm']
                            
                            for (pos, part) in enumerate(self.specs['design']):
                                line[part['name']] = seqInfo[pos]['sequence']
                            
                            writer.writerow(line)
                                         
                        else:
                            line = {'id': str(seqNum) + "-" + str(n), 'oligo' : oligo}
                            writer.writerow(line)
              
        print("INFO: Oligopool Construction Completed. Saved to %s." % self.outputFilename + '.csv')
        return
        
    
    def calcHammingDistanceDistributions(self):
        print("INFO: Calculating Hamming Distance Distributions.")
        hammingDistanceDistributions = []
        for seqList in self.seqSets:
            barcodeList = [seq['barcode'] for seq in seqList]
            distribution = {}
            for n1 in range(len(barcodeList)):
                for n2 in range(n1+1,len(barcodeList)):
                    H = sum( [int(letter1 != letter2) for (letter1, letter2) in zip(barcodeList[n1],barcodeList[n2])] )             
                    if H in distribution:
                        distribution[H] += 1
                    else:
                        distribution[H] = 0
            
            hammingDistanceDistributions.append(distribution)
        
        self.hammingDistanceDistributions = hammingDistanceDistributions      
        print("INFO: Finished Calculating Hamming Distance Distributions.")        
        return    

    def output(self):
        return self.sequenceVariants
        
if __name__ == "__main__":

    
    MAXVARIANTS = 12472
    
    import csv
    
    uniqueSequences = []
    variants = []
    variantCounter = 0
    
    counter = 0
    NUM_NEGATIVES = 30
    while True:
        #Negative Controls
        seq = "".join([random.choice(['A','G','C','T']) for n in range(36)])
        if 'GACGTC' in seq or 'TCTAGA' in seq or 'GAATTC' in seq or 'GGATCC' in seq:
            pass
        else:
            variants.append( {'name' : 'NegControl', 'sequence' : seq})
            uniqueSequences.append(seq)
            variantCounter += 1
            counter += 1
            if counter >= NUM_NEGATIVES: break
                                  
    for (n,filename) in enumerate(filenameList):
        variants.append( {'name' : 'PspCas13b_WT', 'sequence' : 'GTTGTGGAAGGTCCAGTTTTGAGGGGCTATTACAAC' })
        variants.append( {'name' : 'RanCas13b_WT', 'sequence' : 'GTTGGGACTGCTCTCACTTTGAAGGGTATTCACAAC' })
        uniqueSequences.append('GTTGTGGAAGGTCCAGTTTTGAGGGGCTATTACAAC')
        uniqueSequences.append('GTTGGGACTGCTCTCACTTTGAAGGGTATTCACAAC')
        variantCounter += 2
        with open(pathToFiles + filename, 'r') as csvfile:
            reader = csv.reader(csvfile)
            header = reader.__next__()
            sequenceConstraint = header[1].strip()
            structureConstraint = header[3].strip()
            Lmax = int(header[5])
            for row in reader:
                id = row[0].strip()
                seq = row[1].strip().replace('U','T')
                if seq in uniqueSequences: 
                    pass
                else:
                    str_revcomp = revcomp(seq)
                    if 'GACGTC' in str_revcomp or 'TCTAGA' in str_revcomp or 'GAATTC' in str_revcomp or 'GGATCC' in str_revcomp:
                        print('Found undesired restriction site in variant %s' % id)
                        pass
                    else:
                        variants.append( {'name' : 'C%s-#%s' % (str(n+1),id),
                                          'sequence' : str_revcomp,
                                          'sequenceConstraint' : sequenceConstraint,
                                          'structureConstraint' : structureConstraint,
                                          'Lmax' : Lmax
                                         })
                        variantCounter += 1
                        uniqueSequences.append(seq)
                if variantCounter >= MAXVARIANTS: break
    
            
    oligopool = {'design' : [ {'type' : 'restriction_site',     'name' : 'AatII',                   'sequence' : 'GACGTC'},
                              {'type' : 'primer_binding_site',  'name' : '5p cDNA primer',          'sequenceConstraint' : 'NNNNNNNNNNNNNNNCCM', 'Lmax' : 15, 'direction' : 'forward'},
                              {'type' : 'barcode',              'name' : 'barcode',                 'sequenceConstraint' : 'NNNNNNNNNNNNM', 'Lmax' : 9},
                              {'type' : 'restriction_site',     'name' : 'XbaI',                    'sequence' : 'TCTAGA'},
                              {'type' : 'constant',             'name' : 'linker',                  'sequence' : 'ATATAT'},
                              {'type' : 'restriction_site',     'name' : 'EcoRI',                   'sequence' : 'GAATTC'},
                              {'type' : 'constant',             'name' : 'terminator',              'sequence' : 'GCACCGAGTCGGTGCTTTTTTT'},
                              {'type' : 'variable',             'name' : 'test',                    'variants' : variants},
                              {'type' : 'constant',             'name' : 'guide RNA',               'sequence' : 'GGAAGACGGTGCTGTGGAAGGTGAAATCAA'},
                              {'type' : 'restriction_site',     'name' : 'BamHI',                   'sequence' : 'GGATCC'}
                            ],
                'excludeSequences' : ['AAAA','TTTT','CCCC','GGGG'],
                'Tm_target' : [63.0, 65.0],
                'Tm_settings' : {'Na' : 1.0, 'K' : 50.0, 'Mg' : 2.0, 'dNTPs' : 0.20, 'Tris' : 25.0, 'saltcorr' : 7},
                'maxOligoLength' : 170,
                'backgroundSequenceList' : ['CAATGTTGACGCTGACGCTTGAGAAGAGAAAAGAAAACCGCCGATCCTGTCCACCGCATTACTGCAAGGTAGTGGACAAGACCGGCGGTCTTAAGTTTTTTGGCTGAATCATTACTTCATGATGGCGTATTCCCCAAAAGCCTTTTTAATTGACATCGCAATTTCGGGTAACGCCTTAATTTCCACAACGCCCTTGTCCGGATAATTGTTATGGTCGAATGCATTTCGAATTTTTCTCAGTATGTCGCTTTGCTCCTTATTAATATTCTTATTATTCAGGAGGATTTTGAGAATAGATTTAAAATCAACTTTTTCCTCACGGTCAACACGCGCGCTCAGTTCAGGATACGTATCGAAGGCCCATTTTTCCAGGTTGAACACGATGCTGCTTATTTCTGGGCGACATTGATCATATTTATTAAATTCTTCCATTATATCCTCTTTGCTAACGATATCACTCCCTACCAGTTCCAGTAAATTGCCAATGCGTTTGTCCGAGGCCAGTACGAAGAAATCTCCGTAGTTTTTCAATTTCATGCCCTCCGAGGTAATCGTGTATTTTTTGCCGCCTTTTTCAAAAGTGAATGACATGGGCATGATCTCTGACAAGATACCTTTCTCTGCATCCGGCATAATTTCCTTTAGCTTGAAGCGTTCGCCATCGAAATCCGCCAGCTCCGTCAACGTTTTTTTGGCCAAGAGGAATAATAAAGCATCCTGTACCCGATAGCGTCGTATTACCTTCTCGCTTTTCTGATATTCATTTCGTGAGTTAGACAGACGTTTATCCAGGATCGTCTCAATTTCTTCGGATGAGGCATTGCGCATCTGTCTATTCGACCGTATCTTGTTGGATGCCTGCTTGCGGTAACGCTCGGTGCGCGACGCCCTTTCTTTCCATAATCCTTCCCGCTCTTCTACAGAAGTAAAACAATGCTGTAAGGACCCTTTCCGATCGTATTCCCCTTTGAGCATATCCATGTACCGATAATTACGATTCCATTGATAAAACGTTTGAAAATCGTCGTCAAGTACGCGCTTCATGTATTCGGCAATCAAATAGGTGACATTCGCATTGTTGAAATCAATTCCCTCCATCTGAGGTAAAGATTTCAGATGGCTTTTAATCTCATTGTCAAACATCTGTCGCGGCAATTCGACAGGAAGATCTTCAGAATATATGCGGCCCAGCGTTTTCATGGCGGGTGTTTTCCATTTGTTTTGATCTCTACGGATAAAGGGAACATCAACGCGATTACCCTTCTTTATCTCATTGGATAGTCCGGTGAGATAAAATTTGCGTTCAATTAGATAGCGCTCATAAAATTCCACTGCATTTGCCGGTATTGAACGGGCGAAGACCTTATACAAAAATGGATGCGGTTCAGTTGTCCCTTTACCTATGAGTCTTGCCTTTTCGAACATCAGCTTGAACTGTTGTTTGGCCTCATAATCGTCACCCGAATCATAAACAGCGATAGCACTCTGCATGATCCGGTAATTTAAACCCGTAATCTTGTTCTCGCCGTCATTCACCGACGGCTGAAAAAGAACGATATCTTTCGCTAAGAAATCGGCAAGTTTACCGGTGCTAATCTGCTTAAAACCACGTTTCCCCATCTTGTTATCTGCGCTACGTATTGACTTCCTGTCGTCTTTGAAGCGCTTAATCCTTCGTTCGGTATCTGTAAGCATATCATCCACGGTTAAGCGTATAAAGGCATCAACATCTTTGCCGTGAGCGTTACCACTGATCAGGTCCAGTATCTTCTGCGGTAAATCGCTCTCAGCAATGCCGAACGAGGCAATATTCTCGGCGGTGACCTCCTCTTTCTGCATAGCCTGAAATAACCTCTTGTAGCGGTTGTGGACGTCAACGATTAATTTCTCCGTTTTTTTTGAGCCGAAGAGGAACATATGAAAGGCCATTGCGGGAATTTCCAGCGTGGACATTCTGCACGAAGGAATCGTTTTAACAACATACCGATCGTCTTCAATTACCGGCAGCAACGGTGCACTATCTTCTTTATCATTAATAAACATCTCAACTTTATTATTCTCCAGTATATAATGTGTATAAGTATCAACAATGTAGGGGTAATTAGCAGGGTTCGCGTCATCTCGTTTCATGTTCTCAAAGTCTCTAATGCGTATACCGCTGTTGCCAAAGGTGCCGTTTTCTTGCTTACGCATCGTCTCCGCTTCTTCCAGCCGTCCGAAACCGTTAAGGGGCTGTTCTATCACTCGGACTCGCGTCTGTCCGTCGATGCAAGTCTTGTCAGCCTTTAAAAGATACCGCAATTTGCCCATATTCACATGAAATCGGATATGATCGAACAGTTTTCCGTAATCAATGTATTGCAGCAGAAGGGGCACAAAGCGGTCGCTGCTTCTTTTCATAAGGACCTCATTATGATCATCGCTGATTATTCTGAAGCGCGACTGCTTCTCAGCACTCAACGTAGTGAACAGTTCATCGGGACACCGTTTAACTTCATTCAGCATATCCATGGCAACGCTTTTATTAGACTTTTCGCTATGGATTCGATCCTTCGGAAGTTTAATAGAATTGATGCCAAAGCTCCGTATGATAATTCGCCTTTCCTCGCTTTGCGCGTTGTAGGAACTGAAAATTGGAAGGCGAGACAGAAAAATATTAATATATTGTTTATCCAGGAACAGGCAAATTAGGAGAGCAATTCCGACACCAGACAGATGAAGCTTCTTTTGGGTATCGCCATTATAATCCTGAAGGCTTAAAAAGAACCCCGTATTGACTTGACTTTTCTTTTTTCCATAGGCATCCTTGACAAATTTAAAGCGTTTGTCTTGAATGAACGCCAGGTCCTCTGTTTTGTAGCCATACCTCTCATTCATGTTGCGCAAAGCAACAGTATAATAATTATTTATCATACCACTCAGAGGTTGCTCCGTACTCGTCAGAAATTCACAACCGTCATTCAGTTTCTCCTCGTACGTTTTATAGTGATTAGTGAGATCCCTATACATTTTGAGCACACCGAATGCTCGCTTCAGGACTTCAAAGATATCGTTTGAGTTAACTTCCACTCTGTTCTGCTTATACTTGCCATTAGAGTATTCACGCTGATTCTCCGCCATAATCTTCAGAAAAGGAAAATATGACTGCAGGCGCTCTATTATAAACATTGTTTTTTCCGGCTGCTTGTCATATCCGTTCTTTGCGTTGTACAAGTGTGACATAACAGGATGAAACCACAGGTTCTCGTTGTTTTCATTTTGCTCGCCTTCTATATCAGCCACTTTCTGGATATGATCCAGCACCGTCTGTGCGTTCAACATGGCCATGACGCTATAAGTACCGAAGTATTTTTTCTGGTTTTCTACCAAAGCCGGTATATTCATCATATCACTTCCTTGACATGGTGTTTTAGACCTTAACGATACGGTACGTTTCGTATCATGTCAATTGGTAACGAATCAGATTCCACCGTACGTCGTTTTCGCCAGATATCCTCGAGCGCGGAATTCAAAAAAAGCACCGACTCGGTGCGTTGTAATAGCCCCTCAAAACTGGACCTTCCACAACTTGATTTCACCTTCCACAGCACCGTCTTCCGGATCCATTATACTCGAGGCAATGATGTGTGTCAAGACTATGTCTTTAGCATCGGGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGAATAGGGAGGTCGTGATTTGACGTTTTGTTGACCTCTACGCTACAAGGGAGACCAAAAAAAAAAAAACACCCGTTAGGGTGTTCAATAATTGGATTTTAGCCCAGTCCAACCTAAGCCTTTGACACCCGGAAAAGCTCACTTTATAATGACGTCNNNNNNNNNNNNNTCTAGAACCCGCCATATACGGCGGGACATCTGCGGGGGTTCCATATGAAACAGAACAAAGAACAGATCACCAAACAGAACATCAACGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAAGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAAGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGCTGCCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAAGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAATAAACATACATGATCCAGCCCATTCGTGGGCTGGATCTTTTTTTTTCAGGTTTGTCCATTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGCCTATGGACTATGTTTGAAAGGGAGAAATACTAGATGGCACGTACCCCGAGCCGTAGCAGCATTGGTAGCCTGCGTAGTCCGCATACCCATAAAGCAATTCTGACCAGCACCATTGAAATCCTGAAAGAATGTGGTTATAGCGGTCTGAGCATTGAAAGCGTGGCACGTCGCGCCGGTGCAGGCAAACCGACCATTTATCGTTGGTGGACCAACAAAGCAGCACTGATTGCCGAAGTGTATGAAAATGAAATCGAACAGGTACGTAAATTTCCGGATTTGGGTAGCTTTAAAGCCGATCTGGATTTTCTGCTGCATAATCTGTGGAAAGTTTGGCGTGAAACCATTTGTGGTGAAGCATTTCGTTGTGTTATTGCAGAAGCACAGTTGGACCCTGTAACCCTGACCCAACTGAAAGATCAGTTTATGGAACGTCGTCGTGAGATACCGAAAAAACTGGTTGAAGATGCCATTAGCAATGGTGAACTGCCGAAAGATATCAATCGTGAACTGCTGCTGGATATGATTTTTGGTTTTTGTTGGTATCGCCTGCTGACCGAACAGTTGACCGTTGAACAGGATATTGAAGAATTTACCTTCCTGCTGATTAATGGTGTTTGTCCGGGTACACAGTGTTGATAAGCCCCCGGAAGATCATTCCGGGGGCTTTTTTATTGCGGCCGCTCGCTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAAGCCCAATCTGAATAATGTTACAACCAATTAACCAATTCTGATTAGAAAAACTCATCGAGCATCAAATGAAACTGCAATTTATTCATATCAGGATTATCAATACCATATTTTTGAAAAAGCCGTTTCTGTAATGAAGGAGAAAACTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGGTATCGGTCTGCGATTCCGACTCGTCCAACATCAATACAACCTATTAATTTCCCCTCGTCAAAAATAAGGTTATCAAGTGAGAAATCACCATGAGTGACGACTGAATCCGGTGAGAATGGCAAAAGTTTATGCATTTCTTTCCAGACTTGTTCAACAGGCCAGCCATTACGCTCGTCATCAAAATCACTCGCATCAACCAAACCGTTATTCATTCGTGATTGCGCCTGAGCGAGACGAAATACGCGATCGCTGTTAAAAGGACAATTACAAACAGGAATCGAATGCAACCGGCGCAGGAACACTGCCAGCGCATCAACAATATTTTCACCTGAATCAGGATATTCTTCTAATACCTGGAATGCTGTTTTTCCGGGGATCGCAGTGGTGAGTAACCATGCATCATCAGGAGTACGGATAAAATGCTTGATGGTCGGAAGAGGCATAAATTCCGTCAGCCAGTTTAGTCTGACCATCTCATCTGTAACATCATTGGCAACGCTACCTTTGCCATGTTTCAGAAACAACTCTGGCGCATCGGGCTTCCCATACAAGCGATAGATTGTCGCACCTGATTGCCCGACATTATCGCGAGCCCATTTATACCCATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTCGACGTTTCCCGTTGAATATGGCTCATAACACCCCTTGTATTACTGTTTATGTAAGCAGACAGTTTTATTGTTCATGATGATATATTTTTATCTTGTGCAATGTAACATCAGAGATTTTGAGACACGGGCCAGAGCTGCACCCCAGTCTTAGTTCGCC'],
                'enablePadding' : True
                }

    outputFilename = 'Oligopool_8_7_2020'
    job = OligoPoolCalculator(oligopool, pathToFiles + outputFilename)
    job.run()
    output = job.output()

    
    
