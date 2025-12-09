#PSSM (position-specific probability matrix) based on Johnson et al., 2023
#this uses the PSSM ser_thr_all_norm_matrices.txt file which is the normalized values to each position in column.

import sys, os, pprint, regex, statistics, math, json

file_path = os.path.dirname(os.path.realpath(__file__))

class PSSM:
    
    def __init__(self, kinase = None):
        
        #parse the probability matrix.
        
        filehandle = open(file_path + '/pssm_ser_thr_all_norm_matrices.txt', 'r')
        lines = filehandle.readlines() #read line by line...
        self.matrix = {}
        
        self.positions = ( -5, -4, -3, -2, -1, 1, 2, 3, 4 )
        self.aas = ( 'P', 'G', 'A', 'C', 'S', 'T', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'D', 'E', 'pS', 'pT', 'pY' )
        
        #define kinase equivalents (human to yeast):
        
        self.mitotic_kinases = {
        
            'Mps1': 'TTK',
            'Ipl1': 'AURB',
            'Cdc5': 'PLK1',
            'Bub1': 'BUB1',
            'Cdk1': 'CDK1',
            'Alk1/Alk2': 'HASPIN',
            'Yck1/Yck2 (CK1A)': 'CK1A',
            'Yck1/Yck2 (CK1A2)': 'CK1A2',
            'Yck1/Yck2 (CK1G1)': 'CK1G1',
            'Yck1/Yck2 (CK1G2)': 'CK1G2',
            'Yck1/Yck2 (CK1G3)': 'CK1G3',
            'Yck1/Yck2 (CK1D)': 'CK1D',
            'Yck1/Yck2 (CK1E)': 'CK1E',
            'Cka1/Cka2 (CK2A1)': 'CK2A1',
            'Cka1/Cka2 (CK2A2)': 'CK2A2',
            'Cdc15 (MST1)': 'MST1',
            'Cdc15 (MST2)': 'MST2',
            'Dbf2/Dbf20 (NDR1)': 'NDR1',
            'Dbf2/Dbf20 (NDR2)': 'NDR2'
            #'Swe1': '' #THIS IS A TYROSINE KINASE AND THEREFORE WE NEED MORE INFO.
            
        }
        
        self.all_kinases = []
        
        #generate the probability matrix in a nested dictionary:
        
        count = 1
        
        for line in lines:
            
            line = line.strip()
            
            if line != '':
            
                line = line.split('\t')

                if count > 1:

                    self.matrix[line[0]] = {}
                    self.all_kinases.append(line[0])

                    col = 1 #start pulling matrix values out at column 1 (with index being 0), i.e. the second column

                    for position in self.positions:

                        self.matrix[line[0]][position] = {}

                        for aa in self.aas:

                            self.matrix[line[0]][position][aa] = float(line[col])

                            col += 1


                count += 1
            
    def probability(self, ten_aas, kinase):
        
        #get the probability at a position for a particular kinase
        
        #define variables
        self.max_score  = None
        self.min_score  = None
        self.score      = None
        
        if ten_aas is not None:

            self.kinase = kinase
            self.warning = True #a warning that will trigger when there is not enough sequence context (i.e. the S/T is near the beginning or end of a sequence). Default to True and then turn off if no X's
            max_scores = [] #a list that stores the max scores for each aa position
            min_scores = [] #a list that stores the min scores for each aa position.
            mean_scores = [] #a list that stores the mean scores for each aa position. (i.e. this is equal to the probability of phosphorylating a random peptide - if random peptides were sampled infinitely, the probabilities will average to this value)

            #we want to get the min and max scores for this kinase.
            #IMPORTANT: the max score will NEVER be reached if the user supplies any X's (i.e. the phosphosite is close to the N- or C-termini and is therefore missing sequence context.
            #this will skew the probability. therefore, calculate the min and max scores taking into account whether there are X's or not.

            #Let's see how many X's:

            positions_of_xs = []
            searches = regex.finditer('(X)', ten_aas)

            if searches is not None:

                for search in searches:

                    position = search.span()[1]-6 #-6 is relative to the 0 position being the phosphorylated residue.
                    positions_of_xs.append(position)

            if len(positions_of_xs) == 0:

                positions_of_xs = None
                self.warning = False #turn off the trigger warning.

            positions_to_calculate = list(self.positions) #it's a tuple, remember?
            #remove the positions that have X's:

            if positions_of_xs is not None:
                for x in positions_of_xs:
                    positions_to_calculate.remove(x)
                    
            each_residue = [] # extract each residue from ten_aas. remember they can be phos so one residue is S or pS for example.
            residue_search = regex.finditer('([p]?[A-Z])', ten_aas)
            
            if residue_search is not None:
                for residue in residue_search:
                    each_residue.append(residue.group())
                    
            for position in positions_to_calculate:

                position_scores = []

                for score in self.matrix[kinase][position].values():

                    position_scores.append(score)

                max_scores.append(max(position_scores))
                min_scores.append(min(position_scores))
                mean_scores.append(round(statistics.mean(position_scores), 6))
                
            self.max_score = float(1)
            self.min_score = float(1)
            self.probability_phos_random_peptide = float(1)
            
            #get the probability of phosphorylating a random peptide (i.e. the means of each position):
            for mean_score in mean_scores:
                self.probability_phos_random_peptide *= mean_score
                
            for max_sc in max_scores:
                self.max_score *= max_sc #the maximum probability you can obtain with this kinase.
                
            for min_sc in min_scores:
                self.min_score *= min_sc #the minimum probability you can obtain with this kinase.
                
            #self.max_score = self.max_score/self.probability_phos_random_peptide
            #self.min_score = self.min_score/self.probability_phos_random_peptide

            #self.max_score = sum(max_scores)
            #self.min_score = sum(min_scores)

            #calculate the score for this peptide:
            self.raw_score = 1
            for pos in positions_to_calculate:
                actual_pos = pos + 5 #position in the ten_aas string. the pos var is the position relative to the phosphoresidue.
                aa = each_residue[actual_pos]
                self.raw_score *= self.matrix[kinase][pos][aa]
                
            #the final score is how many more times likely is the given peptide going to be phosphorylated compared to a random peptide.
            #calculate an odds ratio (how likely is this kinase to phosphorylate this site over how likely is this kinase likely to phosphorylate a random peptide (converged at the mean, that is... the mean!)
            self.odds_ratio = self.raw_score/self.probability_phos_random_peptide
            self.odds_ratio = round(self.odds_ratio, 2)
            
            self.final_score = self.raw_score/((1/20)**10)
            self.final_score = round(self.final_score, 2)
            
 
        return self

    def rank_kinases_with_known_sites(self, ten_aas, mitotic = True):
        
        #ranks all kinases based on where the raw score sits within the distribution of 85,000 reported phospho-ST residues (percentile rank, as defined by Johnson paper)
        
        if ten_aas is not None:
            
            score_per_kinase = {}
            
            for human_kinase in self.all_kinases:
                #get the score for each kinase:
                
                score_per_kinase[human_kinase] = self.probability(ten_aas, human_kinase).final_score_supp_note_2
                
            #for human_kinase in self.all_kinases:
                #Put all known scores in a list for each kinase
            
        else: 
            return None
            
    def rank_kinases_odds_ratio(self, ten_aas, mitotic = True):
        
        #ranks all kinases based on the the peptide score for that kinase (raw score / probability of PHOSPHORYLATING a random peptide).
        # @mitotic - whether to only look at mitotic kinases or ALL kinases. Default to mitotic kinases if the user does not supply 'false' for second argument.
        
        if ten_aas is not None:
        
            ranks = {}

            if mitotic is True:
            
                for yeast_kinase, human_kinase in self.mitotic_kinases.items():
                    #get the odds-ratio for each kinase, store it in a dictionary.
                    ranks[yeast_kinase] = self.probability(ten_aas, human_kinase).odds_ratio
            else:
                
                for human_kinase in self.all_kinases:
                    #get the percentile score for each kinase, store it in a dictionary.
                    ranks[human_kinase] = self.probability(ten_aas, human_kinase).odds_ratio

            #sort the ranks by highest to lowest value in dictionary:
            rank = sorted(ranks.items(), key=lambda x: x[1], reverse=True)
            rank_sorted = {}

            for kinase_rank_pair in rank:

                    rank_sorted[kinase_rank_pair[0]] = kinase_rank_pair[1]

            return rank_sorted
            
        else:
            return None
        
    def rank_kinases_based_on_score_percentile(self, ten_aas, mitotic = True):
        
        #ranks all kinases based on the raw score from Supplementary Note 2 (Johnson et al, 2023 - Cantley lab)
        # @mitotic - whether to only look at mitotic kinases or ALL kinases. Default to mitotic kinases if the user does not supply 'false' for second argument.
        
        if ten_aas is not None:
            
            #get the rank of the sequence within the distribution of 89,000 human phosphorylation sites per kinase (stored in kinase_score_distributions) folder:
            #first, we need to extract the list of scores from each kinase file.
            
            #a dictionary containing each kinase and the percentile of the site within each kinase distribution:
            ranks = {}
            
            if mitotic is True:
                
                for yeast_kinase, human_kinase in self.mitotic_kinases.items():
                    
                    kinase_file = open(file_path + '/kinase_score_distributions/' + human_kinase + '.json', 'r')

                    distribution = json.load(kinase_file)
                    kinase_file.close()

                    score = self.probability(ten_aas, human_kinase).final_score
                    distribution.append(score) #add the final score to the distribution for that kinase

                    distribution.sort() #re-sort the scores from lowest to highest.
                    indx = distribution.index(score) #where in the list (position) does our final_score sit within the distribution?

                    #get the percentile of this score against all other scores:
                    percentile = (indx/len(distribution))*100
                    percentile = round(percentile, 2)

                    #add the kinase to the list as key, and the site percentile as value.
                    ranks[yeast_kinase] = percentile

                    del distribution #delete this crazy long list to free up memory.

                rank = sorted(ranks.items(), key=lambda x: x[1], reverse=True)
                rank_sorted = {}

                for kinase_rank_pair in rank:

                        rank_sorted[kinase_rank_pair[0]] = kinase_rank_pair[1]
                    
                
            else:
            
                for human_kinase in self.all_kinases:
                    kinase_file = open(file_path + '/kinase_score_distributions/' + human_kinase + '.json', 'r')

                    distribution = json.load(kinase_file)
                    kinase_file.close()

                    score = self.probability(ten_aas, human_kinase).final_score
                    distribution.append(score) #add the final score to the distribution for that kinase

                    distribution.sort() #re-sort the scores from lowest to highest.
                    indx = distribution.index(score) #where in the list (position) does our final_score sit within the distribution?

                    #get the percentile of this score against all other scores:
                    percentile = (indx/len(distribution))*100
                    percentile = round(percentile, 2)

                    #add the kinase to the list as key, and the site percentile as value.
                    ranks[human_kinase] = percentile

                    del distribution #delete this crazy long list to free up memory.

                rank = sorted(ranks.items(), key=lambda x: x[1], reverse=True)
                rank_sorted = {}

                for kinase_rank_pair in rank:

                        rank_sorted[kinase_rank_pair[0]] = kinase_rank_pair[1]

            return rank_sorted
        
        else:
            return None
        
        '''
        if ten_aas is not None:
        
            ranks = {}

            if mitotic is True:
            
                for yeast_kinase, human_kinase in self.mitotic_kinases.items():
                    #get the percentile score for each kinase, store it in a dictionary.
                    ranks[yeast_kinase] = self.probability(ten_aas, human_kinase).final_score_supp_note_2
            else:
                
                for human_kinase in self.all_kinases:
                    #get the percentile score for each kinase, store it in a dictionary.
                    ranks[human_kinase] = self.probability(ten_aas, human_kinase).final_score_supp_note_2

            #sort the ranks by highest to lowest value in dictionary:
            rank = sorted(ranks.items(), key=lambda x: x[1], reverse=True)
            rank_sorted = {}

            for kinase_rank_pair in rank:

                    rank_sorted[kinase_rank_pair[0]] = kinase_rank_pair[1]

            return rank_sorted
            
        else:
            return None
            
        '''
        