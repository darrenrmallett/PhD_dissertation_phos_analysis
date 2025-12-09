import proteins, PSSM, regex, re, kinetochore, time, os

current_output_path = os.path.dirname(os.path.realpath(__file__)) + '/output'

n = 0

scores_above_1 = {}

pssm = PSSM.PSSM()
unix = int(time.time()) #UNIX timestamp.

kinase = 'TTK'

output_file_name = 'all_proteins_' + kinase + '_' + str(unix) + '.csv'

print('Running... please wait...')

for protein in proteins.Protein('STU1').master_list:
    
    protein = proteins.Protein(protein)
    
    name        = protein.name
    sequence    = protein.sequence
    accession   = protein.identifier
    length      = protein.length

    # see if there is a serine or threonine in this protein sequence:
    
    basic_search = re.search('([ST])+', sequence)
    
    if basic_search is not None:
        
        #search all S and T's
        
        match_object = regex.finditer('([ST])', sequence, overlapped=True)
        
        for single_match in match_object:
            
            position_of_ST = single_match.span()[1]
            residue = single_match.group()
            
            residue_and_pos = residue + str(position_of_ST)
            
            if (position_of_ST-5) >= 1 and (position_of_ST+4) <= length:
                
                #this is the peptide we are going to score:                
                peptide = protein.surrounding_seq(position_of_ST)
                
                prob            = pssm.probability(peptide, kinase).final_score
                #rank            = pssm.rank_kinases(peptide, False)
                #rank_list       = list(rank.keys())
                #rank_number     = rank_list.index(kinase)+1
                
                nn = 1
                #rank_list_str = ''
                '''
                for kin in rank_list:
                    
                    rank_list_str += '(' + str(nn) + ') ' + kin + '; '
                    
                    nn+=1'''
                
                if prob > 1:
                    
                    #it is more likely to be phosphorylated than a random peptide is. put it in the scores_above_1 dict.
                    
                    if protein.name in scores_above_1.keys():
                        
                        scores_above_1[protein.name][residue_and_pos] = {
                        
                            'pep': peptide,
                            'o_r': prob #odds ratio
                            #'rank_num': rank_number,
                            #'rank_list_str': rank_list_str
                            
                        }
                        
                    else:
                        
                        scores_above_1[protein.name] = { 
                        
                            residue_and_pos: {
                            
                                'pep': peptide,
                                'o_r': prob #odds ratio
                                #'rank_num': rank_number,
                                #'rank_list_str': rank_list_str
                                
                            }
                            
                        }
    
    n+=1
    
    '''if n == 1:
        break'''
        


with open(current_output_path + '/' + output_file_name, 'w') as filehandle:
    
    filehandle.write('Protein\tSite\tPeptide\tOdds_ratio\tRank_num\tRank_list\n')
    
    for protein_name, info in scores_above_1.items():
        for residue, res_info in info.items():
            
            filehandle.write(f"{protein_name}\t{residue}\t{res_info['pep']}\t{res_info['o_r']}\n")