#export values of the PSSMs for each kinase.

import PSSM, os, time

unix = int(time.time()) #UNIX timestamp.
current_output_path = os.path.dirname(os.path.realpath(__file__)) + '/output'
output_file_name = 'all_kinase_values_' + str(unix) + '.csv'

pssm = PSSM.PSSM()

kinase_dict = {}

for kinase in pssm.all_kinases:
    
    kinase_dict[kinase] = {
    
        'min_score': pssm.probability('AAAAATAAAA', kinase).min_score*(10**9),
        'mean_score': pssm.probability('AAAAATAAAA', kinase).probability_phos_random_peptide*(10**9),
        'max_score':pssm.probability('AAAAATAAAA', kinase).max_score*(10**9)
        
    }
    

with open(current_output_path + '/' + output_file_name, 'w') as filehandle:
    
    filehandle.write('Kinase\tmin_score\tmean_score\tmax_score\n')
    
    for kinase, info in kinase_dict.items():            
        filehandle.write(f"{kinase}\t{info['min_score']}\t{info['mean_score']}\t{info['max_score']}\n")
    