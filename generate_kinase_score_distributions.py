# to generate distributions of 89,784 known sites for percentile ranking (Johnson et al., 2023).
# each kinase will have its own .txt file named [kinase_name].txt in the folder kinase_score_distributions
# it will be a json file with a list containing all raw scores for these peptides ( divided by (1/20)**10 )

import sys, os, pprint, regex, statistics, math, PSSM as pssm, re, json
current_filepath = os.path.dirname(os.path.realpath(__file__))

pssm_data = pssm.PSSM()

#get a list of all human kinases from the study:
all_kinases = pssm_data.all_kinases

#rank the site.

'''

json_file = open(current_filepath + '/kinase_score_distributions/AAK1.json', 'r')

distribution = json.load(json_file)
score = pssm_data.probability('MSFTSPRKSP', 'AAK1').final_score
distribution.append(score)

distribution.sort()
ind = distribution.index(score)

print(ind)


print(len(distribution))

print(ind/(len(distribution))*100)

json_file.close()

'''

#sys.exit()

#get a list of all phosphosites (and surrounding sequences) detected in the study (and make them 10 AAs):

#make a kinase file

k = 1 #kinase number
for kinase in all_kinases:
    
    #skip the first 30 kinases, we already have it in there.
    if k > 200:
    
        sequence_file = open(current_filepath + '/kinase_score_distribution_sequences.txt', 'r')
        write_file = open(current_filepath + '/kinase_score_distributions/' + kinase + '.json', 'w')
        list_of_scores = []

        x = 0 #line number
        for line in sequence_file.readlines():

            if x != 0:
                #skip the header.

                line = line.strip()
                sequence = line[2:12] #only extract 10 AAs since our PSSMs are 10 positions.

                #replace _ characters with X characters (since my program uses X for unknown sequences).
                sequence = re.sub('_', 'X', sequence)

                list_of_scores.append(pssm_data.probability(sequence, kinase).final_score)
                
            int_check = x/1000
            int_check = isinstance(int_check, int)

            if int_check is True:

                print("kinase: " + str(x))

            x+=1

        json_txt = json.dumps(list_of_scores)
        write_file.write(json_txt)   
        write_file.close()
        
        sequence_file.close()

    if k == 303:
        print("Done with " + str(k))
        break

    k += 1



