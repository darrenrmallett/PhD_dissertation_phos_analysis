import os, sys, regex, re, mass_spec, proteins as p, xlsxwriter, time, kinetochore as KT, sgd_phosphorylation as sgd, json
from xlsxwriter.utility import xl_col_to_name, xl_rowcol_to_cell

current_output_path = os.path.dirname(os.path.realpath(__file__)) + '/output'

workbook = xlsxwriter.Workbook(current_output_path + '/kinetochore_predicted_sites'), {'constant_memory': True}) #constant_memory = writes row-by-row to dump memory after writing

# # # third worksheet - contains predicted phosphorylation sites for the whole kinetochore # # # 

worksheet3 = workbook.add_worksheet('All pred. kin. sites') #all kinetochore proteins
worksheet3.freeze_panes(1, 0)

#set column widths:
worksheet3.set_column(0, 0, 15) #column A = 15
worksheet3.set_column(1, 1, 13) #column B
worksheet3.set_column(2, 2, 10) #column C
worksheet3.set_column(3, 3, 12) #column D
worksheet3.set_column(4, 4, 18)                                   
worksheet3.set_column(5, 5, 13)
worksheet3.set_column(6, 6, 28)
worksheet3.set_column(7, 7, 33)
worksheet3.set_column(8, 8, 12)
worksheet3.set_column(9, 9, 100)

#conditional formatting:
worksheet3.conditional_format(f'E1:E1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

    'type': 'cell',
    'criteria': 'equal to',
    'value': '"NO"',
    'format': red_format
})

worksheet3.conditional_format(f'E1:E1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

    'type': 'cell',
    'criteria': 'equal to',
    'value': '"YES"',
    'format': green_format
})

worksheet3.conditional_format(f'I1:I1048576', { #entire column = H1:H1048576 (the last possible row is 1048576).

    'type': 'cell',
    'criteria': 'equal to',
    'value': '"Reported"',
    'format': blue_format
})

#write the header string to the Excel worksheet:

row = 0
col = 0

output_header = [

    'Complex/type',
    'Accession',
    'Protein',
    'Potential phosphosite',
    'Detected?',
    'Predicted kinase',
    'Kinase consensus (regex)',
    'Sequence +/- 5',
    'Literature report? (SGD)',
    'Citation'

]

header_format = workbook.add_format()
header_format.set_align('center')
header_format.set_align('vcenter')
header_format.set_text_wrap()

for header_item in output_header:

    worksheet3.write(row, col, header_item, header_format)

    col += 1 #move to next column every header item

row = 1
col = 0

current_protein = '' #the name of the protein we are currently looping through. Important for the collapsible row feature in Excel.
new_protein = False
l_blue = True #formatting for the collapsable rows.

# To read the JSON file containing all published phosphorylation sites already known:    
json_file = open('published_phos.json', 'r')
read = json_file.read()
json_dict = json.loads(read) #the contents of the json file, converted to a python dictionary.
json_file.close()

for desc, complx_dict in kinetochore.proteins.items():
    for complx, proteins in complx_dict.items():

        for protein in proteins:                
            all_predicted_phos_sites = {} #a dictionary with site (number) as key, and another dictionary containing info as values
            all_predicted_phos_sites_sorted_dict = {} #the same dictionary as before, except sorted by residue order.

            accession = kinetochore.accessions[protein] #use the accession, since there are discrepancies with our common KT protein names and SGD common name
            protein = p.Protein(accession).name #now convert the protein identified from the accession back to the common name.

            if current_protein != protein:
                new_protein = True #we are going to insert a collapsible row in this iteration of the loop.
            else:
                new_protein = False

            current_protein = protein

            #loop through each kinase consensus (regex):
            for kinase, consensus in p.Protein(protein).consensus.items():


                find_sites = p.Protein(protein).find_sites(consensus)
                sites = find_sites.phosphoresidues

                if len(sites) >= 1:
                    #there are matched consensus sequences!

                    for site in sites:

                        #put this site in the all_predicted_phos_sites dictionary
                        #remember a site can be predicted to be phos by more than one kinase!

                        if site not in all_predicted_phos_sites.keys():

                            all_predicted_phos_sites[site] = {

                                'aa': p.Protein(protein).sequence[site-1],
                                'kinase': [kinase],
                                'consensus': [consensus]

                            }

                        else:

                            #the site is already in the dictionary. must match more than one kinase consensus...

                            all_predicted_phos_sites[site]['kinase'].append(kinase)
                            all_predicted_phos_sites[site]['consensus'].append(consensus)

            # now we have looped through all the kinases.
            # the phosphosites are going to be out of order, since we searched by kinase sequence.
            # generate a list of phosphosites, then sort them

            all_predicted_phos_sites_list = list(all_predicted_phos_sites.keys())
            all_predicted_phos_sites_list.sort()


            for site in all_predicted_phos_sites_list:

                all_predicted_phos_sites_sorted_dict[site] = all_predicted_phos_sites[site] #now they are in order.

            formatting = {

                0: workbook.add_format({'align': 'center'}),
                1: workbook.add_format({'align': 'center'}),
                2: workbook.add_format({'align': 'center'}),
                3: workbook.add_format({'align': 'center'}),
                4: workbook.add_format({'align': 'center'}),
                5: '',
                6: workbook.add_format({'font_name': 'Courier'}),
                7: workbook.add_format({'align': 'center', 'font_name': 'Courier'}),
                8: workbook.add_format({'align': 'center'}),
                9: ''
            }

            #create the collapsible row:

            worksheet3.outline_settings(True, False, True, False)

            if new_protein is True:

                #worksheet3.set_row(row, None, None, {'level': 1, 'collapsed': True})

                sum_of_predicted_sites = len(all_predicted_phos_sites_list)

                #etract sum of sites reported:
                #some sites are repeatedly reported. therefore put them as dictionary keys to get the unique number:
                unique_reported_sites = {}
                for resideaux in all_predicted_phos_sites_list:
                    if str(resideaux) in json_dict['proteins'][accession].keys():
                        unique_reported_sites[resideaux] = 0                                           

                sum_of_sites_reported = len(unique_reported_sites.keys())

                collapse_output_item = [
                    complx,
                    accession,
                    protein,
                    sum_of_predicted_sites,
                    '',
                    '',
                    '',
                    '',
                    sum_of_sites_reported,
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '',
                    '' #added extra columns to make the blue a little longer.
                ]

                cc = 0 #column
                for item in collapse_output_item:

                    #header rows formatting

                    light_blue = workbook.add_format({'bg_color': '#E1EFFF', 'align': 'center', 'border': 1, 'border_color': '#A5C9F5'})
                    darker_blue = workbook.add_format({'bg_color': '#D5E8FF', 'align': 'center', 'border': 1, 'border_color': '#A5C9F5'})

                    if l_blue is True:
                        worksheet3.write(row, cc, item, light_blue)
                    else:
                        worksheet3.write(row, cc, item, darker_blue)

                    cc += 1

                if l_blue is True:
                    l_blue = False
                else:
                    l_blue = True

                row += 1 #increase the row.

            #now loop through the sorted list:

            for site, site_info in all_predicted_phos_sites_sorted_dict.items():

                # was the site detected to be phosphorylated in any of the MS files?                            
                detected = 'NO' #default

                if protein in all_phosphosites.keys():
                    for phosphosite in all_phosphosites[protein].keys():

                        if phosphosite == site:
                            detected = 'YES'

                aa_site = site_info['aa'] + str(site)

                predicted_kinases = ''
                ii = 1
                for kinases in site_info['kinase']:
                    predicted_kinases += kinases

                    if ii < len(site_info['kinase']):
                        predicted_kinases += '; '

                kinase_consensuses = ''
                ii = 1
                for consensus in site_info['consensus']:
                    kinase_consensuses += consensus

                    if ii < len(site_info['consensus']):
                        kinase_consensuses += '; '

                #surrounding sequence:                    
                site_minus_1 = site-1
                five_minus = (site-1)-5
                five_plus = (site-1)+6

                if five_minus < 0:
                    five_minus = 0 #we don't want to go below zero, the starting methionine.

                if five_plus > (p.Protein(protein).length-1):
                    five_plus = p.Protein(protein).length-1

                surrounding_seq = p.Protein(protein).sequence[five_minus:site] + '*' + p.Protein(protein).sequence[site:five_plus]
                surrounding_seq = '[' + str(five_minus+1) + ']-' + surrounding_seq + '-[' + str(five_plus) + ']'

                #literature report? :

                literature_report = ''
                citation = ''
                PMID = ''

                for acc, reported_phos_info in json_dict['proteins'].items():

                    if accession == acc:

                        #check whether this site has been reported before:
                        if str(site) in reported_phos_info.keys():

                            literature_report = 'Reported'

                            ii = 0                                
                            for ref in reported_phos_info[str(site)]['refs']:
                                PMID = reported_phos_info[str(site)]['ref_PMIDs'][ii]
                                PMID = str(PMID)

                                citation += ref + ' (PMID: ' + PMID + ')'

                                if ii < len(reported_phos_info[str(site)]['refs'])-1:
                                    citation += '; '
                                ii+=1

                output_items = [

                    '',
                    '',
                    protein,
                    aa_site,
                    detected,
                    predicted_kinases,
                    kinase_consensuses,
                    surrounding_seq,
                    literature_report,
                    citation

                ]


                for output_item in output_items:                        

                    worksheet3.set_row(row, None, None, {'level': 1, 'hidden': True})
                    worksheet3.write(row, col, output_item, formatting[col])

                    col += 1

                row += 1
                col = 0 #reset back to 0 for each site.

                worksheet3.write(row, col, 'Note that this tab only contains sites that have been predicted to be phosphorylated by a kinase provided. This is not a full list of all reported phosphorylation.')
                
# # # fourth worksheet (the regex info):

worksheet4 = workbook.add_worksheet('Regex used') #all kinetochore proteins

formatting = {

    0: workbook.add_format({'align': 'center'}),
    1: workbook.add_format({'align': 'left', 'font_name': 'Courier'})

}

#set column widths:
worksheet4.set_column(0, 0, 20) #column A = 15
worksheet4.set_column(1, 1, 60) #column B

worksheet4.write(0, 0, 'Kinase')
worksheet4.write(0, 1, 'Regular expression used')

row = 1

for kinase, reg_exp in p.Protein('NDC80').consensus.items(): #it just needs a protein name in order to work... I just put Ndc80.

    worksheet4.write(row, 0, kinase, formatting[0])
    worksheet4.write(row, 1, reg_exp, formatting[1])

    row += 1


workbook.close()