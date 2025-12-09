import requests, kinetochore, pprint

filename = open('kt_accessions.py', 'w')

KT = kinetochore.Kinetochore()
KT_proteins = {}

for complx, info in KT.proteins.items():
    
    for sub_complx, list_of_proteins in info.items():
    
        for protein_name in list_of_proteins:
            
            print(protein_name)
    
            r = requests.get(f"https://www.yeastgenome.org/backend/locus/{protein_name}")

            output = r.json()
    
            KT_proteins[protein_name] = output['format_name']
    

    
filename.write(pprint.pformat(KT_proteins))    
filename.close()