'''
Code for extracting score metrics from AlphaFold3 runs
Written by Miles Woodcock-Girard for Drew Lab at UIC
'''

import json
import csv
import os
import sys

import pdockq as pDockQ
import vars_uniprot



with open("metrics.csv", 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ",")
    csv_writer.writerow(["Pair", "ipTM", "pTM", "Ranking Score", "pDockQ", "Disease Variant"])

    for pair in os.listdir():
        if os.path.isdir(pair) == False:
            continue

        curr_json_path = f"./{pair}/{pair}_summary_confidences.json"
        curr_mmcif_path = f"./{pair}/{pair}_model.cif"

        if not os.path.exists(curr_json_path) or not os.path.exists(curr_mmcif_path):
            continue

        prot1 = pair.split("_")[0]
        prot2 = pair.split("_")[1]

        with open(curr_json_path) as curr_json_file:
            json_data = json.load(curr_json_file)
    
            # Obtain ipTM, pTM, Ranking Scores from confidences JSON
            iptm = json_data["iptm"]
            ptm = json_data["ptm"]
            ranking_score = json_data["ranking_score"]

            # Obtain chain coordinates, pLDDTs from mmCIF file
            chain_coords, chain_plddt = pDockQ.read_cif(curr_mmcif_path)
            pdockq, ppv = pDockQ.calc_pdockq(chain_coords, chain_plddt, 8)

            # Obtain number of disease variants
            numVars1 = vars_uniprot.count_disease_var(prot1)
            numVars2 = vars_uniprot.count_disease_var(prot2)

            if numVars1 > 0 or numVars2 > 0:
                disease_var = True
            else:
                disease_var = False

            csv_writer.writerow([pair, iptm, ptm, ranking_score, pdockq, disease_var])

            curr_json_file.close()

    csv_file.close()
