'''
Code for extracting score metrics from AlphaFold3 runs
Written by Miles Woodcock-Girard for Drew Lab at UIC
'''

import json
import csv
import os
import sys

import pdockq_mmcif
import pdockq



with open("metrics.csv", 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ",")
    csv_writer.writerow(["Pair", "ipTM", "pTM", "Ranking Score", "pDockQ"])

    for pair in os.listdir():
        if os.path.isdir(pair) == False:
            continue

        curr_json_path = f"./{pair}/{pair}_summary_confidences_0.json"
        curr_mmcif_path = f"./{pair}/{pair}_model_0.cif"

        if not os.path.exists(curr_json_path) or not os.path.exists(curr_mmcif_path):
            continue

        with open(curr_json_path) as curr_json_file:
            json_data = json.load(curr_json_file)
    
            # Obtain ipTM, pTM, Ranking Scores from confidences JSON
            iptm = json_data["iptm"]
            ptm = json_data["ptm"]
            ranking_score = json_data["ranking_score"]

            # Obtain chain coordinates, pLDDTs from mmCIF file
            chain_coords, chain_plddt = pdockq_mmcif.read_cif(curr_mmcif_path)
            pdockq, ppv = pdockq_mmcif.calc_pdockq(chain_coords, chain_plddt, 8)
            pdockq_old, ppv = pdockq.calc_pdockq(chain_coords, chain_plddt, 8)


            print(f"{pair}:\n  {pdockq_old}\n  {pdockq}")

            csv_writer.writerow([pair, iptm, ptm, ranking_score, pdockq])

            curr_json_file.close()

    csv_file.close()
