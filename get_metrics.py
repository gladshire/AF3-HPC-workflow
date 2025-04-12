'''
(Dirty) Code for extracting score metrics from AlphaFold3 runs

Written by Miles Woodcock-Girard for Drew Lab at UIC
'''

import json
import csv
import os
import sys
import glob
import subprocess

import pdockq as pDockQ
import uniprot_api
#import pae
import utils


with open("metrics.csv", 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ",")
    csv_writer.writerow(["Pair", "ipTM", "pTM", "Ranking Score", "pDockQ", "ipSAE", "Disease Variant", "Novel Interaction"])

    for i, pair in enumerate(os.listdir()):
        if os.path.isdir(pair) == False:
            continue

        print(f"{i} {pair}")

        if utils.usedAF3Server(pair):

            if not glob.glob(f"./{pair}/fold_{pair}_*_summary_confidences_*.json"):
                print(f"Empty directory: {pair} Skipping ...")
                continue

            best_model = utils.getBestModel(pair)

            if not glob.glob(f"./{pair}/fold_{pair}_*_summary_confidences_{best_model}.json"):
                print(f"Empty directory: {pair} Skipping ...")
                continue

            curr_json_scores = glob.glob(f"./{pair}/fold_{pair}_*_full_data_{best_model}.json")[0]
            curr_json_summary = glob.glob(f"./{pair}/fold_{pair}_*_summary_confidences_{best_model}.json")[0]
            curr_mmcif_path = glob.glob(f"./{pair}/fold_{pair}_*_model_{best_model}.cif")[0]

        else:
            curr_json_scores = f"./{pair}/{pair}_confidences.json"
            curr_json_summary = f"./{pair}/{pair}_summary_confidences.json"
            curr_mmcif_path = f"./{pair}/{pair}_model.cif"

        if not os.path.exists(curr_json_summary) or not os.path.exists(curr_mmcif_path):
            continue


        prot1 = pair.split("_")[0]
        prot2 = pair.split("_")[1]

        with open(curr_json_summary) as curr_json_file:
            json_data = json.load(curr_json_file)
    
            # Obtain ipTM, pTM, Ranking Scores from confidences JSON
            iptm = json_data["iptm"]
            ptm = json_data["ptm"]
            ranking_score = json_data["ranking_score"]

            # Obtain chain coordinates, pLDDTs from mmCIF file
            chain_coords, chain_plddt = pDockQ.read_model_file(curr_mmcif_path)

            # Calculate pDockQ score for protein pair
            pdockq, ppv = pDockQ.calc_pdockq(chain_coords, chain_plddt, 8)

            # Calculate ipSAE score for protein pair
            ipsae_output = curr_mmcif_path[:-4:] + "_15_15.txt"
            subprocess.run(["python3", "ipsae.py", curr_json_scores, curr_mmcif_path, str(15), str(15)])
            with open(ipsae_output, 'r') as ipsae_file:
                ipsae_data = ipsae_file.read().splitlines()[1:-1:]
                ipsae_max = float(ipsae_data[3].split()[5])


            # Discern whether disease variants exist for protein
            numVars1 = uniprot_api.count_disease_var(prot1)
            numVars2 = uniprot_api.count_disease_var(prot2)

            if numVars1 > 0 or numVars2 > 0:
                disease_var = True
            else:
                disease_var = False

            # Determine if interaction can be found on UniProt
            known_interact = uniprot_api.interaction_in_uniprot(prot1, prot2)

            # Get average PAE across all pairs of residues
            #avg_pae, pae_data, nres_c1, nres_c2 = pae.load_json(curr_json_scores)

            # Write pair data to row in CSV output
            csv_writer.writerow([pair, iptm, ptm, ranking_score, pdockq, ipsae_max, disease_var, known_interact])

            curr_json_file.close()

    csv_file.close()
