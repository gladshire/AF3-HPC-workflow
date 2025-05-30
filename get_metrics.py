'''
Code for extracting score metrics from AlphaFold3 runs
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

import pandas as pd

dc2_csv = "dc2_af3.csv"
dc2_df = pd.read_csv(dc2_csv)

existing_pairs = set()

# Check if metrics.csv exists and read existing "Pair" values
csv_exists = os.path.exists("metrics.csv")
if csv_exists:
    with open("metrics.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            if row:  # Avoid blank lines
                existing_pairs.add(row[0])

with open("metrics.csv", "a") as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ",")

    if not csv_exists:
        csv_writer.writerow(["Pair", "Gene1", "Gene2", "ipTM", "pTM", "Ranking Score", "DC2 Score", "pDockQ", "ipSAE", "Pathogenic Variant (A)", "Pathogenic Variant (B)", "Path. Var. Publication (A)", "Path. Var. Publication (B)", "Novel Interaction"])

    for i, pair in enumerate(os.listdir()):
        if os.path.isdir(pair) == False:
            continue

        # Determine if pair has been processed
        if pair in existing_pairs:
            print(f"Skipping already processed pair: {pair}")
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

            # Get DirectContacts2 score
            dc2_row = dc2_df[dc2_df["fset"] == f"frozenset({{\'{prot1.upper()}\', \'{prot2.upper()}\'}})"]
            if dc2_row.empty:
                dc2_row = dc2_df[dc2_df["fset"] == f"frozenset({{\'{prot2.upper()}\', \'{prot1.upper()}\'}})"]

            dc2_score_str = dc2_row.iloc[0]["score"]
            dc2_score = float(dc2_score_str)
            

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


            # Get gene names from Uniprot
            geneID1 = uniprot_api.get_gene_id(prot1)
            geneID2 = uniprot_api.get_gene_id(prot2)

            # Discern whether disease variants exist for protein
            numVars1 = uniprot_api.count_disease_var(prot1)
            numVars2 = uniprot_api.count_disease_var(prot2)

            if numVars1 > 0:
                has_var1 = True
            else:
                has_var1 = False

            if numVars2 > 0:
                has_var2 = True
            else:
                has_var2 = False

            # Check if variants have associated publications
            has_pub1 = uniprot_api.has_publication(prot1)
            has_pub2 = uniprot_api.has_publication(prot2)

            # Determine if interaction can be found on UniProt
            known_interact = uniprot_api.interaction_in_uniprot(prot1, prot2)

            # Write pair data to row in CSV output
            csv_writer.writerow([pair, geneID1, geneID2, iptm, ptm, ranking_score, dc2_score, pdockq, ipsae_max, has_var1, has_var2, has_pub1, has_pub2, known_interact])


            # Update metrics in dc2_df\
            fset_key_f = f"frozenset({{\'{prot1.upper()}', '{prot2.upper()}\'}})"
            fset_key_r = f"frozenset({{\'{prot2.upper()}', '{prot1.upper()}\'}})"

            dc2_index = dc2_df[dc2_df["fset"] == fset_key_f].index
            if dc2_index.empty:
                dc2_index = dc2_df[dc2_df["fset"] == fset_key_r].index

            if not dc2_index.empty:
                idx = dc2_index[0]
                dc2_df.at[idx, "AF3_ipTM"] = iptm
                dc2_df.at[idx, "AF3_pTM"] = ptm
                dc2_df.at[idx, "pDockQ"] = pdockq
            else:
                print(f"Warning: could not find fset for {pair}")
            
            curr_json_file.close()


    csv_file.close()

# Dump final dataframe to CSV
dc2_df.to_csv(dc2_csv, index=False)
