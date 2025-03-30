# Code for screening proteins for disease causing variants, queries Rest API with UniProt IDs
# Written by Miles Woodcock-Girard for Drew Lab at UIC


import requests
import json
import os


# Given a UniProt ID, count number of pathogentic variants associated with it
def count_disease_var(uniprotID):
    num_variants = 0

    variant_api_url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprotID}?format=json"
    try:
        response = requests.get(variant_api_url)
        response.raise_for_status()
        variant_json = response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error requesting {uniprotID}: Error {e}")
        exit()

    for elem in variant_json["features"]:
        if "association" in elem:
            num_variants += 1

    return num_variants


# Given UniProt ID, return list of associated pathogentic variants
def get_disease_vars(uniprotID):
    disease_vars = []

    variant_api_url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprotID}?format=json"
    try:
        response = requests.get(variant_api_url)
        response.raise_for_status()
        variant_json = response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error requesting {uniprotID}: Error {e}")
        exit()

    for feat in variant_json["features"]:
        if "association" in feat:
            disease_vars.append(feat["association"])

    return disease_vars


# Given UniProt ID, return list of its interactions
def interaction_in_uniprot(uniprotID_ref, uniprotID_qry):
    api_url_ref = f"https://rest.uniprot.org/uniprotkb/{uniprotID_ref}.json"
    api_url_qry = f"https://rest.uniprot.org/uniprotkb/{uniprotID_qry}.json"
    try:
        response_ref = requests.get(api_url_ref)
        response_ref.raise_for_status()
        response_ref_json = response_ref.json()

        response_qry = requests.get(api_url_qry)
        response_qry.raise_for_status()
        response_qry_json = response_qry.json()
    except requests.exceptions.RequestException as e:
        print(f"Error requesting {uniprotID_ref}: Error {e}")
        exit()

    # Get gene name
    gene_id_ref = response_ref_json["genes"][0]["geneName"]["value"]
    gene_id_qry = response_qry_json["genes"][0]["geneName"]["value"]

    api_url_ref = f"https://rest.uniprot.org/uniprotkb/{uniprotID_ref}.txt"
    try:
        response = requests.get(api_url_ref)
        response.raise_for_status()
        response_txt = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error requesting {uniprotID_ref}: Error {e}")
        exit()


    if gene_id_qry in response_txt:
        return False
    else:
        return True


