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
    except requests.exception.RequestException as e:
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
    except requests.exception.RequestException as e:
        print(f"Error requesting {uniprotID}: Error {e}")
        exit()

    for elem in variant_json["features"]:
        if "association" in elem:
            disease_vars.append(elem["association"])

    return disease_vars

