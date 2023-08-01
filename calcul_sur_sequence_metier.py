import numpy as np
import pandas as pd
import unicodedata
from time import time
import re

# from tabulate import tabulate
import time
import pickle
import string
import base64

# librairie affichage
import matplotlib.pyplot as plt
import csv
import re


dico_AA = pd.read_csv('./data/pka_AA.csv', sep=';',
                      encoding='latin-1', decimal=",")
dico_IONI = pd.read_csv('./data/IONISABLE.csv', sep=';',
                        encoding='latin-1', decimal=",")

dico_AA = dico_AA.copy()
dico_AA.index = dico_AA["Sym"]
dico_AA = dico_AA[["pKa1_c", "pKa2_n", "pKr", "pHi"]]
dico_AA.T.to_dict()

dico_IONI = dico_IONI.copy()
dico_IONI.index = dico_IONI["Sym"]
dico_IONI = dico_IONI[["protonnee", "deprotonnee"]]
dico_IONI.T.to_dict()

data_AA = pd.read_csv("./data/AA.csv", sep=";", encoding="latin-1")
data_AA.PM = data_AA.PM.str.replace(",", ".").astype(float)

# Définition des dicitonnaire de la base de connaissance pour AA

d_ioni = {'NH2+': 'NH', 'COOHc': 'cCOO-', 'SH': 'S-', 'COOH': 'COO-',
          'NH+': 'N', 'NH3+': 'NH2', 'OH': 'O-', 'NH3+n': 'nNH2'}

d_charge = {'NH2': 0, 'NH': 0, 'COOH': 0, 'COO-': -1, 'SH': 0, 'S-': -1, 'NH+': 1, 'N': 0,
            'NH3+': 1, 'NH2+': 1, 'OH': 0, 'O-': -1, 'NH3+n': 1, 'nNH2': 0, 'COOHc': 0, 'cCOO-': -1}


def calculer_phi_AA(dico_AA, seq):
    pka_n = dico_AA.get("pKa2_n")[seq]
    pka_c = dico_AA.get("pKa1_c")[seq]
    phi = (pka_n + pka_c) / 2

    if phi == dico_AA.get("pHi")[seq]:
        print("le calcul de phi est correct, il vaut : " + str(phi))
    else:
        print("Erreur de calcul pour l'acide aminé " + seq)


def calcul_pKa(seq, dico_AA, dico_IONI):

    dico_temp = {}

    for aa in seq:
        if aa in ['R', 'D', 'C', 'E', 'H', 'K', 'Y']:
            tmp = dico_IONI['protonnee'][aa]
            # print(tmp)

            if tmp not in dico_temp.keys():
                i = 1
                dico_temp[tmp] = dico_AA.loc[aa, 'pKr']
            else:
                dico_temp[tmp + "m" + str(i)] = dico_AA.loc[aa, 'pKr']
                i += 1

    dico_temp['NH3+n'] = dico_AA.get("pKa2_n")[seq[0]]
    dico_temp['COOHc'] = dico_AA.get("pKa1_c")[seq[-1]]

    # tri par ordre croissant
    marklist = sorted(dico_temp.items(), key=lambda x: x[1])
    dico_final = dict(marklist)

    return dico_final


def create_dataframe(dico_final, d_ioni, d_charge):
    taille = len(dico_final) + 1
    protonée = np.array([[key]*taille for key in dico_final.keys()])

    déproto = np.array([[d_ioni[re.sub(r"m\d+\Z", "", key)]]
                       * taille for key in dico_final.keys()])

    proto_charge = np.array(
        [[d_charge[re.sub(r"m\d+\Z", "", key)]]*taille for key in dico_final.keys()])
    deproto_charge = np.array(
        [[d_charge[d_ioni[re.sub(r"m\d+\Z", "", key)]]]*taille for key in dico_final.keys()])

    arr = np.char.add(np.tril(protonée), np.triu(déproto, k=1))
    df_visu = pd.DataFrame(arr)
    # print(protonée)
    print(" ")
    # print(déproto)
    print(" ")
    # print(df_visu)
    return df_visu


def compute_charge_and_pH(seq, df_visu, dico_final, d_ioni, d_charge, taille, pH=7, seuil_bas=-1, seuil_haut=1):

    proto_charge = np.array(
        [[d_charge[re.sub(r"m\d+\Z", "", key)]]*taille for key in dico_final.keys()])
    # print(proto_charge)
    deproto_charge = np.array(
        [[d_charge[d_ioni[re.sub(r"m\d+\Z", "", key)]]]*taille for key in dico_final.keys()])
    # print(deproto_charge)
    charge_table = np.tril(proto_charge) + np.triu(deproto_charge, k=1)
    # print(deproto_charge)
    charge_global = charge_table.sum(axis=0)
    # print(deproto_charge)
    indice_zero = np.where(charge_global == 0)[0][0]
    # print(deproto_charge)
    indice_avant = indice_zero - 1
    # print(deproto_charge)

    print(type(df_visu))

    valeur_pka = list(dico_final.values())[indice_avant]
    # print(deproto_charge)
    valeur_pkb = list(dico_final.values())[indice_zero]
    # print(deproto_charge)

    pHi_seq = ((valeur_pka + valeur_pkb)/2)
    # print(deproto_charge)

    zwiterion_0 = df_visu.iloc[:, indice_zero]
    # print(deproto_charge)
    zwiterion = zwiterion_0.values
    # print(deproto_charge)

    indice_bas_soluble = np.where(charge_global <= seuil_bas)[0][0]
    """
    print(charge_global)
    print(np.where(charge_global <= seuil_bas))
    print(np.where(charge_global >= seuil_haut))
    print(list(dico_final.values()))
    """
    indice_haut_soluble = np.where(charge_global >= seuil_haut)[0][-1]

    indice_pH = (np.abs(np.array(list(dico_final.values())) - pH)).argmin()

    charge_pH = charge_global[indice_pH]

    forme_sequence = " ".join(df_visu[indice_pH].to_numpy())

    soluble = True
    borne_inf_pH = list(dico_final.values())[0]
    borne_inf_soluble_pH = list(dico_final.values())[indice_haut_soluble]
    borne_sup_soluble_pH = list(dico_final.values())[indice_bas_soluble]
    borne_sup_pH = list(dico_final.values())[-1]

    print("A pH = {}, la charge global est {}".format(pH, charge_pH))

    print("La forme de la séquence est {}".format(forme_sequence))

    print(f"Le phi de la sequence {seq} est de :{pHi_seq}")
    print(" ")
    print(
        f"La forme des groupe ionisé du zwiterion de la séquence {seq} est :\n{zwiterion}")
    print(" ")
    print(
        f"Les deux pka qui encadre la forme neutre de la sequence {seq} est : {valeur_pka} et {valeur_pkb}")
    # print(charge_global)
    # print(indice_bas_soluble, indice_haut_soluble)

    if pH in [borne_inf_soluble_pH, borne_sup_soluble_pH]:
        print("la mollécule n'est pas soluble à pH = {}".format(pH))
        soluble = False
    print("la mollécule est soluble à pH = {}".format(pH))

    print("La sequence est soluble pour un pH dans l'intervalle [{}, {}] et [{}, {}]".format(
        borne_inf_pH, borne_inf_soluble_pH, borne_sup_soluble_pH, borne_sup_pH))
    print("Pas soluble dans l'intervalle [{}, {}]".format(
        borne_inf_soluble_pH, borne_sup_soluble_pH))
    # print("Indices des seuils : bas {}, haut {}".format(indice_bas_soluble, indice_haut_soluble))

    return pHi_seq, zwiterion, valeur_pka, valeur_pkb, forme_sequence, charge_pH, soluble, borne_inf_soluble_pH, borne_sup_soluble_pH, borne_inf_pH, borne_sup_pH


def verifier_seq_2(seq):
    for caractere in seq:  
        if caractere not in ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]:
            print("False : ", caractere)
            # print("erreur, la séquence ne contient pas uniquement des acides aminés protéinogènes")
            return False
    # print("La séquence contient uniquement des acides aminés protéinogènes")
    return True


def get_PM(chain):
    n = len(chain)
    return sum([data_AA[data_AA["Abr_L"] == x]["PM"].values[0] for x in chain]) - (n-1)*18.01528


def calcul_all(seq, pH=7):
    print(seq)
    if verifier_seq_2(seq):
        print("ok")
        dico_final = calcul_pKa(seq, dico_AA, dico_IONI)
        taille = len(dico_final) + 1
        try:
            df_visu = create_dataframe(dico_final, d_ioni, d_charge)
            pHi_seq, zwiterion, valeur_pka, valeur_pkb, forme_sequence, charge_pH, soluble, borne_inf_soluble_pH, borne_sup_soluble_pH, borne_inf_pH, borne_sup_pH = compute_charge_and_pH(
            seq, df_visu, dico_final, d_ioni, d_charge, taille, pH)

            pM = get_PM(seq)
        

            return pHi_seq, zwiterion, valeur_pka, valeur_pkb, forme_sequence, charge_pH, soluble, borne_inf_soluble_pH, borne_sup_soluble_pH, borne_inf_pH, borne_sup_pH, pM
      
        except TypeError:
            print(seq)
            return np.nan

    else:
        return 0,0,0,0,0,0,0,0,0,0,0,0
