
import streamlit as st
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from calcul_sur_sequence_metier import calcul_all
from inference_models import predict_pI, predict_solubility_classification, predict_solubility_regression

def plot_amino_acid_proportions(sequence):
    # Calcul des fréquences des acides aminés dans la séquence
    amino_acid_counts = Counter(sequence)
    total_amino_acids = len(sequence)
    amino_acids = list(amino_acid_counts.keys())
    proportions = [count / total_amino_acids * 100 for count in amino_acid_counts.values()]
    return proportions, amino_acids


sequence = str(st.text_input('Write a sequence'))
pH = st.text_input("Donner un pH")


pHi_seq = None
if st.button("Calculer les propriétées de la séquence"):
    state_calc = True
    pHi_seq, zwiterion, valeur_pka, valeur_pkb, forme_sequence, charge_pH, soluble, valeur_pkc,valeur_pkb, pM = calcul_all(seq = sequence, pH=float(pH))

    st.text("Le poids moléculaire est {}".format(pM))

    st.text("A pH = {}, la charge global est {}".format(pH, charge_pH))

    st.text("La forme de la séquence est {}".format(forme_sequence))

    st.text(f"Le phi de la sequence {sequence} est de :{pHi_seq}")
    st.text(" ")
    st.text(
        f"La forme des groupe ionisé du zwiterion de la séquence {sequence} est :\n{zwiterion}")
    st.text(" ")
    st.text(
        f"Les deux pka qui encadre la forme neutre de la sequence {sequence} est : {valeur_pka} et {valeur_pkb}")
    if not soluble:
        st.text("la mollécule n'est pas soluble à pH = {}".format(pH))
    else:
        st.text("la mollécule est soluble à pH = {}".format(pH))

    #st.text("La sequence est soluble pour un pH dans l'intervalle [{}, {}] et [{}, {}]".format(
        #valeur_pkc,valeur_pkb))
    
    #st.text("Pas soluble dans l'intervalle [{}, {}]".format(
        #valeur_pkc,valeur_pkb))
    st.text("Soluble si c'est inferieur à {} et superieur à {}".format(valeur_pkc, valeur_pkb))
    st.text("Peu Soluble si pH compris entre [{}, {}]".format(valeur_pka, valeur_pkb))
    

    proportions, amino_acids = plot_amino_acid_proportions(sequence)

    #st.text(amino_acids)
    #st.text(proportions)

    fig, ax = plt.subplots()
    ax.bar(amino_acids, proportions)
    plt.xlabel('Acides aminés')
    plt.ylabel('Pourcentage')
    plt.title('Proportion des acides aminés dans la séquence peptidique')
    plt.grid(False)
    plt.show()

    st.pyplot(fig)

else:
    st.text('En attente d\'un pH')


if st.button("Prédire par machine learning la solubilité à pH 7 et son pI"):

    prediction_solub_binaire = predict_solubility_classification(sequence)
    prediction_solub_regression = predict_solubility_regression(sequence)
    prediction_pI = predict_pI(sequence)

    if pHi_seq != None:
        prediction_pI["val_calc"] = pHi_seq

    st.text("Prediction de la solubilité en 2 méthodes: ")
    
    st.text("Regression Tree {}".format(prediction_solub_regression))
    st.text("Score des métriques de Regression Tree: \n Erreur absolue moyenne : 23.419 \n Erreur relative moyenne : 2.9499")
   
    st.text("Classification SVM {}".format(prediction_solub_binaire))
    st.text("Score des métriques SVM : \n recall = 0.60279 \n accuracy = 0.60245 \n f1 = 0.60243")
    
    
    st.text("Prédiction du Point Isoélétrique : ")

    st.text("Regression MLP {}".format(prediction_pI))
    st.text("Score des métriques de Regression MLP: \n Erreur absolue moyenne : 0.1072 \n Erreur relative moyenne : 0.0155")
    
    









