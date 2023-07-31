
import streamlit as st

import matplotlib.pyplot as plt
from collections import Counter

def plot_amino_acid_proportions(sequence):
    # Calcul des fréquences des acides aminés dans la séquence
    amino_acid_counts = Counter(sequence)
    total_amino_acids = len(sequence)
    amino_acids = list(amino_acid_counts.keys())
    proportions = [count / total_amino_acids * 100 for count in amino_acid_counts.values()]

    # Création du graphique "plot bar"
    plt.bar(amino_acids, proportions)
    plt.xlabel('Acides aminés')
    plt.ylabel('Pourcentage')
    plt.title('Proportion des acides aminés dans la séquence peptidique')
    plt.grid(False)
    plt.show()

sequence = st.text_input('Write a sequence')
st.text(plot_amino_acid_proportions(sequence))

