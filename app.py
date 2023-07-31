
import streamlit as st

from collections import Counter

def plot_amino_acid_proportions(sequence):
    # Calcul des fréquences des acides aminés dans la séquence
    amino_acid_counts = Counter(sequence)
    total_amino_acids = len(sequence)
    amino_acids = list(amino_acid_counts.keys())
    proportions = [count / total_amino_acids * 100 for count in amino_acid_counts.values()]
    return proportions


sequence = st.text_input('Write a sequence')
if st.button('Calculer la proportion des AA'):
    st.text(plot_amino_acid_proportions(sequence))

