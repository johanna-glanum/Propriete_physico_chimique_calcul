
import streamlit as st
import matplotlib
import matplotlib.pyplot as plt
from collections import Counter

def plot_amino_acid_proportions(sequence):
    # Calcul des fréquences des acides aminés dans la séquence
    amino_acid_counts = Counter(sequence)
    total_amino_acids = len(sequence)
    amino_acids = list(amino_acid_counts.keys())
    proportions = [count / total_amino_acids * 100 for count in amino_acid_counts.values()]
    return proportions, amino_acids


sequence = st.text_input('Write a sequence')
proportions, amino_acids = plot_amino_acid_proportions(sequence)

if st.button('Calculer la proportion des AA'):
    st.text(amino_acids)
    st.text(proportions)
    
    chart_data = pd.DataFrame(
    proportions,
    columns=amino_acids)

    st.bar_chart(chart_data)
else:
    st.text('En attente d\'une séquence')

proportions = plot_amino_acid_proportions(sequence)


import pandas as pd
import numpy as np



