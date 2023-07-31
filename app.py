
import streamlit as st
import streamlit_echarts
from streamlit_echarts import st_echarts
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
    st.text(sequence)
    st.text(proportions)
else:
    st.text('En attente d\'une séquence')


proportions = plot_amino_acid_proportions(sequence)

options = {
    "xAxis": {
        "type": "category",
        "data": amino_acids,
    },
    "yAxis": {"type": "value"},
    "series": [
        {"data": proportions}
    ],
}

st_echarts(options=options, height="500px")
