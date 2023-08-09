from joblib import load
import pandas as pd
import numpy as np
import pickle
import sklearn
from conjoint_triad_sequence import conjoint_triad


AA = pd.read_csv("./data/AA.csv", sep=";", encoding="latin-1")
all_AA = AA["Abr_L"].to_numpy()
N_AA = len(all_AA)

pI_ref = pd.read_csv("./data/data_test_compa_10kda.csv", sep=";")
pI_ref.columns = list(pI_ref.iloc[0])
pI_ref = pI_ref.drop(0, axis=0)
pI_ref = pI_ref[[" sequence", " Avg_pI"]]




def preprocess_sequence_frequence(sequence):
    
    features = np.zeros(N_AA)
    for i in range(N_AA):
        features[i] = sequence.count(all_AA[i])
    return features.reshape(1,-1)

def preprocess_sequence_conjoint(sequence):
    return np.array(list(conjoint_triad(sequence=sequence).values())).reshape(1,-1)


#loads models
    #solubility models
with open("./models/solubility/Tree_Solub_one.pkl", "rb") as model_file:
    solubity_model_regression = pickle.load(model_file)

with open("./models/solubility_classif/SVM_solubility_classif_conj.pkl", "rb") as model_file:
    solubility_model_classification = pickle.load(model_file)

with open("./models/pI/MLP_pI_fr.pkl", "rb") as model_file:
    pI_model = pickle.load(model_file)





def predict_solubility_classification(sequence):
    features = preprocess_sequence_conjoint(sequence)
    prediction = solubility_model_classification.predict(features)
    if prediction[0] == 0:
        return "La sequence n'est pas soluble en classification"
    else:
       return "La sequence est soluble en classification"
    


def predict_pI(sequence):
    features = preprocess_sequence_frequence(sequence)
    return pI_model.predict(features)[0]


def predict_solubility_regression(sequence):
    features = preprocess_sequence_frequence(sequence)
    return solubity_model_regression.predict(features)[0]
    





