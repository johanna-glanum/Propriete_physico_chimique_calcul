from joblib import load
import pandas as pd
import numpy as np
import pickle
import sklearn
from conjoint_triad_sequence import conjoint_triad


dico_AA = pd.read_csv('./data/pka_AA.csv',sep=';', encoding='latin-1', decimal = ",")
dico_AA = dico_AA.copy()
dico_AA.index = dico_AA["Sym"]
dico_AA = dico_AA[["pKa1_c","pKa2_n","pKr","pHi"]]
all_AA = dico_AA.index.to_numpy()
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
        return "Non soluble"
    else:
       return "Soluble"
    


def predict_pI(sequence):
    features = preprocess_sequence_frequence(sequence)
    return pI_model.predict(features)[0]


def predict_solubility_regression(sequence):
    features = preprocess_sequence_conjoint(sequence)
    return solubity_model_regression.predict(features)[0]
    





