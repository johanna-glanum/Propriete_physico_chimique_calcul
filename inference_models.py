from joblib import load
import pandas as pd
import numpy as np
import pickle


MODELS = ["Dummy","KNN" ,"SVM" , "Tree", "Random Forest", "MLP" ]
dico_AA = pd.read_csv('./data/pka_AA.csv',sep=';', encoding='latin-1', decimal = ",")
dico_AA = dico_AA.copy()
dico_AA.index = dico_AA["Sym"]
dico_AA = dico_AA[["pKa1_c","pKa2_n","pKr","pHi"]]
all_AA = dico_AA.index.to_numpy()
N_AA = len(all_AA)


def preprocess_sequence(sequence):
    
    features = np.zeros(N_AA)
    for i in range(N_AA):
        features[i] = sequence.count(all_AA[i])
    return features



#loads models
    #solubility models

solubity_models = {}
pI_models = {}
for model_name in MODELS:
    with open("./models/solubility/{}_solubility.joblib".format(model_name), "rb") as f:
        solubity_models[model_name] = pickle.load(f)
    with open("./models/pI/{}_pI.joblib".format(model_name), "rb") as f:
        pI_models[model_name] = pickle.load(f)




def predict_solubility(sequence):
    result = {}
    features = preprocess_sequence(sequence)
    for model in solubity_models.keys():
        prediction = solubity_models[model].predict(features)
        if prediction == 0:
            result[model] = "Non soluble"
        else:
            result[model] = "Soluble"
    
    return pd.DataFrame(result)


def predict_pI(sequence):
    result = {}
    features = preprocess_sequence(sequence)
    for model in solubity_models.keys():
        prediction = pI_models[model].predict(features)
        result[model] = prediction
    
    return pd.DataFrame(result)





