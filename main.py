import pickle
import os
import pandas as pd
import numpy as np
from descriptors import calculateDescriptors

def mean(list_values):
  total = sum(list_values)
  if len(list_values) > 0:
    return total/len(list_values)
  else:
    return np.nan

def load_model(model_name):
    model = pickle.load(open(f"models/{model_name}", 'rb'))
    return model

def predict(descriptors):
    model_names = os.listdir(f"models")
    predictions = {smiles:[] for smiles in descriptors["smiles"]}
    for model_name in model_names:
      if not model_name.startswith('.'):
        model = load_model(model_name)
        na_smiles = []
        non_na_smiles = []
        model_descriptors = list(model.feature_names_in_)
        prediction_descriptors = descriptors[["smiles"] + model_descriptors]
        if prediction_descriptors.isna().any().any():
          na_rows = prediction_descriptors[prediction_descriptors.isna().any(axis=1)]
          na_smiles = list(na_rows["smiles"])
          print(f"At least one value is missing for model {model_name} for molecules {na_smiles}")
          prediction_descriptors.dropna(axis=0,inplace=True)
        non_na_smiles = prediction_descriptors["smiles"]

        prediction_descriptors.drop(columns=["smiles"],inplace=True)
        prediction_descriptors.columns = prediction_descriptors.columns.astype(str)
        prediction_descriptors = prediction_descriptors.rename(str,axis="columns")
        model_prediction = model.predict(prediction_descriptors)
        if len(na_smiles) > 0:
          try:
            na_rows.columns = na_rows.columns.astype(str)
            model_prediction_na_smiles = model.predict(na_rows.drop(columns=["smiles"]))
            for smiles,prediction in zip(na_smiles,model_prediction_na_smiles):
              predictions[smiles].append(prediction)
          except:
            print(f"Prediction failed for {na_smiles}")
        for smiles,prediction in zip(non_na_smiles,model_prediction):
          predictions[smiles].append(prediction)
    return [mean(predictions[smiles]) for smiles in descriptors["smiles"]]

def make_predictions(smiles):
    if len(smiles) == 0:
      return
    unique_smiles = len(smiles) == 1
    if unique_smiles:
      smiles = smiles + ["C"]
    descriptors = calculateDescriptors.calculateDescriptors(smiles)
    if unique_smiles:
      descriptors.drop(descriptors.tail(1).index,inplace = True)
      smiles = smiles[:-1]
    logKoa_prediction = predict(descriptors)
    predictions = pd.DataFrame({"smiles":smiles,"Predicted logKoa":logKoa_prediction})
    return predictions

def main():
  smiles = input("Enter the SMILES:\n")
  print(make_predictions(smiles.split(" ")))

if __name__ == "__main__":
  main()
