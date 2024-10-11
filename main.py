import pickle
import os
import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from descriptors import calculateDescriptors
from parameterReader import ParameterReader

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

def reliability_message(similarity):
  if similarity is None:
    return "Could not calculate similarity, prediction realiabilty uncertain"
  elif similarity >= 0.8:
    return "Very reliable prediction"
  elif 0.8 > similarity >= 0.6:
    return "Reliable prediction"
  elif 0.6 > similarity >= 0.4:
    return "Poor prediction is likely"
  elif 0.4 > similarity:
    return "Poor prediction is very likely"

def calculate_fingerprint(smiles,fingerprint_generator):
  """Calculates the fingerprint of a list of molecules given their smiles representation.

  Args:
    smiles: A list of smiles .
    fingerpring_calculator: A function that calculates RDKit fingerprints.

  Returns:
    A list of the calculated fingerprints, one for each molecule.
  """
  fingerprints = []
  for s in smiles:
    try:
      mol = Chem.MolFromSmiles(s)
      fingerprint = fingerprint_generator.GetFingerprint(mol)
      fingerprints.append(fingerprint)
    except:
      fingerprints.append(None)
  return fingerprints

def mean_similarity_topN(input_smiles, comparison_smiles, fingerprint_generator,N=10):
  """Calculates the mean of the top 10 Tanimoto similarity values for the a given
  fingerprint of each input molecule in comparison with the train molecules.

  Args:
    input_smiles: A list of smiles to get the top 10 Tanimoto similarity values.
    comparison_smiles: A list of smiles to compare the input smiles with.
    fingerpring_generator: An RDKit Fingerprint generator.

  Returns:
    A list of mean top 10 Tanimoto similarity values for the given fingerprint, one for each input molecule.
  """
  results = []
  input_fps = calculate_fingerprint(input_smiles, fingerprint_generator)
  comparison_fps = calculate_fingerprint(comparison_smiles, fingerprint_generator)
  for input_fp in input_fps:
    if input_fp is not None:
      tanimoto_similarities = []
      for comparison_fp in comparison_fps:
        tanimoto_similarities.append(DataStructs.TanimotoSimilarity(input_fp, comparison_fp))
      top_10_values = sorted(tanimoto_similarities, reverse=True)[:N]
      mean_top_10 = np.mean(top_10_values)
      results.append(mean_top_10)
    else:
      results.append(None)
  return results

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
    #train_smiles = pd.read_csv("datasets/train.csv")["smiles"]
    #mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    #similarities = mean_similarity_topN(smiles,train_smiles,mfpgen,3)
    #reliability = [reliability_message(sim) for sim in similarities]
    #predictions = pd.DataFrame({"smiles":smiles,"Predicted logKoa":logKoa_prediction,"similarity":similarities,"Prediction reliability": reliability})
    predictions = pd.DataFrame({"smiles":smiles,"Predicted logKoa":logKoa_prediction})    
    return predictions

def evaluate_similarity(smiles, n_values,experimental_logKoa=None):
    train_smiles = pd.read_csv("datasets/train.csv")["smiles"]
    results = pd.DataFrame()
    results["smiles"] = smiles
    if experimental_logKoa is not None:
      descriptors = calculateDescriptors.calculateDescriptors(smiles)
      logKoa_prediction = predict(descriptors)
      error = np.abs(np.array(logKoa_prediction) - np.array(experimental_logKoa))
      results["experimental"] = experimental_logKoa
      results["predicted"] = logKoa_prediction
      results["error"] = error
    for n in n_values:
      mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
      similarities = mean_similarity_topN(smiles,train_smiles,mfpgen,n)
      results[f"{n}_similarity"] = similarities
    return results


def get_smiles_from_csv(file_path):
  smiles_colum_names = ["smiles","SMILES","Smiles"]
  data = pd.read_csv(file_path)
  for col_name in smiles_colum_names:
    if col_name in data.columns:
      return list(data[col_name])
  print(f"SMILES columns not found in {file_path}")
  return []


def main():
  parameters = ParameterReader().readParameters()
  predictions = []
  if parameters.smiles != None:
    smiles_predictions = make_predictions(parameters.smiles)
    predictions.append(smiles_predictions)
  if parameters.file != None:
    smiles = get_smiles_from_csv(parameters.file)
    file_predictions = make_predictions(smiles)
    predictions.append(file_predictions)
  if predictions:
    predictions = pd.concat(predictions,ignore_index=True)
  if parameters.output != None:
    predictions.to_csv(parameters.output)
  else:
    predictions.to_csv("predictions.csv")



if __name__ == "__main__":
  main()
