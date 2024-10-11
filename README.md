# KoaLA

KoaLA is an AI tool for the prediction of the octanol-air partition coeficient.

KoaLA requires the use of Anaconda to install its dependencies. To install them run the following command:

```
conda env create -f environment.yml
```

To activate the environment created use the following command:

```
conda activate koala
```

To use KoaLA you can execute the main.py file with the following options:

- -f: To indicate the csv file with the SMILES. This file must have a column named smiles.
- -s: To indicate the SMILES you want to use. Following the option type all the SMILES separated by a space.
- -o: To indicate the name of the output file. If not indicated, the default value is results.csv.

Example:

```
python3 main.py -f molecules.csv -s C COC C#N -o output.csv
```

