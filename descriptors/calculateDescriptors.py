import pandas as pd
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from . import descriptorsRDKit as descRDKit
from . import descriptorsJazzy as descJazzy
from . import descriptorsConjugated as descConj
from openbabel import pybel
from rdkit.Chem import CanonSmiles
from parameterReader import ParameterReader

pandas2ri.activate()

from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
rpy2_logger.setLevel(logging.ERROR)


def calculateDescriptorsCDK(smiles):
    rcdk = rpackages.importr("rcdk")
    base = rpackages.importr("base")
    descCategories = rcdk.get_desc_categories()
    getDescNames = rcdk.get_desc_names

    descNames = base.unique(base.unlist(base.sapply(descCategories,getDescNames)))
    descNames = [name for name in descNames if name != "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor"]
    descNames = ro.StrVector(descNames)

    mols = rcdk.parse_smiles(ro.StrVector(smiles))
    descriptorsCDK = rcdk.eval_desc(mols,descNames)

    descriptorsCDK = ro.conversion.rpy2py(descriptorsCDK)
    descriptorsCDK.index = [i for i in range(len(descriptorsCDK.index))]

    return descriptorsCDK




def calculateDescriptorsObabel(smiles):
    descriptorsObabel = pd.DataFrame()
    i = 0
    for smile in smiles:
        mol = pybel.readstring("smi",smile)
        desc = mol.calcdesc()
        descriptorsObabel = pd.concat([descriptorsObabel,pd.DataFrame(desc,index=[i])])
        i+=1
    descriptorsObabel.drop(columns=["cansmi","cansmiNS","formula","title","InChI","InChIKey","smarts"],inplace=True)
    columns = [col + "_Obabel" for col in descriptorsObabel.columns]
    descriptorsObabel.columns = columns
    return descriptorsObabel

def calculateDescriptors(smiles):
    smilesCanon = []
    for smile in smiles:
        try:
            smilesCanon.append(CanonSmiles(smile))
        except:
            print(f"Failed to canonize {smile}")
    descriptorsCDK = calculateDescriptorsCDK(smilesCanon)
    descriptorsObabel = calculateDescriptorsObabel(smilesCanon)
    descriptorsRDKit = descRDKit.calculateDescriptorsRDKit(smilesCanon)
    descriptorsJazzy = descJazzy.calculateDescriptorsJazzy(smilesCanon)
    descriptorsConj = descConj.calculateConjugatedDescriptors(smilesCanon)
    descriptors = pd.concat([descriptorsJazzy,descriptorsCDK,descriptorsRDKit,descriptorsObabel,descriptorsConj],axis=1)
    descriptors.insert(0,"smiles",smilesCanon)
    return descriptors

def createDataset(data,fileName):
    data.sort_index(inplace=True)
    smiles = data["smiles"].to_numpy()
    objetivo = data[data.columns[1]].to_numpy()
    dataset = calculateDescriptors(smiles)
    dataset.insert(1,data.columns[1],objetivo)
    for index in range(2,len(data.columns)):
        dataset.insert(2,data.columns[index],data[data.columns[index]].to_numpy())
    dataset.to_csv("desc_"+fileName,index=False)

def main():
    reader = ParameterReader()
    parameters = reader.readParameters()
    fileName = parameters["data"]
    createDataset(pd.read_csv(fileName),fileName)

if __name__ == "__main__":
    main()
