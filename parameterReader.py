import argparse

class ParameterReader:

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-f","--file",required=False,dest="file")
        self.parser.add_argument("-s","--smiles",required=False,nargs="+",dest="smiles")
        self.parser.add_argument("-o","--output",required=False,dest="output")

    def readParameters(self):
        arguments = self.parser.parse_args()
        return arguments


