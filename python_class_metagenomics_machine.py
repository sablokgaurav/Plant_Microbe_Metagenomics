# a machine learning class for metagenomics from the sequence to multidimensional scaling 
# this class takes a metagenomics sequence and expression and also the other variables 
# associated with the sequences and then prepares the data for the classifier 
# and trains the variables and the expression using the classifier and prepares the 
# model for the expressed model. Applicable to the meta-transcriptomics with defined variables
import os 
import pandas as pd
import numpy as np 
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
class multiDimensionalMetagenomics:
    def __init__(self, expression, fastafile):
        self.expression = expression
        self.fastafile = fastafile 
    def read_fasta(self): 
        self.fasta_read = []
        self.fasta_names = ""
        with open("self.fastafile", "r") as fasta:
            for i in fasta.readlines():
                if i.startswith(">"):
                    self.fasta_read.append(self.fasta_names)
                    self.fasta_names = i.strip()
                else:
                    self.fasta_read += i.strip()
        names = list(filter(None,[i[0] for i in 
                                   ([i.split("\t") for i in self.fasta_read])]))
        sequences = [i[1] for i in (list(filter(lambda n: n!=[''],
                                    [i.split("\t") for i in self.fasta_read])))]
        self.fastadf = pd.DataFrame([(i,j) for i,j in 
                                        zip(self.fasta_read, self.fasta_names)])
    def getSequenceDataframe(self):
        return self.fastadf
    def matrixNormalization(self):
        self.expdata = pd.read_csv("self.expression", sep = ",")
        self.motifscan = []
        while True:
            self.takemotif = input("Enter the name of the motifs to scan:")
            self.motifscan.append(self.takemotif)
            if self.takemotif == "None":
                break
            return self.motifscan
        self.motifname = []
        while True:
            self.name = input("Enter the name of the motifs:")
            self.motifname.append(self.name)
            if self.motifname == "None":
                break 
            return self.motifname 
        self.motifdata = pd.DataFrame([(i,j) for i,j in 
                                zip(self.motifscan, self.motifname)])
        if len(self.expdata) != len(self.motifdata):
            raise ValueError("motifdata and expression data dont match")
        print(f"motifdata and expression data dont match")
        if len(self.expdata) == len(self.motifdata):
            self.motifexp = pd.concat([self.motifdata, self.expdata])
        return self.motifexp
    def getMotifdataframe(self):
            return self.motifexp
        
    def makeSelection(self, selected_datasets):
            if selected_datasets:
                selected_datasets = []
                selected_datasets.append([f"selected_{i}" for i in self.motifname])
                for i in range(len(self.motifscan)):
                    selected_datasets[i] = self.fastadf["sequence"].str.contrains(self.motifscan[i])
                    self.cdmotif = pd.concat(selected_datasets)
                    return self.cdmotif      
    def  motifReadcount(self, count):
        self.count = count
        if self.expression:
            with open("self.expression", "r") as expr:
                self.expcount = pd.read_csv(self.expression, sep = ",")
                self.expcountcat = pd.concat([self.cdmotif,self.expcount])
                self.expcountf = self.expcountcat["count"].apply(lambda n: n >= self.count)
            
    def machineLabels(self):
        self.machinetags = []
        while True:
            taketags = input("Please enter the machine tags:")
            self.machinetags.append(taketags)
            if taketags == "None":
                break
            return self.machinetags
        self.machinesequences = []
        while True:
            takemachine = input("Please enter the machine classifications:")
            self.machinesequences.append(takemachine)
            if takemachine == "None":
                break
            return self.machinesequences
        if len(self.machinesequences) != (self.machinesequences):
            raise ValueError("the number of the machine tags and the number of the machine seqeunces differ")
        print("Please enter the same number of the machine tags and sequences")
        self.cdemtags = self.expcountcat["sequence"].apply(lambda n: self.machinetags[i] if 
                                                self.machinesequences[i] in n else "None"
                                  for i in range(len(self.machinetags)))
        return self.cdemtags
    
    def sequenceXGBOOSTPrep(self, params, columnsD):
        while True:
            takeparams = input("Please enter the path to the params files:")
            self.paramsdf = pd.read_csv(takeparams, sep = ",")
            if takeparams is None:
                break
            return self.paramsdf
        if len(self.paramsdf) != len(self.cdemtags):
            raise ValueError("given parameters cant be attached to the fasta and the expression datasets")
        print("given parameters cant be attached to the fasta and the expression datasets")
        self.cdemtagsparams = pd.concat([self.cdemtags,self.paramsdf])
        self.columnsD = []
        while True:
            takecolumns = input("Please enter the columns to drop before performing XGBOOST:")
            self.columnsD.append(takecolumns)
            if takecolumns == "":
                break
            return columnsD
        self.cdemtagsparams.drop(columns = set(self.columnsD))
    
    def xGBOOSTtrain(self, filename,columns=None,):
        self.filename = filename
        X = self.cdemtagsparams.drop(columns = set(self.columnsD))
        if columns is None:
            y = self.cdemtagsparams
        else:
            y = self.cdemtagsparams["column"]
        self.model = XGBClassifier(booster='gbtree',objective='binary:logistic', random_state=2)
        self.scores = cross_val_score(self.model, X, y, cv=3)
        print('Accuracy:', np.round(self.scores, 2))
        print('Accuracy mean: %0.2f' % (self.scores.mean()))
        with open(self.filename, "r") as f:
            f.write(self.model)
            f.write('Accuracy:', np.round(self.scores, 2))
            f.write('Accuracy mean: %0.2f' % (self.scores.mean()))
