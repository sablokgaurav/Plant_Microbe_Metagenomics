# a complete visualization class for metagenomics abundance
# final check 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
def __init__(self, file, columns, items, dataframe):
    """_summary_
    input the dataframe and the number of the columns
    and which you want to plot and also the path of the
    data frame.
    """
    self.file = input("Enter the path of the metagenomics classification:")
    self.dataframe = pd.read_csv("self.file", sep=",", header=True)
    self.items = int(input("Enter the number of columns:"))
    self.columns = []
    while True:
        generate = input("Enter the columns:")
        self.columns += generate
        if len(self.columns) == self.items:
            break
        print(f'the entered columns are: {self.columns}')
        with open(dataframe, "r+") as data:
            data.write(self.dataframe[self.columns])
            data.close()
def summaryStats(self, filename):
    """_summary_
    this function returns the summary of the 
    dataframe. if filename is given then it
    will write to the file else will print 
    the summary to the console.
    """
    self._summary = self.dataframe.describe
    self._information = self.dataframe.info()
    self._mean = self.dataframe.mean()
    self._median = self.dataframe.median()
    if filename:
        with open(filename, 'r+') as f:
            f.write(self._summary)
            f.write(self._information)
            f.write(self._mean)
            f.write(self._median)
    else:
        print(f'the summary, information, mean, median are:
              {self._summary},
              {self._information},
              {self._mean},
              {self._median}') 
def getPivot(self, col, length, index, write):
    """_summary_
    making a pivot table for the metagenomics
    classification according to the user
    defined columns
    """
    self._length = length
    self._col = []
    while len(self._col) != self._length:
        take = input("Enter the cols for the pivot table:")
        col += take
        if len(self._col) == self._length:
            break
        self.pivot = pd.pivot_table(self.dataframe, 
                                    columns = self._col, 
                                    index=self.index)
        if write:
            with open(self.filename, 'r+') as f:
                f.write(self.pivot)
                f.close()
        else:
            print(f'the pivot table is: {self.pivot}')
def StackPlot(self, colours, loc, id1, id2):
    self._loc = loc
    self._id1 = id1
    self._id2 = id2
    """_summary_
    makes the abundance stack plot for
    the metagenomics dataset with the 
    defined columns. You can define
    your own colour schema. 
    """
    self._colours = []
    while True:
        colour = input("Enter the name of the colour for the palette:")
        self._colours += colour
        if colour is None:
            break
        print(f'the entered colours are:{self._colours}')
    stack_new = self.dataframe[self.columns]
    stack = stack_new(kind = 'barh', 
                       stacked = True,
                       color = self._colour)
    stack.legend(loc="self.loc", bbox_to_anchor = (self.id1, self.id2))

def getStackDistribution(self, column, bins, form):
    """_summary_
    make a stack distribution for the given
    columns for the metagenomics datasets
    """
    self._stackcolumn = column
    self._form = form
    self._bin = int(bins)
    return self.dataframe[self.].hist(normed = 0, 
                                    histype =self._form, 
                                    bins=self._bin)
    
def plotHighAbundance(self, abundance, col, filename):
    """_summary_
    Given the dataframe, plot the high abundance
    species after applying a filter according to
    the given criteria.
    """
    self._abundance = abundance
    self._col = col
    self.abundance_data = self.dataframe[self.dataframe["col"] > self._abundance]
    stack = self.abundance_data(kind = 'barh', 
                       stacked = True,
                       color = self._colour,)
    stack.legend(loc="self.loc", 
                 bbox_to_anchor = (self.id1, self.id2))
    return stack
    
def getPairPlot(self, columns):
    """_summary_
    plot all the correlations and
    the variables and produce a pairplot
    for the defined variables.
    """
    self._columns = []
    while True:
        getCol = input("Enter the name of the columns for the pair plot:")
        self._columns += getCol
        if getCol is None:
            break
    pairplot_data = self.dataframe[self._columns]
    print(f'the correlation between the variables:{new.corr()}')
    print(f'the variance between the variables:{new.var()}')
    return sns.pairplot(pairplot_data)


    
    

