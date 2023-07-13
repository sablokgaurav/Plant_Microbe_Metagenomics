import pandas as pd
import searborn as sns
class MetagenomicsPandas:
    """_summary_
    A motif parsing class for metagenomics, 
    takes a fasta file containing as many sequences
    you want and outputs a dataframe, motifs, and also
    fasta formatted motifs. This class does all you
    want to do on motif, i will also add the seqlogo
    and this also runs the alignment of the motifs.
    """
    def __init__(self, path):
        self.names = []
        self.sequences = []
        self.path = input("Enter the path of the fasta file:")
        with open("self.path", "r") as fasta:
                for i in fasta.readlines():
                    if i.startswith(">"):
                        self.names.append(i.strip().replace(">", "")
                    else:
                        self.sequences.append(i.strip())
            return self.names
            return self.sequences
    def getIDs(self):
        return self.names
    def getSequences(self):
        return self.sequences
    def getDataframe(self):
        self._dataMatrix = [(i,j) for i,j in zip(self.names, self.sequences)]
        self.matrix = pd.DataFrame(self._dataMatrix)
        return self.matrix
    def getMeanLength(self):
        return self.matrix[1].apply(len).mean()
    def getMotifs(self, start, end):
        """_summary_
        Function to check the motifs before
        extraction and for the frequency
        calculations. This functions considers that
        you want to extract the same motif from all
        the sequence at the same position. 
        """
        self._start = start
        self._end = end
        extractMotif = {}
        for k,v in dict([(i,j) for i,j in zip(self.names, self.sequences)]):
            extractMotif[k] = v[self._start:self._end]
            return extractMotif
    def prepareMotifAlignment(self, start, end, filename):
        """_summary_
        Function to prepare the fasta files with the
        extracted motifs for the alignment
        """
        self._start = start
        self._end = end
        self.motifAlignment = list(map(lambda n: n[self._start: self._end]),
                              dict([(i,j) for i,j in zip(self.names, 
                                            self.sequences)]).values())
        self.motifNames = list(self.motifAlignment.keys())
        self.motifValues = list(self.motifAlignment.values())
        for i in range(len(self.motifNames)):
            print(">"+self.motifNames[i], self.motifValues[i], sep = "\n")
            with open(self.filename, "r+") as filename:
                filename.write(">"+self.motifNames[i], self.motifValues[i], sep = "\n")
                filename.close()
        return self.motifAlignment.keys(),list(self.motifAlignment.values())
    def frequencyCalculation(self):
        """_summary_
        Function to calculate the frequency composition of the
        observed motifs.
        """
        self.dataframe_names = pd.DataFrame([(i,j) for i,j in 
                                  zip(self.names,self.sequences)])[0].to_list()
        self.dataframe_sequences = pd.DataFrame([(i,j) for i,j in 
                                  zip(self.names,self.sequences)])[1].to_list()
        self.dataframe_character_count = list(map(lambda n: [(i,n.count(i)) for i in 
                                              list(("A", "T", "G", "C"))],
                                   pd.DataFrame([(i,j) for i,j in zip(self.names,
                                                     self.sequences)])[1].to_list()))
        return [(i,j) for i,j in zip(self.dataframe_names, self.dataframe_sequences)]
    def motifCorrelation(self,features):
        """_summary_
        if you have limited features to find the correlation
        with the motifs use this function. if you have a large
        number of variables to correlate, please use the 
        motifFeatureDataframe instead
        """
        self.feature = []
        while len(feature) != len(pd.DataFrame([(i,j) for i,j in 
                                  zip(self.names,self.sequences)]))
        text = input("Enter the name of the feature variable:")
        feature.append(text)
        if text == None:
            break
        print(f'the provided features are: {self.feature}')
        self.internal = pd.DataFrame([(i,j) for i,j in 
                                               zip(self.names,self.sequences)])
        self.internal["length"] = pd.DataFrame([(i,j) for i,j in 
                                  zip(self.names,self.sequences)])[1].apply(len)
        self.feature_dataframe = pd.DataFrame(self.feature)
        self.dataframe_plot = pd.concat(self.internal,self.feature_dataframe)
        sns.pairplot(self.dataframe_plot)
        plt.show()
    def motifFeatureDataframe(self, dataframe, sep):
        """_summary_
        Function to plot the length variability of the motifs
        for each sequence with the observed variables. 
        """
        self._csv = input("Enter the path of the csv file containing 
                          the variables you want to correlate with the motifs:") 
        self.motif_dataframe = pd.read_csv(self._csv, header=True, sep = str(self.sep))
        self.combine_dataframe = pd.concat(self.motif.dataframe, self.internal)
        sns.pairplot(self.combine_dataframe)
        plt.show()
    def extractMotifUnique(self, seqID, start, end):
        """_summary_
        if you want to extract a motif from a particular sequence
        """
        self._seq_check = seqID
        self._start = start
        self._end = end 
        self.specific_sequence = [i[1] for i in ([(i, j) for i, j in 
                               zip(self.names, self.sequences)]) if i[0]=="self._seq_check"]
        motif_specific = self.specific_sequence[self._start:self._end]
        return motif_specific
    
    def countMotifPresenceAbsence(self, motif):
        """_summary_
        this function will generate a matrix for the
        presence and absence of the specific matrix
        """
        motif_tags = []
        while True:
            take_motif = input("Enter the name of the motif or motifs 
                                   you want to search in the sequences:")
           motif_tags += take_motif
           if take_motif == None:
               break
           print(f'the entered motifs are {motif_tags}')
        
        sequence_iterations = [(i, j) for i, j in zip(self.names, self.sequences)]
        matrix_count = []
        sequence_matrix_count = []
        sequence_matrix_count.append([sequence_iterations[i], 
                                      sequence_iterations[i][1].count[j] 
                       for i in range(len(sequence_iterations)) for j in matrix_tags])
        matrix_count.append([sequence_iterations[i][1].count[j] 
                             for i in range(len(sequence_iterations)) 
                             for j in matrix_tags])
     def countMotifDataframe(self, tags):
         """_summary_
         This function will count the motifs and outputs the sequence
         and the motif count  as a dataframe. 
         """
         motif_tags = []
         while True:
            tags = input("Enter the name of the motif you want 
                                 to return along with the sequence:")
            motif_tags.append(tags)
            if tags == None:
                break
            print(f'the entered motifs are:{motif_tags}')
            for i in range(len(motif_tags)):
                print(pd.DataFrame([(i,j) for i,j in zip(self.names,self.sequences)])[1]
                  [pd.DataFrame([(i,j) for i,j in zip(self.names,self.sequences)])
                                       [1].str.contains(motif_tags[i])])
     def getMotifsCoordinates(self, motif):
         """_summary_
         Get the start and end coordinates of the motifs
         where the motif occurs in the sequences
         """
         self._motif = str(motif)
         self._motif_sequences = pd.DataFrame([(i,j) for i,j in zip(self.names,self.sequences)])[1]
                           [pd.DataFrame([(i,j) for i,j in zip(self.names,self.sequences)])[1]
                                             .str.contains("self._motif")].tolist()
         self._motif_coordinate_start = []
            for i in range(len(self._motif_sequences)):
                motif_coordinate_start.append(self._motif_sequences[i].find("self._motif"))
         self._motif_coordinate_end = [self._motif_coordinate_start[i]+len("self._motif") 
                                       for i in range(len(self._motif_coordinate_start))]
         print(f'the start and end of the motif in the sequences are: 
                       {[(i,j) for i in zip(self._motif_coordinate_start, 
                                                          self._motif_coordinate_end)]}')