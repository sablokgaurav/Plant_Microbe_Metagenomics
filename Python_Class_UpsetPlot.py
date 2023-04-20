
class UpsetPlot:
    """_summary_
    A python classifier for plotting the
    MAGs abundance from the metagenomics
    classification 
    """
    def __init__(self, file1, file2):
        """_summary_

        Args:
            file (_dataframe_): _description_
            a dataframe classification file
            from the metagenomics experiment
            for the classification
        """
        self.file1 = input("file_path_file1")
        self.file2 = input("file_path_file2")
        self.classify1 = pd.read_csv("self.file1", sep=",", headers=TRUE)
        self.classify2 = pd.read_csv("self.file",sep=",", headers=TRUE)
        print('the_summary_information_of_the_file:{self.classify1.info()}')
        print('the_summary_information_of_the_file:{self.classify1.describe()}')
        print('the_summary_information_of_the_file:{self.classify1.describe}')
        print('the_summary_information_of_the_file:{self.classify1.nunique()}')

    def upsetPlot(self, bacterial_column, count_column, plotindex):
        """_summary_

        Args:
            bacterial_column (_type_): _description_
            select the bacterial species abundance column
            count_column (_type_): _description_
            select the count columns
            plotindex (_type_): _description_
            specify the plot type
        """
        bacterial_labels_file = []
        with open(bacterial_file, "r+") as bacterial:
            bacteria_text = bacterial.read()
            bacterial_labels_file.append(re.split(r'"\n"',bacterial_text))
            bactrial.close()
            length_classify = []
        for i in range(len(bacterial_labels_file)):
            length_classify.append([bacterial_label_file[i],
                                  dataframe[dataframe['bacterial_column'].
                             str.contains(bacterial_labels_file[i])==True]])
            print(f'the_number_of_the_observed_bacterial_count_are:{length_clasify}')
        bacterial_labels_1 = []
        bacterial_labels_2 = []
        bacterial_labels_1.append(self.classify1[self.bacterial_column].tolist())
        bacterial_labels2.append(self.classify2[self.bacterial_column].toList())
        bacterial_counts1 = []
        bacterial_counts2 = []
        bacterial_counts1.append(self.classify1[self.count_column].tolist())
        bacterial_counts2.append(self.classify2[self.count_column].tolist())
        final_labels = [bacterial_labels_1,bacterial_labels_2]
        final_counts = [bacterial_counts1, bacterial_counts2]
        plot = ["violin", "box", "swarm"]
        for i in range(len(plot)):
            if plot[i] == "violin" and plotindex == "violin":
                bacterial_plot = UpSetPlotly(samples=final_counts, 
                sample_names=final_labels, 
                order_by='decreasing',plot_type='box')
                usp.plot()
            elif plot[i] == "box" and plotindex == "box":
                bacterial_plot = UpSetPlotly(samples=final_counts, 
                sample_names=final_labels, 
                order_by='decreasing',plot_type='box')
                usp.plot()
            elif plot[i] == "swarm" and plotindex == "swarm":
                bacterial_plot = UpSetPlotly(samples=final_counts, 
                sample_names=final_labels, 
                order_by='decreasing',plot_type='box')
                usp.plot()
            else:
                print("No_plot_type_selected")
    def countBacteria(self):
        """_summary_
        count the abundance and make a list of the
        count and the species for the bar, stack
        plot
        """
        for i in range(len(bacterial_labels_file)): 
            plot = self.classify1[self.classify1['columnB'].
                    str.contains(bacterial_label_file[i])==True]
            count = self.classify1[self.classify1['columnB'].
                      str.contains(bacterial_label_file[i])==True].count()
            return plot,count




