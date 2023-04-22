import venn
def draw(infile,speciesC, conditionC, options, fontsize, filename):
    """_summary_
    a complete function to make the venn plots for all the metagenomics
    analysis from all the abundance analysis and the metagenomics 
    classification.

    Args:
        infile (_type_): _description_
        this argument expects a csv file containing 
        the species column and the condition column
        you dont have to filter the file or make any changes
        after the abundance, this funtion will automatically
        do the same.
        speciesC (_type_): _description_
        species columns containing the abundance
        conditionC (_type_): _description_
        condition column having the experimental conditions
        options (_type_): _description_
        options whether you want numeric or percentage
        fontsize (_type_): _description_
        font size 
        filename (_type_): _description_
        file name to write the length of the species in each dataframe
        and also to write the corresponding bacterial species and 
        also to write the corresponding dataframe for ready to go
        summary tables
    """
    dataframe = pd.read.csv("infile.csv",sep=",",header=T)
    conditions = dataframe["conditionC"].dropna().unique().tolist()
    species = speciesC
    condition = ConditionC
    options = ["numeric", "percentage"]
    plot = {}
    for i in range(len(conditions)):
    try:
        plot[conditions[i]] = set(dataframe[["species","condition"]].
               where(dataframe["condition"] == conditions[i]).
                             dropna().iloc[::,0].tolist())
        final_plot = plot
        for i in range(len(options)):
            if options[i] and options == numeric:
                venn(plot, fontsize=fontsize, legend = "upperleft")
            else:
                if options[i] and options == percentage:
                    venn(plot, fmt= "{percentage:.2f}%", fontsize = 8, legend_loc = "upper left")
        write_columns = []
           for i in range(len(conditions)):
            wtite.append([condition[i], set(dataframe[["species","condition"]].
                                 where(dataframe["condition"] == conditions[i]).
                                                 dropna().iloc[::,0].tolist())])
        length_features = []
            for i in range(len(conditions)):
                length_features.append([conditions[i], len(set(dataframe[["species","condition"]].
                                        where(dataframe["condition"] == conditions[i]).
                                                     dropna().iloc[::,0].tolist()))])
        with open(filename, "r+") as f:
            filename.write(write_columns)
            filename.write(length_features)
            filename.close()
