# a set of functions which can help
#with the preparation of the configuration files
class syntenyPlot:
    """_summary_
    a configuration preparation for 
    NGenomeSyn
    """
    def __init__(self, para, genome1, 
                             genome2, 
                              len1, 
                              len2, 
                              link):
        self.para = ''.join(["SetParaFor", " ", "=", " " "global"])
        self.genome1 = ''.join(["GenomeInfoFile1", " ", "=", "len1"])
        self.genome2 = ''.join(["GenomeInfoFile1", " ", "=", "len2"])
        self.link = input("link_file_path")
        print(f'the_provided_files_are:{self.para},
                        {self.genome1},{self.genome2},{self.link})

    def printTwoGenome(self):
        """_summary_
        prepare for two genomes
        """
        for i in [self.para,self.genome1, self.genome2, self.link]:
            print(i)

    def printMultipleGenome(self,num,length,
                     link, highlightGenome):
        """_summary_
        prepare for multiple genomes and also
        the genome to be highlighted
        Args:
            num (_type_): _description_
            number of genomes to be analyzed
            length (_type_): _description_
            path to the length file
            link (_type_): _description_
            path to the link file
            highlightGenome (_type_): _description_
            name of the genome to be highlighted
        """
        genome_file = [f'GenomeInfoFile{i}' for i in range(1, self.num)]
        length_file = []
        with open(self.length, "r+") as l:
            length_temp = l.read()
            length_file.append(re.split("\n"),length_temp)
            l.close()
        link_file = []
        with open(self.link, "r+") as link:
            link_temp = link.read()
            link_file.append(re.split("\n"),link_temp)
            link.close()
        if len(genome_file) != length_file:
            raise ValueError as file_not_compatible
                    print(f'ValueError'+str(file_not_compatible))
            print(f'')
        configuration = []
        configuration.append([''.join([genome_file[i]," ","="," ", 
                               length_file[i]]) for i in range(len(genome_file))])     
            for i in range(len(genome_file)):
                if i and highlightGenome == genome_file[i]:
                     cratio = "".join(["CanvasHeightRitao"," ","=","," "1.6"])
                     genomehighlight = "".join(["SetParaFor", " ", "="," "genome[i]])
                     zoom = "".join(["ZoomChr", "="," ", 0.5])
                     rotate = "".join(["RotateChr", " ", "="," ",28])
                     shiftX = "".join(["ShiftX"," ","="," ",60])
                     shiftY = "".join(["ShiftY"," ","="," ","-80"])
                    print(f'the_highlight_genome_configuration_file:{self.para},
                    {[i for i in configuration]},\ 
                         {[f'link{i}={link_file[i]}' for i in range(len(link_file))]},
                         {cratio},{genomehighlight},
                         {zoom},{rotate},{shiftX},{shiftY}
                         )




