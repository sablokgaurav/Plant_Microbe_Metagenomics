class ExtractFeatures:
    def __init__(self, file, separator = ","):
        """_summary_
        this class extracts all the genomic and 
        the metagenomic features from the genome
        prediction and prepare them for the BERT
        classification and feature prediction, you
        can also filter the BERT according to the 
        length
        Args:
            file (_type_): _description_
            feature (_type_): _description_
        """
        self.file = input("Enter the path to the file:")
        self.read = pd.read_csv("self.file", sep= ",", header=TRUE)
        genefasta = []
        with open(filename, "r+") as fname:
            for line in fname.readlines():
                genefasta.append(''.join(line.strip().split()))
                return fasta
        if feature == "gene":
            genestart = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="gene").dropna().iloc[::,1].tolist()))))
            geneend = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="gene").dropna().iloc[::,2].tolist()))))
        elif feature == "exon":
            exonstart = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="exon").dropna().iloc[::,1].tolist()))))
            exonend = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="exon").dropna().iloc[::,2].tolist()))))
        elif feature == "intron":
            intronstart = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="intron").dropna().iloc[::,1].tolist()))))
            intronend = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="intron").dropna().iloc[::,2].tolist()))))
        elif feature == tRNA:
            tRNAstart = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="tRNA").dropna().iloc[::,1].tolist()))))
            tRNAend = list(map(int,list(map(int,chrom[["version", "start", "end", "feature"]].
                       where(chrom["feature"]=="tRNA").dropna().iloc[::,2].tolist()))))
        else:
            print(f'No feature selected')

    def extractSequence(self, feature):
        if feature == "gene":
            for i in range(len(genestart)):
                gene_fasta = []
                gene_fasta.append(genefasta[genestart[i]:geneend[i]])
        elif feature == "gene":
            for i in range(len(genestart)):
                intron_fasta = []
                intron_fasta.append(genefasta[intronstart[i]:intronend[i]])
        elif feature == "gene":
            for i in range(len(genestart)):
                exon_fasta = []
                exon_fasta.append(genefasta[exonstart[i]:exonend[i]])
        elif feature == "gene":
            for i in range(len(genestart)):
                tRNA_fasta = []
                tRNA_fasta.append(genefasta[tRNAstart[i]:tRNAend[i]])
        else:
            print(f'No feature provided')

    def makeBERT(self):
        """_summary_
        prepairing the files for the BERT
        token classification and modelling. 
        You can apply the filter accordingly.
        Returns:
            _type_: _description_
        """
        geneBERT = []
        geneBERT.append([gene_fasta[:i]for i in range(len(gene_fasta))]+
             set(list(filter(None,[gene_fasta[i:i+j] for j in range(len(gene_fasta)) 
                 for i in range(len(gene_fasta)-(j-1))])))+set(list(filter(None,
                 [gene_fasta[i+1:] for j in range(len(gene_fasta)) for i in range(len(gene_fasta)-(j-1))]))))
        
        exonBERT = []
        geneBERT.append([exon_fasta[:i]for i in range(len(exon_fasta))]+
             set(list(filter(None,[exon_fasta[i:i+j] for j in range(len(exon_fasta)) 
                 for i in range(len(exon_fasta)-(j-1))])))+set(list(filter(None,
                 [gene_fasta[i+1:] for j in range(len(gene_fasta)) for i in range(len(exon_fasta)-(j-1))]))))

        intronBERT = []
        geneBERT.append([intron_fasta[:i]for i in range(len(intron_fasta))]+
             set(list(filter(None,[inton_fasta[i:i+j] for j in range(len(intron_fasta)) 
                 for i in range(len(intron_fasta)-(j-1))])))+set(list(filter(None,
                 [gene_fasta[i+1:] for j in range(len(intron_fasta)) for i in range(len(intron_fasta)-(j-1))]))))
        
        def filterBERT(self, feature, num):
            """_summary_
            filter the token generated from the features
            according to the length classification
            """
            if feature == "gene":
                filtergeneBERT = filter(lambda n: len(n) > self.num, geneBERT)
                return filtergeneBERT
            elif feature == "exon":
                exongeneBERT = filter(lambda n: len(n) > self.num, geneBERT)
                return exongeneBERT
            elif feature == "intron":
                introngeneBERT = filter(lambda n: len(n) > self.num, geneBERT)
                return introngeneBERT
            else:
                print("No_filter_applied")
