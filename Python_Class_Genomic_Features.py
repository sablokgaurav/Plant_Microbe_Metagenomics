# streamline functions for the training of the genomics
# data from the genomic predictions features 
# few minutes and then all the exonic predictions and ready and the code is ready 
# done and now check this and post it to git few minutes more
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
class Extract:
    """_summary_
    a python class to take the GFF3 files and prepare the 
    features plot after genome predictions and it will also
    extract those features such as gene, exon, tRNA, mRNA
    CDS and other features. It will make the features plot
    for all of them and also the visualization plot for the 
    genome pedictions
    """
    def __init__self(self, filename, outfile,feature_col):
        self.filename = input("Enter the genome prediction file")
        with open(self.filename, "r+") as fname:
            for lines in fname.readlines():
                if line[0] == "#":
                    continue
            else:
                with open(outfile, "r+") as text:
                    text.write(line)
                    text.close()     
                    self.classify = pd.read_csv("self.outfile", 
                                          sep = "\t", header = T)
def exonDifference(self):
    exonstart = list(map(int,self.classify[["version", "start", 
                                        "end", "feature"]].
                        where(self.classify["feature"]=="exon").
                                   dropna().iloc[::,1].tolist()))
    exonend = list(map(int,self.classify[["version", "start", 
                                      "end", "feature"]].
                         where(self.classify["feature"]=="exon").
                                   dropna().iloc[::,2].tolist()))
    exondiff = []
    for i in range(len(exonstart)):
        exondiff.append(exonend[i] - exonstart[i])
    sns.distplot(pd.Series(exondiff))
    sns.histplot(pd.Series(exondiff))
    plt.show()
def geneDifference(self):
    genestart = list(map(int,list(map(int,chrom[["version", "start", 
                                                 "end", "feature"]].
                                     where(chrom["feature"]=="gene")
                                       .dropna().iloc[::,1].tolist()))))
    geneend = list(map(int,list(map(int,chrom[["version", "start", 
                                              "end", "feature"]].
                                      where(chrom["feature"]=="gene").
                                        dropna().iloc[::,2].tolist()))))
    genediff = []
    for i in range(len(genestart)):
        genediff.append(geneend[i] - genestart[i])
    sns.distplot(pd.Series(genediff))
    sns.histplot(pd.Series(genediff))
    plt.show()
                    
def intronDifference(self):
    intronstart = list(map(int,chrom[["version", "start", 
                        "end", "feature"]].where(chrom["feature"]=="intron")
                                              .dropna().iloc[::,1].tolist()))
    intronend = list(map(int,chrom[["version", "start", 
                         "end", "feature"]].where(chrom["feature"]=="intron")
                                             .dropna().iloc[::,1].tolist()))
    introndiff = []
    for i in range(len(intronstart)):
        introndiff.append(intronend[i] - intronstart[i])
    sns.distplot(pd.Series(introndiff))
    sns.histplot(pd.Series(introndiff))
    plt.show()         
def tRNADifference(self):
    tRNAtart = list(map(int,chrom[["version", "start", 
                             "end", "feature"]].
                        where(chrom["feature"]=="tRNA")
                                         .dropna().iloc[::,1].tolist()))
    tRNAend = list(map(int,chrom[["version", "start", 
                                          "end", "feature"]].
                       where(chrom["feature"]=="tRNA")
                                           .dropna().iloc[::,1].tolist()))
    tRNAdiff = []
    for i in range(len(tRNAtart)):
        tRNAdiff.append(tRNAend[i] - tRNAstart[i])
    sns.distplot(pd.Series(tRNAdiff))
    sns.histplot(pd.Series(tRNAdiff))
    plt.show()
   
    number	version	feature	start	end	Unnamed: 5	annotations
0	1	TAIR10	chromosome	1	30427671	.	ID=chromosome:1;Alias=Chr1,CP002684.1,NC_003070.9
1	1	araport11	gene	10942648	10944727	-	ID=gene:AT1G30814;Name=AT1G30814;biotype=prote...
2	1	araport11	gene	29677904	29680757	-	ID=gene:AT1G78930;Name=AT1G78930;biotype=prote...
3	1	araport11	gene	26964087	26966688	+	ID=gene:AT1G71695;Name=AT1G71695;biotype=prote...
4	1	araport11	mRNA	21805932	21807496	-	ID=gene:AT1G58983;Name=AT1G58983;biotype=prote...
5	1	araport11	gene	4429718	    4430965	    +	ID=gene:AT1G12980;Name=ESR1;biotype=protein_co...
6	1	araport11	exon	17141192	17141771	+	ID=gene:AT1G45223;Name=AT1G45223;biotype=prote...
7	1	araport11	gene	21061522	21063047	-	ID=gene:AT1G56250;Name=PP2-B14;biotype=protein...
8	1	araport11	gene	24938500	24938832	-	ID=gene:AT1G66852;Name=AT1G66852;biotype=prote...
9	1	araport11	intron	26276948	26279255	-	ID=gene:AT1G69810;Name=WRKY36;biotype=protein_...
10	1	araport11	gene	27273968	27276562	-	ID=gene:AT1G72450;Name=JAZ6;biotype=protein_co...
11	1	araport11	tRNA	28617577	28622777	-	ID=gene:AT1G76280;Name=AT1G76280;biotype=prote...
12	1	araport11	gene	24406461	24409068	-	ID=gene:AT1G65630;Name=DEG3;biotype=protein_co...
13	1	araport11	gene	26906401	26908975	+	ID=gene:AT1G71390;Name=RLP11;biotype=protein_c...


