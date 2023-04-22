from collections import Counter
def sequenceCounter(sequence):
    """_summary_
    for estimating the sequence collections
    and frequency for generating the defined
    collection types
    Args:
        sequence (_type_): _description_
        accepts a fasta sequence 
    """
    abundance_composition = []
    for fasta in SeqIO.parse("file.fasta", "fasta"):
        sequence = [" ".join(fasta.seq.upper())]
        sequence.append(fasta.name)
        abundance_composition.append(sequence)
        print(Counter(abundance_composition))
def sequenceIterator(sequence):
    """_summary_
    a sequence combinator which will generate the 
    permutations of the fasta sequence from the 
    vcf2fasta and can be comined with the gatk read 
    calling to the estimate the probability of the 
    linked variants. 

    Args:
        sequence (_type_): _description_
        takes a fasta sequence or a file of 
        fasta sequence if IO parser is defined

    Returns:
        _:str_: _description_
        returns all the iter permutations
        of the fasta for the remapping and
        the read calling.
    """
    index=sequence[0]
    for i in range(1,len(sequence)):
        indexList=[]
        for a in index:
            indexList.append(a+sequence[i]) 
            indexList.extend(a[:j] + sequence[i] + 
                       a[j:] for j in range(len(a)))
        index=indexList
    return indexList