import argparse
               
parser = argparse.ArgumentParser("""Prints tab separated GC content in windows of user set size""")
parser.add_argument('--fa' , type=str, required=True,help="the fasta file")
parser.add_argument('--index' , type=int,default=0,help="the coordinate of the first base of each chromosome , usually 0 or 1(default = 0)")
parser.add_argument('--size' , type=int, default=1000, help="the size of the bins")
parser.add_argument('--replace_chr' , action="store_true", help="remove chr prefix(i.e chr1 -> 1)")
parser.add_argument('--n_mask', action="store_true" , help="regions having more than 20%% content set to N or n is set to -1")
args = parser.parse_args()


#read the fasta file
reference={}
chromosomes=[]
with open(args.fa, 'r+') as f:
    sequence=f.read()


split_reference=sequence.split(">")
del sequence
del split_reference[0]
#store the reference as a dictionary
for chromosome in split_reference:
    content=chromosome.split("\n",1)
    reference[content[0]]=content[1].replace("\n","")
    chromosomes.append(content[0])
        
#calculate the coverage of each bin
for chromosome in chromosomes:
    start=0
    end=args.size
    while start < len(reference[chromosome]):
        printed_chromosome=chromosome
        if args.replace_chr:
            printed_chromosome=chromosome.replace("chr","")
        region=reference[chromosome][start:end]
        G=region.count("c")+region.count("C")+region.count("g")+region.count("G")
        A=region.count("A")+region.count("a")+region.count("t")+region.count("T")
        N=region.count("N")+region.count("n")
        if A+G > 0:
            GC=float(G)/float(A+G)
        elif args.n_mask and N > 0.2:
            GC=-1
        else:
            GC=0
        print("{}\t{}\t{}\t{}".format(printed_chromosome,start+args.index,end+args.index, round(GC,2) ))
        start=end
        if end+args.size < len(reference[chromosome]):
            end += args.size
        else:
            end=len(reference[chromosome])
