import common
import numpy
def main(Data,args):
    average_coverage={}
    pooled_avg=0
    total_length=0
    for chromosome in Data["chromosomes"]:
        if not "Un_" in chromosome and not "random" in chromosome and not "hap" in chromosome and not "M" in chromosome:
            average_coverage[chromosome]=numpy.average(Data[chromosome]["coverage"])
            pooled_avg += average_coverage[chromosome]*len(Data[chromosome]["coverage"])
            total_length += len(Data[chromosome]["coverage"])
    pooled_avg = pooled_avg/total_length
    sample=args.coverage.replace(".tab","")
    sample=sample.split("/")[-1]
    
    output= "{}\t{}".format(sample,pooled_avg)
    for chromosome in sorted(Data["chromosomes"]):
        if not "Un_" in chromosome and not "random" in chromosome and not "hap" in chromosome and not "M" in chromosome:
            output += "\t" + str(average_coverage[chromosome])
    print output 
            
            
