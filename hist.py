import common
import numpy
def main(Data,args):
    average_coverage={}
    pooled_avg=0
    total_length=0
    for chromosome in Data["chromosomes"]:
        if not "Un_" in chromosome and not "random" in chromosome and not "hap" in chromosome and not "M" in chromosome:
            total_cov=0
            bins=0
            for bin in Data[chromosome]["coverage"]:
                if bin >= 0:
                    total_cov += bin
                    bins += 1
                    
            total_length += bins 
            average_coverage[chromosome]=total_cov/bins
            pooled_avg += total_cov
            
    pooled_avg = pooled_avg/total_length
    sample=args.coverage.replace(".tab","")
    sample=sample.split("/")[-1]
    
    output= "{}\t{}".format(sample,pooled_avg)
    for chromosome in sorted(Data["chromosomes"]):
        if not "Un_" in chromosome and not "random" in chromosome and not "hap" in chromosome and not "M" in chromosome:
            output += "\t" + str(average_coverage[chromosome])
    print output 
            
            
