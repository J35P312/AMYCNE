import call
import numpy
from scipy.stats import norm
from scipy import stats
import scipy

def smooth(y, box_pts):
    box = numpy.ones(box_pts)/box_pts
    y_smooth = numpy.convolve(y, box, mode='same')
    return y_smooth

def bootstrap_CI(ratio,n,sz,CN):
    dist=numpy.random.choice(ratio,size=(n,sz))
    medians=[]
    for d in dist:
        medians.append(numpy.median(d))
    p=numpy.percentile(medians,[0.1,99.9])
    return round(abs(p[0]*2-p[1]*2),2)

def main(Data,GC_hist,args):
    for chromosome in Data["chromosomes"]:
        Data[chromosome]["ratio"]=[]
        Data[chromosome]["len"]=0
        for i in range(0,len(Data[chromosome]["coverage"])):
           if len(Data[chromosome]["quality"]):
                if Data[chromosome]["quality"][i] < args.Q and Data[chromosome]["coverage"][i]:
                    continue
 
           if GC_hist[ round(Data[chromosome]["GC"][i],2) ][0] > 0 and not Data[chromosome]["GC"][i]== -1:
                    Data[chromosome]["ratio"].append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])
                    Data[chromosome]["len"] += 1
        Data[chromosome]["ratio"]=numpy.array(Data[chromosome]["ratio"])

    ratio_hist=call.chromosome_hist(Data,args.Q)
    chrom_mean=0
    chrom_std=0 
    chrom=[]   
    for chromosome in Data["chromosomes"]:
        if not "Un_" in chromosome and not "random" in chromosome and not "GL" in chromosome and not "hap" in chromosome and not "M" in chromosome and not "X" in chromosome and not "Y" in chromosome:                
          chrom.append(ratio_hist[chromosome][0])
    chrom_mean=numpy.median(chrom)
    chrom_std=numpy.std(chrom)

    print("chromosome\tcount\taccuracy\tgain/loss/normal\tnormalised_median\tstd_normalised_coverage\traw_coverage\tstd_raw_coverage")
    for chromosome in sorted(Data["chromosomes"]):
        if not len (Data[chromosome]["ratio"]):
            continue
        if not "Un_" in chromosome and not "random" in chromosome and not "GL" in chromosome and not "hap" in chromosome and not "M" in chromosome:
            gain="normal"
            count=args.ploidy*round(ratio_hist[chromosome][0],1)
            
            if ratio_hist[chromosome][0] < chrom_mean - 3*chrom_std and count < args.ploidy and ratio_hist[chromosome][0] < 1-args.d:
                gain="loss"
            elif ratio_hist[chromosome][0] > chrom_mean + 3*chrom_std and count > args.ploidy and ratio_hist[chromosome][0] > 1+args.d:
                gain="gain"
                
            #CI="({},{})".format(round(count-ratio_hist[chromosome][1],round(count+ratio_hist[chromosome][1],2) )
            CI=bootstrap_CI(Data[chromosome]["ratio"],10000,10000,count)            
            print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chromosome,count,CI,gain,ratio_hist[chromosome][0],ratio_hist[chromosome][1],ratio_hist[chromosome][2],ratio_hist[chromosome][3])
            
            
