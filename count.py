import call
import numpy
from scipy.stats import norm
from scipy import stats
import scipy

def smooth(y, box_pts):
    box = numpy.ones(box_pts)/box_pts
    y_smooth = numpy.convolve(y, box, mode='same')
    return y_smooth

def filter(Data,minimum_bin):

	
    start=-1
    for chromosome in Data["chromosomes"]:
        filtered_list=[]
        for i in range(0,len(Data[chromosome]["ratio"])):
            if start == -1 and Data[chromosome]["ratio"][i] != -1:
                start = i
            
            if start !=-1 and Data[chromosome]["ratio"][i] != -1:	
                filtered_list.append(Data[chromosome]["ratio"][i])
                
                         
                
            if start !=-1 and Data[chromosome]["ratio"][i] == -1:
                if minimum_bin % 2 == 0:
                    minimum_bin += 1
                
                
                median=scipy.signal.medfilt(filtered_list,minimum_bin)
                wiener = scipy.signal.wiener(median,minimum_bin)
                Data[chromosome]["ratio"][start:start+len(filtered_list)]=wiener
                filtered_list=[]
                start=-1
                
        if start != -1:
            Data[chromosome]["ratio"][start:start+len(filtered_list)]=scipy.signal.medfilt(filtered_list)
            filtered_list=[]
            start=-1 
                
    return(Data)



def main(Data,GC_hist,args):
    for chromosome in Data["chromosomes"]:
        Data[chromosome]["ratio"]=[]
        Data[chromosome]["len"]=0
        for i in range(0,len(Data[chromosome]["coverage"])):
           if GC_hist[Data[chromosome]["GC"][i]][0] > 0 and not Data[chromosome]["GC"][i]== -1:
                    Data[chromosome]["ratio"].append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])
                    Data[chromosome]["len"] += 1
           else:
                Data[chromosome]["ratio"].append(-1)
    
    #Data=filter(Data,101)
    ratio_hist=call.chromosome_hist(Data,args.Q)
    chrom_mean=0
    chrom_std=0 
    chrom=[]   
    for chromosome in Data["chromosomes"]:
        if not "Un_" in chromosome and not "random" in chromosome and not "GL" in chromosome and not "hap" in chromosome and not "M" in chromosome and not "X" in chromosome and not "Y" in chromosome:                
          chrom.append(ratio_hist[chromosome][0])
    chrom_mean=numpy.mean(chrom)
    chrom_std=numpy.std(chrom)
    error=1-numpy.median(chrom)
    
    print("chromosome\tcount\tCI\tgain/loss/normal\tP_val\tnormalised_coverage\tstd_normalised_coverage\traw_coverage\tstd_raw_coverage")
    for chromosome in sorted(Data["chromosomes"]):
        if not "Un_" in chromosome and not "random" in chromosome and not "GL" in chromosome and not "hap" in chromosome and not "M" in chromosome:
            gain="normal"
            count=args.ploidy*round(ratio_hist[chromosome][0]+ error,2)
            
            P_val = stats.t.sf(  abs(chrom_mean-ratio_hist[chromosome][0])/(chrom_std/numpy.sqrt(len(chrom))) , len(chrom) -1)*2
            if ratio_hist[chromosome][0] < chrom_mean - 3*chrom_std and count < args.ploidy:
                gain="loss"
            elif ratio_hist[chromosome][0] > chrom_mean + 3*chrom_std and count > args.ploidy:
                gain="gain"
                

            SEM=args.ploidy*ratio_hist[chromosome][1]/numpy.sqrt( Data[chromosome]["len"] )
            CI="({},{})".format(round(count-SEM*norm.ppf(args.p),2),round(count+SEM*norm.ppf(args.p),2) )            
            print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chromosome,count,CI,gain,P_val,ratio_hist[chromosome][0],ratio_hist[chromosome][1],ratio_hist[chromosome][2],ratio_hist[chromosome][3])
            
            
