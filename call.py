import common
import scipy.signal
import scipy.stats
import bottleneck
import numpy
import math
import sys
import os
import time

#claibrate the ratio values of the sex chromosomes
def calibrate_sex(Data):
    female=True
    y_ratio=0
    for chromosome in Data["chromosomes"]:
        if chromosome == "Y" or "chrY" == chromosome:
            y_ratio=numpy.median( Data[chromosome]["ratio"][ numpy.where(Data[chromosome]["ratio"] >= 0 ) ])
            continue

    if y_ratio > 0.3:
        female = False

    if female:
        for chromosome in Data["chromosomes"]:
            if chromosome == "Y" or "chrY" == chromosome:
                 #add one to all the ratios on the y chromosome, any coverage found will become a duplication
                for i in range(0,len(Data[chromosome]["ratio"])):
                    if not -1 == Data[chromosome]["ratio"][i]:
                        Data[chromosome]["ratio"][i] += 1
    else:
        for chromosome in Data["chromosomes"]:
            if chromosome == "Y" or "chrY" == chromosome or "X" == chromosome or "chrX" == chromosome:
                #multiply the ratios on the x and y chromosomes by 2
                for i in range(0,len(Data[chromosome]["ratio"])):
                    if not Data[chromosome]["ratio"][i] == -1:
                        Data[chromosome]["ratio"][i] = 2*Data[chromosome]["ratio"][i]
    
    return(Data)

def retrieve_phred(bins,ratio,percentiles):
    pvals=[]
    for b in bins:
        percentile=0
        for val in percentiles:

            if val > b:
                break
            percentile+=0.001

        if ratio > 1:
            percentile=1-percentile
            if percentile < 0:
                percentile=0

        if not percentile:
            percentile=0.000001
        pvals.append(percentile)

    p= scipy.stats.combine_pvalues(pvals)[1]
    p=p*len(pvals)
    if not p or p < 0:
        return(1000)
    phred=int(round(-math.log10(p)))
    if phred > 1000:
        return 1000

    return(phred)

def retrieve_phred_non_param(nbins,ratio,data,ratio_hist):
    pvals=[]
    n=10000
    p=0
    sampled_chromosomes=[]
    #pick chromosomes big enough to sample from
    for chromosome in data["chromosomes"]:
        if len(data[chromosome]["ratio"]) > 10*nbins and not "X" in chromosome and not "Y" in chromosome and abs(ratio_hist[chromosome][0]-1) < 0.1:
            sampled_chromosomes.append(chromosome)

    chromosomes=list(sorted(numpy.random.choice(sampled_chromosomes,size=n)))
    simulated_positions=[]
    #simulate 
    for chromosome in sorted(sampled_chromosomes):
        simulated_positions+=list(numpy.random.randint(0,high=len(data[chromosome]["ratio"])-nbins,size=chromosomes.count(chromosome)))
    failed=0
    for i in range(0,n):
        chromosome=chromosomes[i]
        pos=simulated_positions[i]
        sim_bins=data[chromosome]["ratio"][pos:pos+nbins]

        if list(sim_bins).count(-1)/float(len(sim_bins)) >= 0.6:
            failed+=1
            continue

        sim_ratio=bottleneck.median(sim_bins[numpy.where(sim_bins >= 0)])
        if ratio > 1 and sim_ratio >= ratio:
            p+=1
        elif ratio < 1 and sim_ratio <= ratio:
            p +=1

    if failed == n:
        return 1000

    p=p/float(n-failed)
    if not p:
        return(int(-10*math.log10(1/float(n-failed))))
    #normalise between 1000 and 1
    phred=int(-10*math.log10(p))
    return(phred)


#divide all bins into segments, these segments are dup,del,normal, or filtered
def segmentation(Data,minimum_bin):
    variants={}
    
    for chromosome in Data["chromosomes"]:
        start_pos=-1;
        end_pos=-1;
        variant_type=None
        past_variant_type=-1
        for i in range(0,len(Data[chromosome]["var"])):
            variant_type=Data[chromosome]["var"][i]
            
            if past_variant_type == -1:
                start_pos=i
                end_pos = i+1
                past_variant_type=variant_type
            elif past_variant_type == variant_type:
                end_pos +=1
            else:
                if not chromosome in variants:
                    variants[chromosome] = []
                ratio_list=Data[chromosome]["ratio"][start_pos:end_pos+1]
                ratio_list=ratio_list[numpy.where(ratio_list >= 0)]

                variants[chromosome].append({"start":start_pos,"end":end_pos,"type":past_variant_type,"ratio":numpy.average(ratio_list),"bins":end_pos-start_pos,"ratio_list":list(ratio_list)})
                
                ratio_list=[]
                past_variant_type=variant_type       
                start_pos=i
                end_pos=start_pos+1
                
    return(variants)

def mad(arr):
    #copied from:https://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
    arr = numpy.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = numpy.median(arr)
    return numpy.median(numpy.abs(arr - med))

def chromosome_hist(Data,Q):
    ratio_hist={}
    
    for chromosome in Data["chromosomes"]:
        bins=[]
        cov=[]
        for i in range(0,len(Data[chromosome]["ratio"])):
            if len(Data[chromosome]["quality"]):
                #if Data[chromosome]["ratio"][i] > 0 and Data[chromosome]["quality"][i] > Q and Data[chromosome]["GC"][i] > 0:
                if Data[chromosome]["ratio"][i] >= 0 and Data[chromosome]["GC"][i] > 0:
                   bins.append(Data[chromosome]["ratio"][i])
                   cov.append(Data[chromosome]["coverage"][i])
            else:
                if Data[chromosome]["ratio"][i] >= 0 and Data[chromosome]["GC"][i] > 0:
                   bins.append(Data[chromosome]["ratio"][i])
                   cov.append(Data[chromosome]["coverage"][i])

        if len(bins) > 0:        
            ratio_hist[chromosome]=[numpy.median(bins),0,numpy.median(cov),0]
            n=len(bins)
            for i in range(0,len(bins)):
                ratio_hist[chromosome][1]+=(ratio_hist[chromosome][0]-bins[i])*(ratio_hist[chromosome][0]-bins[i])/(n)
                ratio_hist[chromosome][3]+=(ratio_hist[chromosome][2]-cov[i])*(ratio_hist[chromosome][2]-cov[i])/(n)
            
            ratio_hist[chromosome][1]=mad(bins)
            ratio_hist[chromosome][3]=mad(cov)
        else:
            ratio_hist[chromosome]=[0,0,0,0]
    return(ratio_hist)
    
#merge segments    
def merge(variants,min_bins):
    merged_variants={}
    #segments that are separated by nothing will be merged
    for chromosome in variants:
        past_variant_type=-1
        merged_variant=[]
        
        for i in range(0,len(variants[chromosome])):
            variant_type = variants[chromosome][i]["type"]
            #print variant_type
            if past_variant_type == -1:
                past_variant_type= variant_type
                merged_variant.append( i )
            elif past_variant_type == variant_type:
                merged_variant.append( i )
            
            else:
                if not chromosome in merged_variants:
                    merged_variants[chromosome] = []

                variant={}
                variant["start"]=variants[chromosome][ merged_variant[0] ]["start"]
                variant["end"]=variants[chromosome][merged_variant[-1]]["end"]
                variant["type"]=past_variant_type
                bins=[]
                for var in merged_variant:
                    bins += variants[chromosome][var]["ratio_list"]
                variant["ratio_list"]=bins
                variant["bins"]=len(variant["ratio_list"])
                variant["ratio"]=numpy.median(variant["ratio_list"])
                merged_variants[chromosome].append(variant)
                merged_variant=[i]

                past_variant_type= variant_type
    return(merged_variants)

def merge_similar(variants):
    merged_variants={}
    #closely located segments are merged
    for chromosome in variants:
        i=0
        index_list=[]
        while i < len(variants[chromosome])-1:

            if i == 0:
                merged_variant=variants[chromosome][i]
                index_list.append(i)
                i += 1
                continue
                
            variant_span= variants[chromosome][i]["end"] - merged_variant["start"]
            variant_dist= variants[chromosome][i]["start"] - merged_variant["end"]
            
            if (variants[chromosome][i]["type"] == merged_variant["type"]) and  0.1 >= variant_dist/float(variant_span) and variant_dist < 20000:
                 merged_variant["end"]=variants[chromosome][i]["end"]
                 merged_variant["ratio_list"]+=variants[chromosome][i]["ratio_list"]
                 index_list.append(i)
            else:
                if not chromosome in merged_variants:
                    merged_variants[chromosome] = []
                variant={}   
                variant["start"]=merged_variant["start"]
                variant["end"]=merged_variant["end"]
                variant["type"]=merged_variant["type"]
                variant["ratio_list"]=merged_variant["ratio_list"]
                variant["bins"]=len(variant["ratio_list"])
                variant["ratio"]=numpy.median(variant["ratio_list"])
                merged_variants[chromosome].append(variant)
                merged_variant=variants[chromosome][i]
                index_list=[i]
            
            i +=1
    return(merged_variants)
    
#filter the data
def filter(Data,minimum_bin):
    if minimum_bin % 2 == 0:
        minimum_bin += 1
	
    start=-1
    for chromosome in Data["chromosomes"]:
        filtered_list=[]
        for i in range(0,len(Data[chromosome]["ratio"])):
            if start == -1 and Data[chromosome]["ratio"][i] != -1:
                start = i
            
            if start !=-1 and Data[chromosome]["ratio"][i] != -1:	
                filtered_list.append(Data[chromosome]["ratio"][i])
                                
            if start !=-1 and Data[chromosome]["ratio"][i] == -1:
                if len(filtered_list) < minimum_bin:
                    remove=numpy.zeros(len(filtered_list))-1
                    Data[chromosome]["ratio"][start:start+len(filtered_list)]=remove
                else:
                    small=scipy.signal.medfilt(filtered_list,minimum_bin)
                    #Data[chromosome]["ratio"][start:start+len(filtered_list)]=scipy.signal.wiener(small,minimum_bin)
                    Data[chromosome]["ratio"][start:start+len(filtered_list)]=small
                filtered_list=[]
                start=-1
                
        if start != -1:
            Data[chromosome]["ratio"][start:start+len(filtered_list)]=scipy.signal.medfilt(filtered_list,minimum_bin)
            filtered_list=[]
            start=-1 
                
    return(Data)

def coverage_hist(Data,ratio_hist):

    hist=[]
    for chromosome in Data["chromosomes"]:
        if ratio_hist[chromosome][0] < 0.5:
            continue
            
        for i in range(0,len(Data[chromosome]["ratio"])):
            if Data[chromosome]["ratio"][i] < 0:
                continue
            if not i % 100 and not numpy.isnan(Data[chromosome]["ratio"][i]):
                hist.append(Data[chromosome]["ratio"][i])

    hist=numpy.array(hist)
    return(hist)

def main(Data,GC_hist,args):
    #compute the scaled coverage
    print("finished reading the coverage data")
    bin_size=Data["bin_size"]
    args.min_bins=int(args.nbins/2)


    if not args.min_bins:
        print "Error: the minimum variant size is smaller than the bin sie of the input data!"
        quit()
			
    for chromosome in Data["chromosomes"]:
        Data[chromosome]["ratio"]=[]
        for i in range(0,len(Data[chromosome]["coverage"])):
           if GC_hist[Data[chromosome]["GC"][i]][0] > 0 and not Data[chromosome]["GC"][i]== -1:
                    if Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0] < args.max:
                        Data[chromosome]["ratio"].append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])
                    else:
                        Data[chromosome]["ratio"].append(-1)                        
           else:
                Data[chromosome]["ratio"].append(-1)
        Data[chromosome]["ratio"]=numpy.array(Data[chromosome]["ratio"])
                
    Data=calibrate_sex(Data)
    #filter the bins
    print("applying filters")
    Data=filter(Data,args.nbins*2)
    print("computing coverage histogram")
    ratio_hist=chromosome_hist(Data,args.Q)

    hist=coverage_hist(Data,ratio_hist)
    percentiles=numpy.percentile(hist,numpy.array(range(0,1001))/10.0)
    overall_sd=numpy.std(hist[ numpy.where(hist <= 2) ])
  
    print("segmentation")
    for chromosome in Data["chromosomes"]:
        Data[chromosome]["var"]=numpy.repeat( "NEUTRAL",len(Data[chromosome]["ratio"]) );
        Data[chromosome]["ratio"]=numpy.array(Data[chromosome]["ratio"])

        ratio_indexes=[]
        ratios=[]
        for i in range(1,len(Data[chromosome]["ratio"])):
            if Data[chromosome]["ratio"][i] >= 0:
                ratio_indexes.append(i)
                ratios.append(Data[chromosome]["ratio"][i])
        differences=[]
        for i in range(1,args.nbins+1):
            tmp=[]
            for j in range(0,len(ratios)-args.nbins):
                tmp.append( abs(ratios[j]-ratios[i+j]))
            differences.append(tmp)
        differences=numpy.array(differences)

        change_points=[]
        #print len(ratios)

        lim=overall_sd*2
        for i in range(0,len(ratios)-args.nbins):
            changes=differences[:,i]
            #print "{} {}".format(lim,numpy.min(changes))
            if bottleneck.median(changes) > lim and numpy.std(changes[1:]) < overall_sd:
                #print "{} {}".format(lim,numpy.min(changes))
                change_points.append(ratio_indexes[i])

        segments=[]
        change_points.append( len( Data[chromosome]["ratio"] ) )        
        for i in range(0,len(change_points)):
            if i == 0:
                segments.append(range(0,change_points[i]))
            elif i != len(change_points)-1:
                segments.append(range(change_points[i-1],change_points[i]))
            else:
                segments.append(range(change_points[i-1],len(Data[chromosome]["ratio"])))

        for segment in segments:
            segment_intensities= Data[chromosome]["ratio"][segment]
            non_filt_bins=segment_intensities[numpy.where(segment_intensities >= 0)]
            TYPE="NEUTRAL"
            med=bottleneck.median(non_filt_bins)

            if len(non_filt_bins) < args.min_bins: 
                TYPE="FILT"
            elif med <= 1-0.5/args.plody:
                TYPE="DEL"
            elif med >= 1+0.5/args.plody:
                TYPE="DUP"             
            Data[chromosome]["var"][segment]=TYPE

    print("merging")
    variants=segmentation(Data,args.min_bins)
    
    size_filtered_variants={}   
    for chromosome in variants:
        for variant in variants[chromosome]:
            if variant["bins"] >= args.min_bins:
                if not chromosome in size_filtered_variants:
                    size_filtered_variants[chromosome] = []
                size_filtered_variants[chromosome].append(variant)
  
    variants=merge(size_filtered_variants,args.min_bins)
    
    CNV_filtered={}
    for chromosome in variants:
        for variant in variants[chromosome]:
            if variant["type"] == "DUP" or variant["type"] == "DEL":
                if not chromosome in CNV_filtered:
                    CNV_filtered[chromosome] = []
                CNV_filtered[chromosome].append(variant)    
    
    #read the bam header
    args.contigs={}
    args.contig_order=[]
    if args.bam:
        with os.popen("samtools view -H {}".format(args.bam)) as pipe:
            for line in pipe:
                if line[0] == "@":
                    if "SN:" in line:
                        content=line.strip().split()
                        chromosome=content[1].split("SN:")[-1]
                        length=content[2].split("LN:")[-1]
                        args.contigs[chromosome]=length
                        args.contig_order.append(chromosome)
                    elif "\tSM:" in line and not args.sample:
                        args.sample=line.split("\tSM:")[-1].split("\t")[0].strip()

    #print the variants
    print("computing statistics")

    vals=[]
    counts={}
    for chromosome in Data["chromosomes"]:
        if chromosome in variants:
          for variant in variants[chromosome]:
            if variant["type"] == "DUP" or variant["type"] == "DEL" or 1 == 2:
                phred_non_param=retrieve_phred_non_param(variant["bins"],variant["ratio"],Data,ratio_hist)
                if not phred_non_param in counts:
                    vals.append(phred_non_param)
                    counts[phred_non_param]=0
                counts[phred_non_param]+=1
                variant["pred_non_param"]=phred_non_param

    args.scoren=0
    n=0
    for val in sorted(vals,reverse=True):
        #print "{} {}".format(val,counts[val])
        if n+counts[val] < args.Evar:
            args.scoren=val
            n+=counts[val]
        else:
            break              

    f=open(args.output,"w")

    f.write("##fileformat=VCFv4.1\n")
    f.write("##source=AMYCNE\n")
    f.write("##ALT=<ID=DEL,Description=\"Deletion>\n")
    f.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    f.write("##INFO=<ID=RDR,Number=1,Type=Float,Description=\"Average coverage/reference ratio\">\n")
    f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"The end position of the variant\">\n")
    f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"The length of the variant\">\n")
    f.write("##INFO=<ID=BINS,Number=1,Type=Integer,Description=\"The number of bins used to call the variant\">\n")
    f.write("##INFO=<ID=SCOREF,Number=1,Type=Integer,Description=\"The variant score produced from Fishers method\">\n")
    f.write("##INFO=<ID=SCOREN,Number=1,Type=Integer,Description=\"The variant score produced from non-parametric sampling method\">\n")
    f.write("##INFO=<ID=QUAL,Number=1,Type=Float,Description=\"The fraction of low quality bins\">\n")
    f.write("##INFO=<ID=FAILED_BINS,Number=1,Type=Float,Description=\"The fraction of filtered bins\">\n")
    f.write("##INFO=<ID=ratio,Number=1,Type=Float,Description=\"Normalised coverage across the chromosome\">\n")
    f.write("##INFO=<ID=ratioMAD,Number=1,Type=Float,Description=\"normalised Median absolute deviation across the chromosome\">\n")
    f.write("##INFO=<ID=coverage,Number=1,Type=Float,Description=\"Median coverage of the chromosome\">\n")
    f.write("##INFO=<ID=coverageMAD,Number=1,Type=Float,Description=\"Median absolute deviation of the coverage across the chromosome\">\n")
    if args.contig_order:
        for contig in args.contig_order:
            f.write("##contig=<ID={},length={}>\n".format(contig,args.contigs[contig]))
    f.write("##FILTER=<ID=LowBinQual,Description=\"More than 90% of the bins have less than {} mapping quality\">\n".format(args.Q))
    f.write("##FILTER=<ID=RegionFilter,Description=\"More than 90% of the bins are flagged extremed GC and/or mapping quality\">\n")
    f.write("##FILTER=<ID=RatioFilter,Description=\"The RD ratio is less than 2 sd of the RD, or RDR higher than ratiolim\">\n")
    f.write("##FILTER=<ID=LowScore,Description=\"Low variant score\">\n")
    f.write("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n")
    f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    f.write("##nbins={} RDstdev={} ScoreNLimit={}\n".format(args.nbins,overall_sd,args.scoren))
    f.write("##AMYCNEcmd=\"{}\"\n".format(" ".join(sys.argv)))
    if not args.sample:
        args.sample=args.coverage.split("/")[-1].split(".")[0]
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(args.sample))
    format_column="GT:CN"
    id_tag=0;
    for chromosome in Data["chromosomes"]:
        if chromosome in variants:
          for variant in variants[chromosome]:
            if variant["type"] == "DUP" or variant["type"] == "DEL" or 1 == 2:
                id_tag +=1
                filt="PASS"
                info_field="END={};SVLEN={};RDR={};BINS={}".format(bin_size*variant["end"],(variant["end"]-variant["start"]+1)*bin_size,variant["ratio"],variant["bins"] )
                CN=int(round(variant["ratio"]*args.plody))
                if "quality" in Data[chromosome]:
                    failed_bins=0
                    for i in range(variant["start"],variant["end"]):
                        if Data[chromosome]["quality"][i] < args.Q and Data[chromosome]["GC"][i] > 0 and Data[chromosome]["ratio"][i] > 0:
                            failed_bins += 1
                    if failed_bins/float(variant["end"]-variant["start"]) > 0.9:
                        filt="LowBinQual"
                    info_field +=";QUAL={}".format( failed_bins/float(variant["end"]-variant["start"]) )

                phred=retrieve_phred(variant["ratio_list"],variant["ratio"],percentiles)
                phred_non_param=variant["pred_non_param"]
                info_field+=";SCOREF={};SCOREN={}".format(phred,phred_non_param)
                #info_field+=";SCOREF={}".format(phred)

                if phred < args.scoref or phred_non_param < args.scoren:
                    filt="LowScore"                      
                
                failed_bins=0
                for i in range(variant["start"],variant["end"]):
                    if Data[chromosome]["ratio"][i] < 0:
                        failed_bins += 1
                if failed_bins/float(variant["end"]-variant["start"]) > 0.9:
                    filt="RegionFilter"
                if abs(variant["ratio"]-1) <= overall_sd*2 or abs(variant["ratio"]) > args.ratioLim:
                    filt="RatioFilter"

                info_field +=";FAILED_BINS={}".format( failed_bins/float(variant["end"]-variant["start"]) )                
                
                mean=numpy.average(variant["ratio_list"])
                SEM=numpy.std(variant["ratio_list"])/numpy.sqrt( len(variant["ratio_list"]) )
                ci="({},{})".format(round(mean-SEM*3,2),round(mean+SEM*3,2))

                firstrow = "{}\t{}\tAMYCNE_{}\tN\t<{}>\t{}\t{}".format(chromosome,bin_size*variant["start"],id_tag,variant["type"],phred_non_param,filt)
                info_field+=";ratio={};ratioMAD={};coverage={};coverageMAD={}".format(ratio_hist[chromosome][0],ratio_hist[chromosome][1],ratio_hist[chromosome][2],ratio_hist[chromosome][3])
                alt=abs((CN-args.plody))
                if alt > args.plody:
                    alt=args.plody
                ref=args.plody-alt
                genotype="/".join(["0"]*ref+["1"]*alt)
                format_field="{}\t{}:{}".format(format_column,genotype,CN)
                f.write("\t".join([firstrow,info_field,format_field])+"\n")
    f.close()
