import common
import scipy.signal
import numpy
import math

#claibrate the ratio values of the sex chromosomes
def calibrate_sex(Data):
    female=True
    y_ratio=[]
    for chromosome in Data["chromosomes"]:
        if "Y" in chromosome or "y" in chromosome:
            for i in range(0,len(Data[chromosome]["ratio"])):
                if not -1 == Data[chromosome]["ratio"][i]:
                    y_ratio.append(Data[chromosome]["ratio"][i])


    if sum(y_ratio)/len(y_ratio) > 0.5:
        female = False

    if female:
        for chromosome in Data["chromosomes"]:
            if "Y" in chromosome or "y" in chromosome:
                 #add one to all the ratios on the y chromosome, any coverage found will become a duplication
                for i in range(0,len(Data[chromosome]["ratio"])):
                    if not -1 == Data[chromosome]["ratio"][i]:
                        Data[chromosome]["ratio"][i] += 1
                
    
    else:
        for chromosome in Data["chromosomes"]:
            if "Y" in chromosome or "y" in chromosome or "X" in chromosome or "x" in chromosome:
                #multiply the ratios on the x and y chromosomes by 2
                for i in range(0,len(Data[chromosome]["ratio"])):
                    if not Data[chromosome]["ratio"][i] == -1:
                        Data[chromosome]["ratio"][i] = 2*Data[chromosome]["ratio"][i]
    
    return(Data)

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
                    
                variants[chromosome].append({"start":start_pos,"end":end_pos,"type":past_variant_type,"ratio":sum(Data[chromosome]["ratio"][start_pos:i])/float(end_pos-start_pos),"bins":end_pos-start_pos})
                
                past_variant_type=variant_type       
                start_pos=i
                end_pos=start_pos+1
                
    return(variants)

def chromosome_hist(Data,Q):
    ratio_hist={}
    
    for chromosome in Data["chromosomes"]:
        bins=[]
        cov=[]
        for i in range(0,len(Data[chromosome]["ratio"])):
            if Data[chromosome]["ratio"][i] > 0 and Data[chromosome]["quality"][i] > Q and Data[chromosome]["GC"][i] > 0:
                bins.append(Data[chromosome]["ratio"][i])
                cov.append(Data[chromosome]["coverage"][i])
        if len(bins) > 0:        
            ratio_hist[chromosome]=[sum(bins)/len(bins),0,sum(cov)/len(cov),0]
            n=len(bins)
            for i in range(0,len(bins)):
                ratio_hist[chromosome][1]+=(ratio_hist[chromosome][0]-bins[i])*(ratio_hist[chromosome][0]-bins[i])/(n)
                ratio_hist[chromosome][3]+=(ratio_hist[chromosome][2]-cov[i])*(ratio_hist[chromosome][2]-cov[i])/(n)
            
            ratio_hist[chromosome][1]=math.sqrt(ratio_hist[chromosome][1])
            ratio_hist[chromosome][3]=math.sqrt(ratio_hist[chromosome][3])
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
                #compute the number of bins within this variant
                variant["bins"]=0
                for j in merged_variant:
                    variant["bins"] += variants[chromosome][j]["end"]-variants[chromosome][j]["start"]
                #compute the mean ratio of the variant
                variant["ratio"]= 0
                for j in merged_variant:
                    variant["ratio"] += variants[chromosome][j]["ratio"]*variants[chromosome][j]["bins"]
                variant["ratio"]=variant["ratio"]/variant["bins"]
                merged_variants[chromosome].append(variant)
                merged_variant=[i]
                past_variant_type= variant_type
    #deletions or duplications separated by small filtered regions or neutral regions are merged
    return(merged_variants)

def merge_similar(variants):
    merged_variants={}
    #segments that are separated by nothing will be merged
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
            sizes= [variants[chromosome][i]["end"]-variants[chromosome][i]["start"], merged_variant["end"]-merged_variant["start"]]
            size_ratio=max(sizes)/float( min(sizes))
            
            if (variants[chromosome][i]["type"] == merged_variant["type"]) and  0.1 >= variant_dist/float(variant_span) and variant_dist < 50000 and size_ratio > 0.3:
                 merged_variant["end"]=variants[chromosome][i]["end"]
                 index_list.append(i)
            else:
                if not chromosome in merged_variants:
                    merged_variants[chromosome] = []
                variant={}   
                variant["start"]=merged_variant["start"]
                variant["end"]=merged_variant["end"]
                variant["type"]=merged_variant["type"]
                #compute the number of bins within this variant
                variant["bins"]=0
                for j in index_list:
                    variant["bins"] += variants[chromosome][j]["bins"]
                #compute the mean ratio of the variant
                variant["ratio"]= 0
                for j in index_list:
                    variant["ratio"] += variants[chromosome][j]["ratio"]*variants[chromosome][j]["bins"]
                variant["ratio"]=variant["ratio"]/variant["bins"]
                        
                merged_variants[chromosome].append(variant)
                merged_variant=variants[chromosome][i]
                index_list=[i]
            
            i +=1
    return(merged_variants)
    
#filter the data
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
                
                
                small=scipy.signal.medfilt(filtered_list,minimum_bin)
                medium = scipy.signal.medfilt(small,minimum_bin*2+1)
                medium = scipy.signal.wiener(medium,minimum_bin*2+1)
                #Data[chromosome]["ratio"][start:start+len(filtered_list)]=scipy.signal.medfilt(medium,minimum_bin*3)
                Data[chromosome]["ratio"][start:start+len(filtered_list)]=medium
                filtered_list=[]
                start=-1
                
        if start != -1:
            Data[chromosome]["ratio"][start:start+len(filtered_list)]=scipy.signal.medfilt(filtered_list)
            filtered_list=[]
            start=-1 
                
    return(Data)

def main(Data,GC_hist,args):
    #compute the scaled coverage
    print("finished reading the coverage data")
    bin_size=Data["bin_size"]
    args.min_bins=int(args.min_var/bin_size)
    args.min_filt=int(args.filter/bin_size)
    for chromosome in Data["chromosomes"]:
        Data[chromosome]["ratio"]=[]
        for i in range(0,len(Data[chromosome]["coverage"])):
           if GC_hist[Data[chromosome]["GC"][i]][0] > 0 and not Data[chromosome]["GC"][i]== -1:
                    Data[chromosome]["ratio"].append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])
           else:
                Data[chromosome]["ratio"].append(-1)
                
    Data=calibrate_sex(Data)
    #filter the bins
    print("applying filters")
    Data=filter(Data,args.min_filt)
    print("computing coverage histogram")
    ratio_hist=chromosome_hist(Data,args.Q)
    print("Intrachromosomal scaling")
    for chromosome in Data["chromosomes"]:
        if ratio_hist[chromosome][0] == 0:
            continue
            
        for i in range(0,len(Data[chromosome]["ratio"])):
            Data[chromosome]["ratio"][i]=Data[chromosome]["ratio"][i]/ratio_hist[chromosome][0]
            
    deleted_bins={}
    duplicated_bins={}
    print("segmentation")
    for chromosome in Data["chromosomes"]:
        Data[chromosome]["var"]=[];
        for i in range(0,len(Data[chromosome]["ratio"])):
            if Data[chromosome]["coverage"][i] >= 0 and Data[chromosome]["ratio"][i] >= 0 :
                #print "{}\t{}\t{}\t{}".format(chromosome,i*bin_size,(i+1)*bin_size,Data[chromosome]["ratio"][i]*GC_hist[Data[chromosome]["GC"][i]][0])
                if Data[chromosome]["ratio"][i] <= 0.75:
                    Data[chromosome]["var"].append("DEL")
                elif Data[chromosome]["ratio"][i] >= 1.25:
                    Data[chromosome]["var"].append("DUP")
                else:
                    Data[chromosome]["var"].append("NEUTRAL")
                    
            else:
                 Data[chromosome]["var"].append("FILT")
                 #print "{}\t{}\t{}\t{}".format(chromosome,i*bin_size,(i+1)*bin_size,-1)
    
    variants=segmentation(Data,args.min_bins)
    

    #merge segments separated by weak signal variants
    #for chromosome in variants:
    #    i=0
    #    while i+2 < len(variants[chromosome]):
    #        merge_var= False
    #        if variants[chromosome][i]["type"] == variants[chromosome][i+2]["type"] and variants[chromosome][i+1]["type"] =="NEUTRAL":
    #            if variants[chromosome][i]["type"] == "DEL":
    #                if(variants[chromosome][i+1]["ratio"] <= 0.8) :
    #                    merge_var = True
    #            
    #            elif variants[chromosome][i]["type"] == "DUP":
    #                if(variants[chromosome][i+1]["ratio"] >= 1.2) :
    #                    merge_var = True
    #              
    #            
    #            if merge_var:
    #                variants[chromosome][i]["end"]=variants[chromosome][i+2]["end"]
    #                nbins=sum([variants[chromosome][i]["bins"],variants[chromosome][i+2]["bins"],variants[chromosome][i+1]["bins"]])
    #                rp1=variants[chromosome][i]["ratio"]*variants[chromosome][i]["bins"]
    #                rp2=variants[chromosome][i+1]["ratio"]*variants[chromosome][i+1]["bins"]
    #                rp3=variants[chromosome][i+2]["ratio"]*variants[chromosome][i+2]["bins"]
    #                variants[chromosome][i]["ratio"] = sum([rp1,rp2,rp3])/nbins
    #                variants[chromosome][i]["bins"] = nbins
    #                del variants[chromosome][i+2]
    #                del variants[chromosome][i+1]
    #                i += -1
    #                
    #                
    #        i += 1 
    
    
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
            if variant["type"] == "DUP" or variant["type"] == "DEL" or 1 == 2:
                if not chromosome in CNV_filtered:
                    CNV_filtered[chromosome] = []
                CNV_filtered[chromosome].append(variant)    
    
    variants=merge_similar(CNV_filtered)
        
    #print the variants
    print("done!")
    f=open(args.prefix,"w")

    f.write("##fileformat=VCFv4.1\n")
    f.write("##source=AMYCNE\n")
    f.write("##ALT=<ID=DEL,Description=\"Deletion>\n")
    f.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    f.write("##INFO=<ID=RDR,Number=1,Type=float,Description=\"Average coverage/reference ratio\">\n")
    f.write("##INFO=<ID=END,Number=1,Type=float,Description=\"The end position of the variant\">\n")
    f.write("##INFO=<ID=SVLEN,Number=1,Type=float,Description=\"The length of the variant\">\n")
    f.write("##INFO=<ID=BINS,Number=1,Type=float,Description=\"The number of bins used to call the variant\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    
    
    id_tag=0;
    for chromosome in Data["chromosomes"]:
        if chromosome in variants:
          for variant in variants[chromosome]:
            if variant["type"] == "DUP" or variant["type"] == "DEL" or 1 == 2:
                id_tag +=1
                firstrow= "{}\t{}\tAMYCNE_{}\tN\t<{}>\t.\tPASS".format(chromosome,bin_size*variant["start"],id_tag,variant["type"]) 
                info_field="END={};SVLEN={};RDR={};BINS={}".format(bin_size*variant["end"],(variant["end"]-variant["start"])*bin_size,variant["ratio"],variant["bins"] )
                
                CN=0
                cn, gc, length,ref,bins,used_bins,bin_list =common.regional_cn_est( Data ,GC_hist, [chromosome,variant["start"]*bin_size,variant["end"]*bin_size],args.Q )
                info_field += ";CN={}".format( int(round(cn*args.plody)))
                if int(round(cn*args.plody)) == args.plody:
                    firstrow= firstrow.replace("PASS","NOCNV")
                    
                if "quality" in Data[chromosome]:
                    failed_bins=0
                    for i in range(variant["start"],variant["end"]):
                        if Data[chromosome]["quality"][i] < args.Q and Data[chromosome]["GC"][i] > 0 and Data[chromosome]["ratio"][i] > 0:
                            failed_bins += 1
                    if failed_bins/float(variant["end"]-variant["start"]) > 0.9:
                        firstrow= firstrow.replace("PASS","FAIL")
                    info_field +=";QUAL={}".format( failed_bins/float(variant["end"]-variant["start"]) )
                    
                failed_bins=0
                for i in range(variant["start"],variant["end"]):
                    if Data[chromosome]["GC"][i] < 0:
                        failed_bins += 1
                if failed_bins/float(variant["end"]-variant["start"]) > 0.9:
                    firstrow= firstrow.replace("PASS","FILTER")
                info_field +=";FAILED_BINS={}".format( failed_bins/float(variant["end"]-variant["start"]) )                
                
                bins=[]
                for i in range(variant["start"],variant["end"]):
                    if Data[chromosome]["GC"][i] > 0 and Data[chromosome]["ratio"][i] > 0:
                        bins.append(Data[chromosome]["ratio"][i])
                if len(bins) > 0:
                    mean=sum(bins)/len(bins)
                    SEM=numpy.std(bins)/numpy.sqrt( len(bins) )
                    ci="({},{})".format(round(mean-SEM*1.96,2),round(mean+SEM*1.96,2))
                    info_field+=";CI={}".format(ci)
                else:
                    firstrow= firstrow.replace("PASS","LowQual")
                    ci="(-inf,inf)"
                info_field+=";CI={}".format(ci) 
                info_field+=";CHR_ratio={};CHR_ratio_std={};CHR_coverage={};CHR_coverage_std={}".format(ratio_hist[chromosome][0],ratio_hist[chromosome][1],ratio_hist[chromosome][2],ratio_hist[chromosome][3])
                f.write("\t".join([firstrow,info_field])+"\n")
    f.close()
