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
                    
                variants[chromosome].append([start_pos,end_pos,past_variant_type,sum(Data[chromosome]["ratio"][start_pos:i])/float(end_pos-start_pos)])
                
                past_variant_type=variant_type       
                start_pos=i
                end_pos=start_pos+1
                
    return(variants)

#merge segments    
def merge(variants,min_bins):
    merged_variants={}
    #segments that are separated by nothing will be merged
    for chromosome in variants:
        past_variant_type=-1
        merged_variant=[]
        
        for i in range(0,len(variants[chromosome])):
            variant_type = variants[chromosome][i][2]
            #print variant_type
            if past_variant_type == -1:
                past_variant_type= variant_type
                merged_variant.append( i )

            elif past_variant_type == variant_type:
                merged_variant.append( i )
            
            else:
                if not chromosome in merged_variants:
                    merged_variants[chromosome] = []
                merged_variants[chromosome].append([variants[chromosome][ merged_variant[0] ][0],variants[chromosome][merged_variant[-1]][1],past_variant_type,variants[chromosome][merged_variant[0] ][-1]])
                merged_variant=[i]
                past_variant_type= variant_type
    #deletions or duplications separated by small filtered regions or neutral regions are merged
    return(merged_variants)

def merge_similar(variants):
    merged_variants={}
    #segments that are separated by nothing will be merged
    for chromosome in variants:
        merged_variant=[]
        i=0
        while i < len(variants[chromosome])-1:
        
        
            if i == 0:
                merged_variant.append( i )
            elif (variants[chromosome][i][2] == variants[chromosome][merged_variant[-1]][2]) and  0.1 >= (variants[chromosome][i][0] - variants[chromosome][merged_variant[-1]][1])/float(variants[chromosome][i][1] - variants[chromosome][merged_variant[-1]][0] ):
                merged_variant.append( i )
            else:
                if not chromosome in merged_variants:
                    merged_variants[chromosome] = []
                merged_variants[chromosome].append([variants[chromosome][ merged_variant[0] ][0],variants[chromosome][merged_variant[-1]][1],variants[chromosome][merged_variant[-1]][2],variants[chromosome][merged_variant[0] ][-1]])
                merged_variant=[i]
            
            
            i +=1
    return(merged_variants)



#filter the data
def filter(Data,minimum_bin):
    #box = numpy.ones(10)/5
   
	
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

                    
                    
                #wiener_filtered = numpy.convolve(scipy.signal.wiener(filtered_list,5), box, mode='same')
                wiener_filtered = scipy.signal.wiener(filtered_list,minimum_bin*3)
                bin_sized_median_filt=scipy.signal.medfilt(wiener_filtered,minimum_bin)
                bin_sized_median_filt=scipy.signal.medfilt(bin_sized_median_filt,minimum_bin*2+1)
                Data[chromosome]["ratio"][start:start+len(filtered_list)]=scipy.signal.medfilt(bin_sized_median_filt,minimum_bin*3)
                filtered_list=[]
                start=-1
                
        if start != -1:
            Data[chromosome]["ratio"][start:start+len(filtered_list)]=scipy.signal.medfilt(filtered_list)
            filtered_list=[]
            start=-1 
                
    return(Data)

def main(Data,GC_hist,args):
    #compute the scaled coverage
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
    Data=filter(Data,args.min_filt)
    
    deleted_bins={}
    duplicated_bins={}         
    for chromosome in Data["chromosomes"]:
        Data[chromosome]["var"]=[];
        for i in range(0,len(Data[chromosome]["ratio"])):
            if Data[chromosome]["coverage"][i] >= 0 and Data[chromosome]["ratio"][i] >= 0 :
                #print "{}\t{}\t{}\t{}".format(chromosome,i*bin_size,(i+1)*bin_size,Data[chromosome]["ratio"][i]*GC_hist[Data[chromosome]["GC"][i]][0])
                if Data[chromosome]["ratio"][i] <= 0.7:
                    Data[chromosome]["var"].append("DEL")
                elif Data[chromosome]["ratio"][i] >= 1.3:
                    Data[chromosome]["var"].append("DUP")
                else:
                    Data[chromosome]["var"].append("NEUTRAL")
                    
            else:
                 Data[chromosome]["var"].append("FILT")
                 #print "{}\t{}\t{}\t{}".format(chromosome,i*bin_size,(i+1)*bin_size,-1)
    
    variants=segmentation(Data,args.min_bins)
    size_filtered_variants={}
    
    for chromosome in variants:
        for variant in variants[chromosome]:
            if variant[1]-variant[0] >= args.min_bins:
                if not chromosome in size_filtered_variants:
                    size_filtered_variants[chromosome] = []
                size_filtered_variants[chromosome].append(variant)
                
    variants=merge(size_filtered_variants,args.min_bins)
    
    CNV_filtered={}
    for chromosome in variants:
        for variant in variants[chromosome]:
            if variant[2] == "DUP" or variant[2] == "DEL":
                if not chromosome in CNV_filtered:
                    CNV_filtered[chromosome] = []
                CNV_filtered[chromosome].append(variant)    
    
    variants=merge_similar(CNV_filtered)
    #print the variants
    print("##fileformat=VCFv4.1")
    print("##source=AMYCNE")
    print("##ALT=<ID=DEL,Description=\"Deletion>")
    print("##ALT=<ID=DUP,Description=\"Duplication\">")
    print("##INFO=<ID=RDR,Number=1,Type=float,Description=\"Average coverage/reference ratio\">")
    print("##INFO=<ID=END,Number=1,Type=float,Description=\"The end position of the variant\">")
    print("##INFO=<ID=SVLEN,Number=1,Type=float,Description=\"The length of the variant\">")
    print("##INFO=<ID=BINS,Number=1,Type=float,Description=\"The number of bins used to call the variant\">")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    
    id_tag=0;
    for chromosome in Data["chromosomes"]:
        if chromosome in variants:
          for variant in variants[chromosome]:
            if variant[2] == "DUP" or variant[2] == "DEL":
                id_tag +=1
                firstrow= "{}\t{}\tAMYCNE_{}\tN\t<{}>\t.\tPASS".format(chromosome,bin_size*variant[0],id_tag,variant[2]) 
                info_field="END={};SVLEN={};RDR={}".format(bin_size*variant[1],(variant[1]-variant[0])*bin_size,variant[3] )
                print("\t".join([firstrow,info_field]))
