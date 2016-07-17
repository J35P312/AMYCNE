import glob
import os
import math
import numpy

#read the coverage file, store the coverage for each bin together with the gc of each bin
def coverage_tab(file,data):

    for chromosome in data["chromosomes"]:
        data[chromosome]["coverage"] =[]

    for line in open(file):
        if not line[0] == "#":
            content=line.strip().split("\t")
            content[1]=int(content[1])
            content[3]=float(content[3])
            if content[0] in data:
                data[content[0]]["coverage"].append(content[3])
                
    return data


#read the gc content-tab file
def gc_tab(file):
    data={}
    bin_size=-1
    data["chromosomes"] = []
    
    for line in open(file):
        if not line[0] == "#":
            content=line.strip().split("\t")
            content[3]=float(content[3])
            if not "bin_size" in data:
                data["bin_size"] = int(content[2])-int(content[1])
            if not content[0] in data["chromosomes"]:
                data["chromosomes"].append(content[0])
                
            if not content[0] in data:
                data[content[0]] ={}
                data[content[0]]["GC"] =[]
           
            data[content[0]]["GC"].append(content[3])

    return data

#compute a gc histogram
def gc_hist(data,coverage_cutoff,size_cutoff):
    gc_dictionary={}
    #mean_values={}
    #std={}
    #for chromosome in data["chromosomes"]:
    #    mean_values[chromosome] = numpy.mean(data[chromosome]["coverage"])
        
    for chromosome in data["chromosomes"]:
        for i in range(0,len(data[chromosome]["coverage"])):
            if "Y" in chromosome or "y" in chromosome or "X" in chromosome or "x" in chromosome:
                pass
            else:
                if not data[chromosome]["GC"][i] in gc_dictionary:
                    gc_dictionary[ data[chromosome]["GC"][i]  ]=[]
                
                if data[chromosome]["coverage"][i] > 0:
                    gc_dictionary[ data[chromosome]["GC"][i]  ].append(data[chromosome]["coverage"][i])
	                    
    hist={}
    for gc in gc_dictionary:
        hist[gc]=[-1,-1]
        if len(gc_dictionary[gc]) > size_cutoff:
            bin_coverage= sum(gc_dictionary[gc])/len(gc_dictionary[gc])
            if bin_coverage < coverage_cutoff:
                hist[gc]=[bin_coverage,len(gc_dictionary[gc])]
        #print(str(gc) + " " + str(hist[gc])+ " " + str(len(gc_dictionary[gc])))
        
    return(hist)


#retrieve the data from the selected region
def regional_cn_est(Data,GC_hist,region):
    chromosome=region[0]
    start=int(region[1])
    end=int(region[2])
    if not chromosome in Data:
        chromosome=chromosome.replace("chr","")

    CN_list=[]
    GC_list=[]
    REF_list=[]

    bins = 0
    used_bins=0
    
    bin_size=Data["bin_size"]
    
    element = 0
    pos= int(math.floor(start/float(bin_size))*bin_size)
    nextpos = pos + bin_size
    while(nextpos > start and pos < end and pos/bin_size < len(Data[chromosome]["coverage"]) ):
            bins += 1
            element=int(pos/bin_size);
            bin_GC=Data[chromosome]["GC"][element]
            bin_coverage=Data[chromosome]["coverage"][element]
            if not bin_GC == -1 and not GC_hist[bin_GC][0] == -1:
                used_bins += 1
                CN_list.append(bin_coverage/GC_hist[bin_GC][0])
                GC_list.append(bin_GC)
                REF_list.append(GC_hist[bin_GC][0])

            pos+=bin_size;
            nextpos=pos+bin_size;

    if used_bins > 0:
        CN=sum(CN_list)/len(CN_list)
        GC=sum(GC_list)/len(GC_list)
        ref=sum(REF_list)/len(REF_list)
    else:
        CN=-1
        GC=-1
        ref=-1
        
        ci=-1
        
    return([CN,GC,end-start,ref,bins,used_bins,CN_list])


    
