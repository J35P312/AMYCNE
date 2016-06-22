import glob
import os
import math

#read the coverage file, store the coverage for each bin together with the gc of each bin
def coverage_tab(file,gc):
    data={}

    for line in open(file):
        if not line[0] == "#":
            content=line.strip().split("\t")
            content[1]=float(content[1])
            content[2]=float(content[2])
            content[3]=float(content[3])
            if content[0] in data:
                data[content[0]].append(content[1:])
            else:
                data[content[0]]=[content[1:]]

    coverage_gc_data={}
    for chromosome in data:
        i=0;
        for bin in data[chromosome]:
            if chromosome in gc:
                if not chromosome in coverage_gc_data:
                    coverage_gc_data[chromosome]=[]
                coverage_gc_data[chromosome].append(data[chromosome][i])
                coverage_gc_data[chromosome][-1].append(gc[chromosome][i][-1])
            
            i+=1

    return coverage_gc_data


#read the gc content-tab file
def gc_tab(file):
    data={}


    for line in open(file):
        if not line[0] == "#":
            content=line.strip().split("\t")
            content[1]=float(content[1])
            content[2]=float(content[2])
            content[3]=float(content[3])

            if content[0] in data:
                data[content[0]].append(content[1:])
            else:
                data[content[0]]=[content[1:]]

    return data

#compute a gc histogram
def gc_hist(data):
    gc_dictionary={}
    for chromosome in data:
        for bin in data[chromosome]:
            #print bin
            if not bin[-1] in gc_dictionary:
                gc_dictionary[bin[-1]]=[]
            gc_dictionary[bin[-1]].append(bin[-2])

    hist={}
    for gc in gc_dictionary:
        hist[gc]=-1
        hist[gc]=[sum(gc_dictionary[gc])/len(gc_dictionary[gc]),len(gc_dictionary[gc])]
        #print(str(gc) + " " + str(hist[gc])+ " " + str(len(gc_dictionary[gc])))
    return(hist)


#retrieve the data from the selected region
def regional_cn_est(Data,GC_hist,region):
    chromosome=region[0]
    start=int(region[1])
    end=int(region[2])


    CN_list=[]
    GC_list=[]
    REF_list=[]

    bin_size=Data[chromosome][0][1]-Data[chromosome][0][0]
    bins = 0
    element = 0
    pos= int(math.floor(start/float(bin_size))*bin_size)
    nextpos = pos + bin_size
    while(nextpos > start and pos < end and pos/bin_size < len(Data[chromosome]) ):
            bins += 1
            element=int(pos/bin_size);
            content=Data[chromosome][element]
            
            CN_list.append(content[2]/GC_hist[content[3]][0])
            GC_list.append(content[3])
            REF_list.append(GC_hist[content[3]][0])

            pos+=bin_size;
            nextpos=pos+bin_size;


    CN=sum(CN_list)/len(CN_list)
    GC=sum(GC_list)/len(GC_list)
    ref=sum(REF_list)/len(REF_list)
    
    return([CN,GC,end-start,ref])


    
