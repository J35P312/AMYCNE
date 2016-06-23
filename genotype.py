import glob
import os
import common
import numpy

def retrieve_regions(line):
    mode=line.split("(")[0]
    regions=line.split("(")[-1]
    regions=regions.split(")")[0]
    regions=regions.split("|")
    
    region_list=[]
    for region in regions:
        chromosome=region.split(":")[0]
        pos=region.split(":")[-1]
        pos=pos.split("-")
        
        region_list.append([chromosome]+pos)
    
    return(mode,region_list)    

def main(Data,GC_hist,args):
    #get the coverage across the region

    coverage_list=[]
    len_list=[]
    gc_list=[]
    ref_list=[]
    total_bin_list=[]
    bin_count = 0
    used_bin_count = 0
    
    for line in open(args.region):
        mode,regions = retrieve_regions(line)
        for region in regions:
            cn, gc, length,ref,bins,used_bins,bin_list =common.regional_cn_est( Data ,GC_hist, region )
            
            if cn > -1:
                total_bin_list+=bin_list
                len_list.append(length)
                if mode == "sum":
                    nregions=len(regions)
                    coverage_list.append(cn*args.plody*nregions)
                elif mode == "avg":
                    coverage_list.append(cn*args.plody)
                gc_list.append(gc)
                ref_list.append(ref)
            used_bin_count += used_bins
            bin_count += bins

        scaled_cov=[]
        scaled_gc=[]
        scaled_ref=[]

        for i in range(0,len(coverage_list)):
            scaled_cov.append(coverage_list[i]*len_list[i])
            scaled_gc.append(gc_list[i]*len_list[i])
            scaled_ref.append(ref_list[i]*len_list[i])

        regional_coverage=sum(scaled_cov)/float(sum(len_list))
        regional_gc=sum(scaled_gc)/float(sum(len_list))
        reference_coverage=sum(scaled_ref)/float(sum(len_list))

        #print the results
        sample=args.coverage.replace(".tab","")
        sample=sample.split("/")[-1]
  
        if mode == "sum":
            for i in range(0,len(total_bin_list)):
                total_bin_list[i]=total_bin_list[i]*args.plody*nregions
        elif mode == "avg":
            for i in range(0,len(total_bin_list)):
                total_bin_list[i]=total_bin_list[i]*args.plody
        SEM=(numpy.std(total_bin_list)/numpy.sqrt(used_bin_count))
        ci="({},{})".format(round(regional_coverage-SEM*1.96,2),round(regional_coverage+SEM*1.96,2))
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}").format(sample,bin_count,used_bin_count/float(bin_count),reference_coverage,round(regional_coverage,2),ci,int(round(regional_coverage)),line.strip())

