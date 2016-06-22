import glob
import os
import common

def main(Data,GC_hist,args):
    args.plody=2
    #get the coverage across the region

    coverage_list=[]
    len_list=[]
    gc_list=[]
    ref_list=[]
    nregions=0;
    for region in open(args.region):
        nregions+=1
    for region in open(args.region):
        cn, gc, length,ref =common.regional_cn_est( Data ,GC_hist, region.split("\t") ) 
        len_list.append(length)
        coverage_list.append(cn*args.plody*nregions)
        gc_list.append(gc)
        ref_list.append(ref)

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
    print("{}\t{}\t{}\t{}\t{}").format(sample,regional_gc,round(regional_coverage,1),reference_coverage,int(round(regional_coverage)))

