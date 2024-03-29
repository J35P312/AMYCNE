import argparse
import glob
import os
import common
import genotype
import annotate
import call
import time
import sys
import hist
import count

import numpy
import scipy


parser = argparse.ArgumentParser("AMYCNE a copy number estimation toolkit",add_help=False)
parser.add_argument('--genotype' , action="store_true" ,help="compute the copy number in selected regions")
parser.add_argument('--annotate' , action="store_true" ,help="add copy number estimates to structural variant VCF entries")
parser.add_argument('--call' , action="store_true" ,help="perform CNV calling")
parser.add_argument('--hist' , action="store_true" ,help="compute the coverage across each chromosome, return a tab file describing the average coverage, as well as average coverage per contig")
parser.add_argument('--count' , action="store_true" ,help="estimate the copy number of each chromosome")
parser.add_argument('--ff' , action="store_true" ,help="Predict fetal fraction based on the number of reads of chromosome Y")
parser.add_argument('--filt' , action="store_true" ,help="filters the input coverage tab file, prints the filtered and gc corrected version to stdout")
args, unknown = parser.parse_known_args()


if args.ff:

    parser = argparse.ArgumentParser("""AMYCNE-ff: Estimate fetal fraction based on read depth across chromosome Y""")
    parser.add_argument('--ff' , action="store_true" ,help="Predict fetal fraction based on the number of reads of chromosome Y")
    parser.add_argument('--coverage' , type=str, help="the tab file containing coverage tab files")
    parser.add_argument('--Q' , type=int,default=15,help="Minimum average mapping quality of the bins used for copy number estimation default = 15")
    parser.add_argument('--scaling' , type=float,default=2,help="The slope of FFY (default=2)")
    parser.add_argument('--intercept' , type=float,default=0,help="Set the intercept of FFY (default=0)")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--refQ' , type=int,default=10,help="Minimum average mapping quality of the bins used for constructing the reference = 10")
    parser.add_argument('--c_cutoff' , type=int,default=100,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy number esitmation")
    args = parser.parse_args()

    Data= common.gc_tab(args.gc)
    Data=common.coverage_tab(args.coverage,Data)
    GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
    autosomal_bins=[]
    y_bins=[]
    x_bins=[]
 
    for chromosome in Data["chromosomes"]:
        if "21" in chromosome or "13" in chromosome or "9" in chromosome or "18" in chromosome or len(chromosome) > 5:
           continue

        for i in range(0,len(Data[chromosome]["coverage"])):
           if not Data[chromosome]["GC"][i] in GC_hist or Data[chromosome]["GC"][i] == -1 or GC_hist[Data[chromosome]["GC"][i]][0] == -1:
               continue

           if Data[chromosome]["coverage"][i] > 0 and Data[chromosome]["quality"][i] < args.Q:  
               continue

           if chromosome == "Y" or chromosome == "chrY": 
               if Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0] > 0.5:
                   continue
               y_bins.append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])
           elif chromosome == "X" or chromosome == "chrX":
               x_bins.append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])

           else:
               autosomal_bins.append(Data[chromosome]["coverage"][i]/GC_hist[Data[chromosome]["GC"][i]][0])

    sex="female"
    ff="NA"

    if len(y_bins):
        medY=numpy.median(y_bins)
    else:
        medY=0
    if len(x_bins):
        medX=numpy.median(x_bins)
    else:
        medX=0

    medA=numpy.median(autosomal_bins)

    if medY > 0.005:
        sex= "male"

    print ("sample sex medAutosomal FFY FFX")
    if sex == "male":
        print ("{} {} {} {} {}".format(args.coverage.split("/")[-1].split(".")[0],sex, medA, args.scaling*medY+args.intercept, args.scaling*(1-medX) ) )
    else:
        print ("{} {} {} {} {}".format(args.coverage.split("/")[-1].split(".")[0],sex, medA, args.scaling*medY, args.scaling*(1-medX) ) )

elif args.genotype:
    parser = argparse.ArgumentParser("""AMYCNE-genotype:compute the copy number of selected regions based on the region input file""")
    parser.add_argument('--genotype' , action="store_true" ,help="compute the copy number in selected regions")
    parser.add_argument('--region' , type=str,help="the file containing selected regions")
    parser.add_argument('--R' , type=str,help="used instead of --regions, select a single target region via command line in the following format chr:start-end")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--coverage' , type=str, help="the tab file containing coverage tab files")
    parser.add_argument('--folder' , type=str,help="use every .tab file in the folder as coverage file")
    parser.add_argument('--refQ' , type=int,default=30,help="Minimum average mapping quality of the bins used for constructing the reference = 30")
    parser.add_argument('--Q' , type=int,default=10,help="Minimum average mapping quality of the bins used for copy number estimation default = 10")
    parser.add_argument('--c_cutoff' , type=int,default=100,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy nmber esitmation")
    parser.add_argument('--plody' , type=int,default=2,help="The plody of the organism")
    args = parser.parse_args()
    #get the gc content
    if not args.region and not args.R:
        print("select target regions using --R or --region")
        sys.exit(0)   
    Data= common.gc_tab(args.gc)
    
    print("sample\tbins\tused_bin_ratio\tref_coverage\tCN_raw\t95%CI\tCN_rounded\tregion_command")
    if args.coverage:
        #compute a gc content histogram
        Data=common.coverage_tab(args.coverage,Data)
        GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
        genotype.main(Data,GC_hist,args)

    elif(args.folder):
        tab_folder = glob.glob(os.path.join(args.folder,"*.tab"));
        for tab in tab_folder:
            #compute a gc content histogram
            args.coverage=tab
            Data=common.coverage_tab(args.coverage,Data)
            GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
            genotype.main(Data,GC_hist,args)
    else:
        print("coverage data is required, use either the coverage or folder option to select the input. read the manual for more info on how to generate coverage files")

elif args.annotate:
    parser = argparse.ArgumentParser("""AMYCNE-annotate:annotate the intrachromosomal variants of a structural variation vcf""")
    parser.add_argument('--annotate' , action="store_true" ,help="compute the copy number in selected regions")
    parser.add_argument('--vcf' , type=str,help="a structural variation vcf file")
    parser.add_argument('--bed' , type=str,help="a bed file containing genomic regions")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--coverage' , type=str,required= True, help="the tab file containing coverage")
    parser.add_argument('--refQ' , type=int,default=30,help="Minimum average mapping quality of the bins used for constructing the reference = 30")
    parser.add_argument('--Q' , type=int,default=10,help="Minimum average mapping quality of the bins used for copy number estimation default = 10")
    parser.add_argument('--c_cutoff' , type=int,default=100,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy nmber esitmation")
    parser.add_argument('--plody' , type=int,default=2,help="The plody of the organism")
    args = parser.parse_args()
    
    #get the gc content
    Data= common.gc_tab(args.gc)
    #compute a gc content histogram
    Data=common.coverage_tab(args.coverage,Data)
    GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
    annotate.main(Data,GC_hist,args)
    
elif args.call:
    parser = argparse.ArgumentParser("""AMYCNE-call: Detect copy number variants and print them to a vcf""")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--coverage', required=True, type=str,default=None, help="the tab file containing coverage")
    parser.add_argument('--c_cutoff' , type=int,default=100,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy number esitmation")
    parser.add_argument('--plody' , type=int,default=2,help="The plody of the organism")
    parser.add_argument('--nbins' , type=int,default=20,help="number of bins used to detect the changepoints, default = 20")
    parser.add_argument('--refQ' , type=int,default=30,help="Minimum average mapping quality of the bins used for constructing the reference = 30")
    parser.add_argument('--Q' , type=int,default=30,help="Minimum average mapping quality of the bins used for copy number estimation default = 30")
    parser.add_argument('--max' , type=int,default=6,help="Maximum ratio = 6")
    parser.add_argument('--scoref' , type=int,default=30,help="The log(P) threshold of the Fisher method test ")
    parser.add_argument('--scoren' , type=int,default=17,help="minimum non parametric phred score (default = 17)")
    parser.add_argument('--ratioLim' , type=float,default=2.5,help="Maximum variant rdr")
    parser.add_argument('--call' , action="store_true" ,help="perform CNV calling")
    parser.add_argument('--bam' ,type=str ,help="the bam file (AMYCNE will only extract the header for sample and reference information)")
    parser.add_argument('--folder' , type=str,default=None, help="a folder containing coverage files, each file will be analysed")
    parser.add_argument('--output' , type=str, help="the name of the vcf file")
    parser.add_argument('--sample' , type=str, help="the sample id(default= extracted from the SM tag of the bam field, else the  filename of the coverage file)")    
    args = parser.parse_args()

    #get the gc content    
    if args.coverage:
        if not args.output:
            args.output=args.coverage.replace(".tab",".vcf").replace(".bed",".vcf")
            if args.output == args.coverage:
                args.output=args.coverage + ".vcf"

        Data = common.gc_tab(args.gc)
        Data =common.coverage_tab(args.coverage,Data)
        #compute a gc content histogram
        GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
        call.main(Data,GC_hist,args)

elif args.hist:
    parser = argparse.ArgumentParser("""AMYCNE-genotype:compute the copy number of selected regions based on the region input file""")
    parser.add_argument('--hist' , action="store_true" ,help="compute the coverage across each chromosome, return a tab file descri bing the average coverage, as well as average coverage per contig")
    parser.add_argument('--coverage' , type=str, help="the tab file containing coverage tab files")
    parser.add_argument('--folder' , type=str,help="use every .tab file in the folder as coverage file")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    args = parser.parse_args()
    
    Data = common.gc_tab(args.gc)
    header= "sample_ID\tAverage"
    for chromosome in sorted(Data["chromosomes"]):
        if not "Un_" in chromosome and not "random" in chromosome and not "hap" in chromosome and not "M" in chromosome:
            header+= "\t" + chromosome
    print(header)
    
    if args.coverage:
        #compute a gc content histogram
        Data=common.coverage_tab(args.coverage,Data)
        hist.main(Data,args)
    elif(args.folder):
        tab_folder = glob.glob(os.path.join(args.folder,"*.tab"));
        for tab in tab_folder:
            #compute a gc content histogram
            args.coverage=tab
            Data=common.coverage_tab(args.coverage,Data)
            hist.main(Data,args)
    else:
        print("coverage data is required, use either the coverage or folder option to select the input. read the manual for more info on how to generate coverage files")

elif args.count:
    parser = argparse.ArgumentParser("""AMYCNE-genotype:compute the copy number of selected regions based on the region input file""")
    parser.add_argument('--count' , action="store_true" ,help="compute the coverage across each chromosome, return a tab file descri bing the average coverage, as well as average coverage per contig")
    parser.add_argument('--coverage' , required=True,type=str, help="the tab file containing coverage tab files")
    parser.add_argument('--ploidy' , type=int,default = 2,help="the ploidy of the organism")
    parser.add_argument('--d' , type=float,default=0.1,help="minimum ratio deviation to call chromosomal abberation")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--c_cutoff' , type=int,default=200,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy nmber esitmation")
    parser.add_argument('--refQ' , type=int,default=30,help="Minimum average mapping quality of the bins used for constructing the reference = 30")
    parser.add_argument('--Q' , type=int,default=30,help="Minimum average mapping quality of the bins used for copy number estimation default = 30")
        
    args = parser.parse_args()
    
    Data = common.gc_tab(args.gc)
    #compute a gc content histogram
    Data=common.coverage_tab(args.coverage,Data)
    GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
    count.main(Data,GC_hist,args)

elif args.filt:
    parser = argparse.ArgumentParser("""AMYCNE-filter: filter the coverage tab file, prints it to stdout for later use""")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--coverage' , type=str,required= True,default=None, help="the tab file containing coverage")
    parser.add_argument('--c_cutoff' , type=int,default=100,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy nmber esitmation")
    parser.add_argument('--filter' , type=int,default=2000,help="size of the filters, default = 2000")
    parser.add_argument('--refQ' , type=int,default=30,help="Minimum average mapping quality of the bins used for constructing the reference = 30")
    parser.add_argument('--Q' , type=int,default=10,help="Minimum average mapping quality of the bins used for copy number estimation default = 10")
    parser.add_argument('--filt' , action="store_true" ,help="perform CNV calling")

    args = parser.parse_args()

    Data = common.gc_tab(args.gc)
    Data=common.coverage_tab(args.coverage,Data)

    for chromosome in Data["chromosomes"]: 
        if Data[chromosome]["quality"]:
            print("#chromosome\tstart\tend\tcoverage\tQ")
        else:
            print("#chromosome\tstart\tend\tcoverage")
        break

    if args.filter % 2 == 0:
        args.filter += 1
    for chromosome in Data["chromosomes"]:
        median_filtered=scipy.signal.medfilt(Data[chromosome]["coverage"],args.filter)
        wiener_filter = scipy.signal.wiener(median_filtered,args.filter)
        for i in range(0,len(wiener_filter)):
            if Data[chromosome]["quality"]:
                print("{}\t{}\t{}\t{}\t{}".format(chromosome,i*Data["bin_size"],(i+1)*Data["bin_size"],round(wiener_filter[i],1),Data[chromosome]["quality"][i]))
            else:    
                print("{}\t{}\t{}\t{}".format(chromosome,i*Data["bin_size"],(i+1)*Data["bin_size"],round(wiener_filter[i],1) ))
else:
   parser.print_help()



