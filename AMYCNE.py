import argparse
import glob
import os
import common
import genotype
import annotate
import call
import time
import sys

parser = argparse.ArgumentParser("AMYCNE a copy number estimation toolkit",add_help=False)
parser.add_argument('--genotype' , action="store_true" ,help="compute the copy number in selected regions")
parser.add_argument('--annotate' , action="store_true" ,help="add copy number estimates to structural variant VCF entries")
parser.add_argument('--call' , action="store_true" ,help="perform CNV calling")
args, unknown = parser.parse_known_args()


if args.genotype:

    parser = argparse.ArgumentParser("""AMYCNE-genotype:compute the copy number of selected regions based on the region input file""")
    parser.add_argument('--genotype' , action="store_true" ,help="compute the copy number in selected regions")
    parser.add_argument('--region' , type=str,help="the file containing selected regions")
    parser.add_argument('--R' , type=str,help="used instead of --regions, select a single target region via command line in the following format chr:start-end")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--coverage' , type=str, help="the tab file containing coverage")
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
    parser.add_argument('--vcf' , type=str, required=True,help="a structural variation vcf file")
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

    parser.add_argument('--coverage' , type=str,required= True, help="the tab file containing coverage")
    parser.add_argument('--c_cutoff' , type=int,default=100,help="bins having coverage higher than the cut off value are excluded from the ref calculations")
    parser.add_argument('--s_cutoff' , type=int,default=50,help="bins that have less than the s_cutoff value similar bins are discarded from copy nmber esitmation")
    parser.add_argument('--plody' , type=int,default=2,help="The plody of the organism")
    parser.add_argument('--min_var' , type=int,default=2000,help="smallest variant, given in bases, default = 2000")
    parser.add_argument('--filter' , type=int,default=2000,help="size of the filters, default = 2000")
    parser.add_argument('--refQ' , type=int,default=30,help="Minimum average mapping quality of the bins used for constructing the reference = 30")
    parser.add_argument('--Q' , type=int,default=10,help="Minimum average mapping quality of the bins used for copy number estimation default = 10")
    parser.add_argument('--call' , action="store_true" ,help="perform CNV calling")
    args = parser.parse_args()

    #get the gc content
    Data = common.gc_tab(args.gc)
    Data =common.coverage_tab(args.coverage,Data)
    #compute a gc content histogram
    GC_hist=common.gc_hist(Data,args.c_cutoff,args.s_cutoff,args.refQ)
    call.main(Data,GC_hist,args)
    


else:
    parser.print_help()
