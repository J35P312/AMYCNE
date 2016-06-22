import argparse
import glob
import os
import common
import genotype

parser = argparse.ArgumentParser("AMYCNE a copy number estimation toolkit",add_help=False)
parser.add_argument('--genotype' , action="store_true" ,help="compute the copy number in selected regions")
parser.add_argument('--annotate' , action="store_true" ,help="add copy number estimates to structural variant VCF entries")
parser.add_argument('--call' , action="store_true" ,help="perform variant calling(not implemented)")
args, unknown = parser.parse_known_args()


if args.genotype:

    parser = argparse.ArgumentParser("""AMYCNE-genotype:compute the copy number of selected regions based on the region input file""")
    parser.add_argument('--genotype' , action="store_true" ,help="compute the copy number in selected regions")
    parser.add_argument('--region' , type=str, required=True,help="the selected regions")
    parser.add_argument('--gc' , type=str,required= True, help="the tab file containing gc content")
    parser.add_argument('--coverage' , type=str, help="the tab file containing coverage")
    parser.add_argument('--folder' , type=str,help="use every .tab file in the folder as coverage file")
    args = parser.parse_args()

    #get the gc content
    GC= common.gc_tab(args.gc)

    
    print("sample\tregional_GC\tcoverage_ratio\tref_coverage\tCN")
    if args.coverage:
        #compute a gc content histogram
        Data=common.coverage_tab(args.coverage,GC)
        GC_hist=common.gc_hist(Data)
        genotype.main(Data,GC_hist,args)
    elif(args.folder):
        tab_folder = glob.glob(os.path.join(args.folder,"*.tab"));
        for tab in tab_folder:
            #compute a gc content histogram
            args.coverage=tab
            Data=common.coverage_tab(args.coverage,GC)
            GC_hist=common.gc_hist(Data)
            genotype.main(Data,GC_hist,args)
    else:
        print("coverage data is required, use either the coverage or folder option to select the input. read the manual for more info on how to generate coverage files")

elif args.annotate:
    pass
elif args.call:
    pass
else:
    parser.print_help()
