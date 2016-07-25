import common
import readVCF
import numpy

def main(Data,GC_hist,args):

    for line in open(args.vcf):
        if line[0] == "#" and line[1] == "#":
            print line.strip()
        elif line[0] == "#":
            print "##INFO=<ID={},Number=1,Type=Float,Description=\"estimated copy number\">".format("AMYCNE")
            print "##INFO=<ID={},Number=2,Type=Float,Description=\"95% confidence interval around the estimated CN\">".format("AMYCNECI")
            print "##INFO=<ID={},Number=1,Type=Float,Description=\"ratio of bins used for CN estimation\">".format("BIN_RATIO")
            print "##INFO=<ID={},Number=1,Type=Float,Description=\"mean coverage of the reference bins\">".format("REFCOV")
            print line.strip()
        else:
            
            chrA,posA,chrB,posB,event_type,INFO,FORMAT = readVCF.readVCFLine(line)
            if chrA == chrB:
                if chrA in Data:
                    cn, gc, length,ref,bins,used_bins,bin_list =common.regional_cn_est( Data ,GC_hist, [chrA,posA,posB],args.Q )
                    ci="(0,0)"
                    CNE = round( cn*args.plody )
                    if CNE < 0:
                        CNE = -1
                    else:    
                        for i in range(0,len(bin_list)):
                            bin_list[i]=bin_list[i]*args.plody
                        SEM=(numpy.std(bin_list)/numpy.sqrt(used_bins))
                        ci="{},{}".format(round(cn*args.plody-SEM*1.96,2),round(cn*args.plody+SEM*1.96,2))
                        
                    content=line.split("\t")
                    if bins == 0:
                        bins =-1
                    content[7] += ";AMYCNE=" + str(int(CNE)) + ";BIN_RATIO=" + str(used_bins/float(bins)) + ";REFCOV=" + str(round(ref,2)) + ";AMYCNECI=" + ci
                    new_line="\t".join(content)
                    print new_line.strip()
            else:
                print line.strip()
