import common
import readVCF

def main(Data,GC_hist,args):
    args.plody=2

    for line in open(args.vcf):
        if line[0] == "#" and line[1] == "#":
            print line.strip()
        elif line[0] == "#":
            print "##INFO=<ID={},Number=1,Type=Float,Description=\"estimated copy number\">".format("AMYCNE")
            print line.strip()
        else:
            
            chrA,posA,chrB,posB,event_type,INFO,FORMAT = readVCF.readVCFLine(line)
            if chrA == chrB:
                if chrA in Data:
                    cn, gc, length,ref =common.regional_cn_est( Data ,GC_hist, [chrA,posA,posB] )
                     
                    CNE = round( cn*args.plody )
                    content=line.split("\t")
                    content[7] += ";AMYCNE=" + str(int(CNE))
                    new_line="\t".join(content)
                    print new_line.strip()
            else:
                print line.strip()
