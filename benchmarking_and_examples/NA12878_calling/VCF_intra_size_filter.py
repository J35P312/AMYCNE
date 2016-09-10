import sys
import argparse

parser = argparse.ArgumentParser("""removes all vcf entries that are not marked pass or .""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
args, unknown = parser.parse_known_args()
description={}
for line in open(args.vcf):
    
    if(line[0] == "#"):
        print(line.strip())
    else:
        try:
            if "END=" in line:
                content=line.split("\t")
                
                INFO=content[7].split(";");
                for tag in INFO:
                    tag=tag.split("=")
                    if(len(tag) > 1):
                        description[tag[0]]=tag[1];
                
                END=int(description["END"])
                if( abs( END-int(content[1]) ) >= 4000):
                    print(line.strip())
        except:
            pass
