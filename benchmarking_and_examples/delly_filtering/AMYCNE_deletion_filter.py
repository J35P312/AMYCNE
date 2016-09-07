import sys

for line in open(sys.argv[1]):
    if not line[0] == "#":
        content=line.split("\t")
        info=content[7]
        CNV=info.split("AMYCNE=")[-1].split(";")[0]
        CI=info.split("AMYCNECI=")[-1].split(";")[0]
        Bins=float(info.split("BIN_RATIO=")[-1].split(";")[0])
        CI=CI.split(",")
        CI_test=float(CI[0]) <= 1.75 and float(CI[1]) <= 2.25
        if int(CNV) < 2 or CI_test:
           if Bins > 0.5:
              print line.strip()
    else:
        print(line.strip())
