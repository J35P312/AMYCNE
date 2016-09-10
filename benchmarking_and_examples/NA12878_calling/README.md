AMYCNE variant calling
extract and copy the NA12878 tab file and the GC tab file from the delly_filtering folder.
Run AMYCNE variant calling according the manual. In the AMYCNE paper, AMYCNE was run using default settings

To benchmark AMYCNE against CNVnator, install and run cnvnator according to the CNVnator manual. Or use one of the cnvnator vcf files.
Compute precision and specificity using the benchmarking_NA.sh script:
./benchmarking_NA.sh input_VCF.vcf
