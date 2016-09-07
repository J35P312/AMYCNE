This folder contains instructions on how to annotate a vcf file using AMYCNE. And gives an example on how to use this annotation as a CNV filter.
First download Delly and the NA12878 bam, and call deletions using delly. In the AMYCNE article, delly v0.7.2 was run using default settings:

delly -t DEL -o $3.DEL.vcf -g $2 $1

Thereafter, compute a tab coverage file using TIDDIT

Then annotate the Delly vcf(NA12878.DEL.vcf) using AMYCNE, and use the AMYCNE_deletion_filter.py.
python AMYCNE --annotate --vcf NA12878.DEL.vcf --coverage NA12878_100bp.tab.tab --gc 100_GC_hg19.tab > annotated_deletions.vcf

The deletion filter will remove all deletion calls whose copy number is either higher than one, or whose confidence interval clearly is not
separated from 1 copy.
 
python AMYCNE_deletion_filter.py annotated_deletions > cnv_filtered_del.vcf

Validated deletions of the NA12878 sample is stored in this file NA12878_500.db. These deletions are equal to or larger than 500 bases.
A similar file is made by downloading the following fileftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.xlsx, rearanging the columns to the same format as the NA12878_500.db, saving the resulting
file as a tab delimited file and lastly removing all variants that are shorter/longer than a certain size.

Lastly, compute the precision and specificity of delly using the delly filter, the test filter() and the two filters combined:

First run the non-filtered vcf 
./benchmarking_NA.sh NA12878.DEL.vcf

then run the filtered vcf(where all seemingly non-deletion variatns are removed)
./benchmarking_NA.sh cnv_filtered_del.vcf

Note that the two filters are similar, and that the maximum precision is attained by combining the filters.


