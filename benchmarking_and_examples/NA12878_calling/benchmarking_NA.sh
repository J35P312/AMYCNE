#first argument: the vcf path
#second argument: the name of the caller
echo "the number of validated deletions:"
wc -l 4000_tab
echo "the number of hits:"
python query/query_db.py --variations $1 --files 4000_tab > $2_query.vcf
grep "OCC=1" $2_query.vcf | wc -l
echo "the number of false call deletions(larger han 500bp):"
python VCF_intra_size_filter.py --vcf $2_query.vcf > $2_size_filt.vcf
grep "DEL" $2_size_filt.vcf > $2_size_filt_deletions.vcf
grep "OCC=0" $2_size_filt_deletions.vcf | wc -l

echo " the number of hits after using the filter tag as filter:"
python cleanVCF.py --vcf $2_query.vcf > $2_qual.vcf
grep "OCC=1" $2_qual.vcf | wc -l
echo "the number of false called deletions after filtering:"
python VCF_intra_size_filter.py --vcf $2_qual.vcf > $2_size_filt_qual.vcf
grep "DEL" $2_size_filt_qual.vcf > $2_size_filt_deletions_qual.vcf
grep "OCC=0" $2_size_filt_deletions_qual.vcf | wc -l
