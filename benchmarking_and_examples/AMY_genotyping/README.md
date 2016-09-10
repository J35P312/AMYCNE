AMYCNE genotyping
this folder contains a few examples on how to comput the copy number of AMY1, AMY2A and AMY2B using AMYCNE. 
you vould either compute the tab files and gc content file on your own, or use the NA12878 tab files in the delly_filtering
example folder.

the following command will compute the copy number of AMY1
	python AMYCNE.py --genotype --gc GC_tab_file.tab --coverage coverage_tab_file.tab --region AMY1.txt --Q 0
	
the following command will compute the copy number of AMY2A
	python AMYCNE.py --genotype --gc GC_tab_file.tab --coverage coverage_tab_file.tab --region AMY2A.txt --Q 0
	
the following command will compute the copy number of AMY2B
	python AMYCNE.py --genotype --gc GC_tab_file.tab --coverage coverage_tab_file.tab --region AMY2B.txt --Q 0
