The following folder contains instructions on how to generate a coverage tab file using TIDDIT.

	1 download TIDDIT via this link:
		git clone https://github.com/SciLifeLab/TIDDIT.git

	2 compile according to the instructions:
		cd TIDDIT
		mkdir build
		cd build
		cmake ..
		make

	3 run coverage analysis:
		cd ..
		cd bin
		./TIDDIT --cov -b /path/to/bam/file.bam --bin_size x

	where x is the size of the bins, for variant calling x = 1000 is recomended. For genotyping and annotation the corect bin size depends on the size
	of the regions to be annotated/genotyped.

	
