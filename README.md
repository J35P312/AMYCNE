# AMYCNE

AMYCNE is a copy number estimation toolkit, designed for WGS data. It contains modules for genotyping of copy number, counting the number of chromosomes, annotating vcf files, and calling CNVs. AMYCNE require coverage tab files as input, these files may e produced using TIDDIT.
 
# Installation

AMYCNE requires bottleneck, scipy and numpy. These two packages may be installed using pip:

	pip install numpy
	pip install scipy
	pip install bottleneck

AMYCNE has been tested on python 2.7.11, but might run on older versions of python as well.
To improve the performance of AMYCNE, the code of AMYCNE may be compiled using cython:

python setup.py build_ext --inplace


# Run
Type the following for a list of modules:

    python AMYCNE.py

The following section decribes the basic commands for running AMYCNE, for more info, use the --help flag for each module.
Each module requires a coverage file and a gc content file, having the same bin size. 

    Genotype: estimate the copy number in one or more target region
      Use the genotype module by typing:
        python AMYCNE.py --genotype
        
      Use the following command to genotype a specified region(chr1:100-10000):
      python AMYCNE.py --genotype --gc GC_tab_file.tab --coverage coverage_tab_file.tab --R chr1:100-10000
      
      To genotype all coverage files in a folder type the followin command:
      python AMYCNE.py --genotype --gc GC_tab_file.tab --folder /path/to/folder --R chr1:100-10000
      
      Multiple regions could be genotyped using a region text file instead of the --R flag:
      python AMYCNE.py --genotype --gc GC_tab_file.tab --coverage coverage_tab_file.tab --region region.txt
      python AMYCNE.py --genotype --gc GC_tab_file.tab --folder /path/to/folder --region region.txt
      
    The region file consists of operations. Each line within the region text file describes one operation. THe supported operations are sum(sum) and average(avg). The operations are written in the following format:
    
    sum(1:104198143-104207173|1:104230040-104238912|1:104292279-104301311)
    avg(1:104198143-104207173|1:104230040-104238912|1:104292279-104301311)
    
    Each region is separated by |, any number of regions within one operation is supported(except 0), and any number of operations within a region file is allowed.
    
    Anotate: estimate the copy number across structural variants stored in a vcf file
        Use the anotate module by typing:
          python AMYCNE.py --anotate
          
        The annotate module requires a coverage file, a gc content file, as well as the structural variant vcf:
         python AMYCNE.py --anotate --gc gc_content_file.tab --coverage coverage.tab --vcf sv.vcf > annotated.sv.vcf
    
    Call: Perform CNV calling
      Use the call module by typing:
          python AMYCNE.py --call
          
      The call module requires a coverage file, and a gc content file, and prints vcf file to --output
         python AMYCNE.py --call --gc gc_content_file.tab --coverage coverage.tab --output out.vcf
    
# Generate GC content file

  The Generate_GC_tab.py script may used to generate gc content files:
  python Generate_GC_tab.py --fa reference.fa --size bin_size > gc_content.tab

  note that AMYCNE requires the same bin size for the coverage

# Generate coverage Files

Coverage tab files may be generated using TIDDIT.
The files should be given in the following format:

#chromosome	start	end	coverage	quality

chr1	0	100	23	10

chr1	100	200	23	10

chr1	200	300	23	10

chr2	0	100	23	10

chr2	100	200	23	10

chr2	200	300	23	10

chrX	0	100	23	10

chrY	0	100	23	10

The quality column is optional, the header is also optional, and is not read by the software. The bins need to cover the entire genome.
