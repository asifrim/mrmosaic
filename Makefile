default:
	make cnv_baf_data
	make install_r_package

cnv_baf_data: 
	wget -N https://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
	tar -vjxf samtools-1.3.tar.bz2
	cp cnv_baf_data.c samtools-1.3
	cd samtools-1.3; make clean; make;
	cd samtools-1.3; gcc -O3 -Wall -Isrc -I. -Ihtslib-1.3/ -Ihtslib-1.3/htslib -rdynamic -o cnv_baf_data -Lhtslib-1.3/ -L. cnv_baf_data.c -lhts -lbam -lpthread -lz -lm
	mkdir -p bin
	mv samtools-1.3/cnv_baf_data bin

install_r_package:
	Rscript install_package.R

clean:
	rm samtools-1.3.tar.bz2 
	rm -rf samtools-1.3
	rm -rf bin