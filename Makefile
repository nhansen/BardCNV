INC = -I printCompCounts

subdirs ::
	cd c && make all
	mkdir -p blib/bin && cp c/bardcnv blib/bin/
