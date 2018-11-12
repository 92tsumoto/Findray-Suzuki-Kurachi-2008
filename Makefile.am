#INCLUDES = syspara.h

bin_PROGRAMS = FSK 
FSK_CFLAGS = -O2 -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg
#FSK_CFLAGS = -O2 -axsse4.2 -parallel -Wall

FSK_INCLUDE = -I./include
FSK_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ./lib/xhplot.o

FSK_SOURCES = ./src/fsk.c\
	./src/ion_concent1.c ./src/ion_concent2.c ./src/win_out.c\
	./src/eular.c ./src/input_para.c ./src/system.c\
	./src/eventloop.c ./src/exptable.c ./src/syspara.h

