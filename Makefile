LIBS    = -lm -lefence
BINDIR  = ./bin
#CFLAGS  = -Wall -Wshadow -g
CFLAGS  = -O -Wall -Wshadow -g 
#CFLAGS  = -fast -inline speed
#CFLAGS  = -O3 -arch host -fast
CC        = gcc $(CFLAGS)
#CC        = cc $(CFLAGS)
DCC     = gcc -g -Wall -DDMALLOC_FUNC_CHECK -ansi -pedantic
PLUSCC    = g++ $(CFLAGS)
LIBS   = -lm -L/usr/local/lib
#DLIBS   = -lm -L/usr/local/lib -ldmalloc
CGCC    = checkergcc -g -Wall
CGLIBS  = -lm
# -O
PROGS   =           recoal

recoal : recombine.o jdrop.o rec_modellike.o jworld.o getdata.o \
	traitlike.o getmsatdata.o
	$(CC) -o recoal recombine.o jdrop.o rec_modellike.o \
	jworld.o getdata.o traitlike.o getmsatdata.o $(LIBS)

jworld.o : jworld.c
	$(CC) -c jworld.c

jdrop.o : jdrop.c
	$(CC) -c jdrop.c

rec_modellike.o : rec_modellike.c
	$(CC) -c rec_modellike.c

recombine.o : recombine.c
	$(CC) -c recombine.c 

getdata.o : getdata.c
	$(CC) -c getdata.c

getmsatdata.o : getmsatdata.c
	$(CC) -c getmsatdata.c

traitlike.o : traitlike.c
	$(CC) -c traitlike.c

rectreedna : rectreedna.c
	$(CC) -o rectreedna rectreedna.c $(LIBS)

hapdna : hapdna.c
	$(CC) -o hapdna hapdna.c $(LIBS)

mhapdna : mhapdna.c
	$(CC) -o mhapdna mhapdna.c $(LIBS)

drecombine.o : recombine.c
	$(CC) -o drecombine.o -c recombine.c

djdrop.o : jdrop.c
	$(CC)  -o djdrop.o -c jdrop.c

drec_modellike.o : rec_modellike.c
	$(CC) -o drec_modellike.o -c rec_modellike.c

djworld.o : jworld.c
	$(CC) -o djworld.o -c jworld.c

dgetdata.o : getdata.c
	$(CC) -o dgetdata.o -c getdata.c

dgetmsatdata.o : getmsatdata.c
	$(CC) -o dgetmsatdata.o -c getmsatdata.c

drecombine:  drecombine.o djdrop.o drec_modellike.o djworld.o dgetdata.o memdebug.o memalpha.o memfree.o \
        dgetmsatdata.o
	$(CC) drecombine.o djdrop.o drec_modellike.o djworld.o dgetdata.o memdebug.o memalpha.o memfree.o \
	dgetmsatdata.o $(LIBS) -o drecombine

dclean: 
	rm drecombine.o djdrop.o drec_modellike.o djworld.o dgetdata.o dgetmsatdata.o

clean :
	rm -f *.o

gentrees.o : gentrees.c
	$(CC) -c gentrees.c

segtre_mig.o : segtre_mig.c
	$(CC) -c segtre_mig.c

rectree : gentrees.o segtre_mig.o
	$(CC) gentrees.o segtre_mig.o $(LIBS) -o rectree

jrectree : simrectree.c
	$(CC) simrectree.c $(LIBS) -o jrectree

snp : snp.c
	$(CC) snp.c $(LIBS) -o snp

panelmaker : panelmaker.c
	$(CC) panelmaker.c $(LIBS) -o panelmaker

variance : variance.c
	$(CC) variance.c $(LIBS) -o variance

distree : distree.c
	$(CC) distree.c $(LIBS) -o distree
