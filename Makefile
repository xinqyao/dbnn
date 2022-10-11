LIBS	     = -lm
CC           = cc 
LD           = $(CC)

all: nn dbndatagen nndatagen combine charpred scripts 
#nn nndatagen dbndatagen combine charpred scripts
combine:	combine.o tools.o
	$(LD) -o $@ combine.o tools.o $(LIBS) 
charpred:	charpred.o tools.o
	$(LD) -o $@ charpred.o tools.o $(LIBS)
nndatagen:	nndatagen.o tools.o
	$(LD) -o $@ nndatagen.o tools.o $(LIBS)
nn:		nn.o tools.o
	$(LD) -o $@ nn.o tools.o $(LIBS)
dbndatagen:	dbndatagen.o tools.o
	$(LD) -o $@ dbndatagen.o tools.o $(LIBS)	
scripts:
	chmod 777 dbntrain dbnpred nntrain nnpred rundbnn
clean:
	rm -f *.o
	rm -f dbndatagen nndatagen nn combine charpred
	rm -f ~/.dbnn.check.ok .dbn.ok 
