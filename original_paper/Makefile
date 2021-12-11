SRC= anf.c
AUX= README makefile sample.inp sample.out

demo: ./anf sample.inp
	./anf < sample.inp > out
	echo "comparing result with older result - should be similar "
	echo "     but NOT necessarily identical"
	paste out sample.out
	echo "end of comparison"

anf: anf.c
	gcc -O anf.c -o anf -lpthread -lm

anf.tar: $(SRC)  $(AUX)
	tar cvf anf.tar $(SRC) $(H) $(AUX)

uu: anf.tar
	gzip anf.tar
	uuencode anf.tar.gz anf.tar.gz > uu

clean:
	\rm -f *.o *~ *.bak anf.tar uu out anf.tar.gz

spotless: clean
	\rm -f anf
