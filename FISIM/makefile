
INCLUDE = /usr/include/python2.7
CC = gcc

all: ext_clustering.so

ext_clustering.o: ext_clustering.c
	$(CC) -fPIC -I. -I$/home/kipkurui/anaconda2/include/python2.7 -I/home/kipkurui/anaconda2/include/python2.7 -c ext_clustering.c -o ext_clustering.o

ext_clustering.so: ext_clustering.o
	$(CC)  -shared -o ext_clustering.so ext_clustering.o

clean:
	rm *.o *.so
