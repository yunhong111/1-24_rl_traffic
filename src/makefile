CC=g++
CFLAGS = -O3     # optimize code
DFLAGS =         # common defines
OPENMP = -fopenmp

OBJ =   otherFun.o Trie_cuckoo_noise_07_21.o trie.o cuckoo_filter.o aggregation_add_cuckoo.o cuckoo_table.o hash_function.o sha1.o RL.o
PROG = trieNoiseMain

all:	$(PROG)

otherFun.o:otherFun.cpp otherFun.h
	$(CC) -c $(CFLAGS) otherFun.cpp 

Trie_cuckoo_noise_07_21.o:Trie_cuckoo_noise_07_21.cpp Trie_cuckoo_noise_07_21.h
	$(CC) -c $(OPENMP) $(CFLAGS) Trie_cuckoo_noise_07_21.cpp
 
trie.o:trie.cpp trie.h
	$(CC) -c $(CFLAGS) trie.cpp  

cuckoo_filter.o:cuckoo_filter.cpp cuckoo_filter.h
	$(CC) -c $(CFLAGS) cuckoo_filter.cpp 
 
aggregation_add_cuckoo.o:aggregation_add_cuckoo.cpp aggregation_add_cuckoo.h
	$(CC) -c $(CFLAGS) aggregation_add_cuckoo.cpp
  
cuckoo_table.o:cuckoo_table.cpp cuckoo_table.h
	$(CC) -c $(CFLAGS) cuckoo_table.cpp
  
hash_function.o:hash_function.c hash_function.h
	$(CC) -c $(CFLAGS) hash_function.c 

sha1.o:sha1.c sha1.h
	$(CC) -c $(CFLAGS) sha1.c

RL.o:RL.cc RL.h
	$(CC) -c $(CFLASS) RL.cc

$(PROG): $(OBJ)
	$(CC) $(CFLAGS) -fopenmp -o $@ $^

.PHONY: clean
clean:
	rm -f *.o $(PROG)



