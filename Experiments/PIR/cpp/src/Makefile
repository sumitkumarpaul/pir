PREFIX = /usr/local
dbg = -g
CC=g++
#CFLAGS=-fomit-frame-pointer -maes -std=c++11 -msse4.2 -fno-strict-aliasing -pedantic -Wall -Wextra -Wunreachable-code -Wmissing-declarations -Wunused-function -Wno-overlength-strings -Wno-deprecated-declarations -O3 -s -DNDEBUG -I../libtfhe -I../libfss -I/usr/local/include/openfhe/core -I/usr/local/include/openfhe/binfhe -I/usr/local/include/openfhe -I/usr/local/include/openfhe/pke
CFLAGS=-fomit-frame-pointer -maes -std=c++11 -msse4.2 -fno-strict-aliasing -pedantic -O3 -s -DNDEBUG -DAESNI -I../libfss -I/usr/local/include/openfhe/core -I/usr/local/include/openfhe/binfhe -I/usr/local/include/openfhe -I/usr/local/include/openfhe/pke -fopenmp
LDFLAGS= -L/usr/lib/x86_64-linux-gnu -Wl,-Bstatic -lcrypto -Wl,-Bdynamic -ldl -lssl -lgmpxx -lgmp -lpthread ../libfss/libfss.a -L/usr/local/lib -lOPENFHEpke -lOPENFHEcore -fopenmp
SOURCES= fss-client.cpp fss-server.cpp fss-common.cpp pir-test.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=pir-test
LIB=libfss.a
ifeq ($(AESNI), 1)
    CFLAGS+=-DAESNI
endif

all: $(SOURCES) $(EXECUTABLE) $(LIB)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(LIB): fss-client.o fss-server.o fss-common.o
	ar rcs $@ $^

.cpp.o:
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: clean
clean:
	rm -f *.o $(LIB) $(OBJECTS) $(EXECUTABLE)

.PHONY: install
install: 
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	mkdir -p $(DESTDIR)$(PREFIX)/include
	cp $(alib) $(DESTDIR)$(PREFIX)/lib/$(alib)
	cp *.h $(DESTDIR)$(PREFIX)/include/
