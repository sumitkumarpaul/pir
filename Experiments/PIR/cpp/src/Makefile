PREFIX = /usr/local
dbg = -g
CC=g++
#CFLAGS=-fomit-frame-pointer -maes -std=c++11 -msse4.2 -fno-strict-aliasing -pedantic -Wall -Wextra -Wunreachable-code -Wmissing-declarations -Wunused-function -Wno-overlength-strings -Wno-deprecated-declarations -O3 -s -DNDEBUG -I../libtfhe -I../libfss -I/usr/local/include/openfhe/core -I/usr/local/include/openfhe/binfhe -I/usr/local/include/openfhe -I/usr/local/include/openfhe/pke
CFLAGS=-fomit-frame-pointer -mavx -mavx2 -mavx512f -maes -std=c++17 -msse4.2 -fno-strict-aliasing -pedantic -O3 -s -DNDEBUG -DAESNI -I../libfss -I/usr/local/include/openfhe/core -I/usr/local/include/openfhe/binfhe -I/usr/local/include/openfhe -I/usr/local/include/openfhe/pke -fopenmp -fPIC -pthread -m64 -Wa,--noexecstack -Wall -Wa,--noexecstack -g -fdebug-prefix-map=/build/openssl-zfps2k/openssl-1.1.1f=. -fstack-protector-strong -Wformat -Werror=format-security -DOPENSSL_TLS_SECURITY_LEVEL=2 -DOPENSSL_USE_NODELETE -DL_ENDIAN -DOPENSSL_PIC -DOPENSSL_CPUID_OBJ -DOPENSSL_IA32_SSE2 -DOPENSSL_BN_ASM_MONT -DOPENSSL_BN_ASM_MONT5 -DOPENSSL_BN_ASM_GF2m -DSHA1_ASM -DSHA256_ASM -DSHA512_ASM -DKECCAK1600_ASM -DRC4_ASM -DMD5_ASM -DAESNI_ASM -DVPAES_ASM -DGHASH_ASM -DECP_NISTZ256_ASM -DX25519_ASM -DPOLY1305_ASM -DNDEBUG -Wdate-time -D_FORTIFY_SOURCE=2
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
