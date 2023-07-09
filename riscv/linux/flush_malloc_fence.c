
#include "cflush.h"
#include <stdlib.h>
#include <stdio.h>

#define U32 *(volatile unsigned int *)

int main() {

    void *x = malloc(2048);
    for(unsigned int i=0;i<256;++i) {
        U32(x + i*4) = i*4;
    }

    unsigned start = read_csr(cycle);
    asm volatile("fence iorw, iorw");
    CBO_FLUSH_FN(x);
    CBO_FLUSH_FN(x + 0x40);

    unsigned store_cycles = read_csr(cycle);
    asm volatile("fence iorw, iorw");
    U32(x) = 97;
    U32(x + 0x40) = 98;

    asm volatile("fence iorw, iorw");
    unsigned partial = read_csr(cycle);
    asm volatile("fence iorw, iorw");

    CBO_FLUSH_FN(x);
    CBO_FLUSH_FN(x + 0x40);
    unsigned end = read_csr(cycle);
    asm volatile("fence iorw, iorw");

    printf("Start %u store %u partial %u end %u\n", start, store_cycles, partial, end);
    printf("Full flush cycles %u\n", end - start);
    printf("Partial flush cycles %u\n", partial - start);
    printf("Full flush cycles %u\n", partial - store_cycles);

    return 0;
}
