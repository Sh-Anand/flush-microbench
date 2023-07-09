#include <stdio.h>
#include <cflush.h>
#include <stdlib.h>

int bytes_per_thread;
volatile uint64_t total_cycles = 0;

#define U32 *(volatile unsigned int *)
#define MEMORY_MEM_ADDR 0x80000000UL

void bench(int cid, void *x) {
    printf("Core %d got 0x%x\n", cid , x);

    int bytes = bytes_per_thread;

    // dirty each line
    for (int i = 0; i < bytes; i+=8) {
        *((uint64_t *) (x+i)) = i; 
    }

    // //flush
    uint64_t start = read_csr(cycle);
    
    for (int i = 0; i < bytes; i += 64) {
        printf("Core %d flushing 0x%x\n", cid, x+i);
        CBO_CLEAN_FN(x+i);
    }

    asm volatile ("fence iorw, iorw");
    
    uint64_t end = read_csr(cycle);
    uint64_t elapsed = end - start;
    printf("Elapsed core %d, %lu cycles\n", cid, elapsed);

}

void thread_entry(int cid, int nc) {

    printf("threads, bytes, mean, median, stddev\n");
    int reps = 10;
    int numThreads = 2;
    for(int bytes=128;bytes<=16384;bytes*=2) {
        bytes_per_thread = bytes/2;
        printf("Core %d handling %d long flushes\n", cid, bytes);
        for(int i=0;i<reps;i++) {
            total_cycles = 0;
            if(cid == 0)
                bench(cid, (void *)MEMORY_MEM_ADDR);
            else
                bench(cid, (void *)(MEMORY_MEM_ADDR+(bytes/2)));
        }
        printf("\n");

    }

    exit(0);
}