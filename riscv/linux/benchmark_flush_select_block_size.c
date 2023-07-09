#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <cflush.h>

int bytes_per_thread;
volatile uint64_t total_cycles = 0;

void* threadFunc(void* arg) {
    void *x = arg;

    printf("This thread got 0x%x\n", (uint64_t *)x);

    int bytes = bytes_per_thread;

    // dirty each line
    for (int i = 0; i < bytes; i+=8) {
        *((uint64_t *) (x+i)) = i; 
    }

    // //flush
    uint64_t start = read_csr(cycle);
    
    for (int i = 0; i < bytes; i += 64) {
        CBO_FLUSH_FN(x+i);
    }

    asm volatile ("fence iorw, iorw");
    
    uint64_t end = read_csr(cycle);
    uint64_t elapsed = end - start;
    total_cycles += elapsed;

    return NULL;
}

void bench(int numThreads, int bytes) {
// Initialize thread identifiers
    pthread_t* threads = malloc(numThreads * sizeof(pthread_t));

    void *x = malloc(bytes);
    bytes_per_thread = bytes/numThreads;


    printf("x starts at 0x%x\n", x);
    printf("each thread gets %d bytes\n", bytes_per_thread);

    // Create threads
    for (int i = 0; i < numThreads; i++) {
        if (pthread_create(&threads[i], NULL, threadFunc, (x + i*bytes_per_thread)) != 0) {
            printf("Error: Failed to create thread %d.\n", i);
            return;
        }
    }

    // Join threads
    for (int i = 0; i < numThreads; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            printf("Error: Failed to join thread %d.\n", i);
            return;
        }
    }

    // Clean up
    free(threads);
    free(x);
}

int compare_uint64_t(const void* a, const void* b) {
    uint64_t value1 = *(uint64_t*)a;
    uint64_t value2 = *(uint64_t*)b;

    if (value1 < value2) {
        return -1;
    } else if (value1 > value2) {
        return 1;
    } else {
        return 0;
    }
}

uint64_t find_median(uint64_t* array, int N) {
    // Sort the array in ascending order
    qsort(array, N, sizeof(uint64_t), compare_uint64_t);

    // Find the median
    if (N % 2 == 0) {
        uint64_t mid1 = array[N / 2 - 1];
        uint64_t mid2 = array[N / 2];
        return (mid1 + mid2) / 2;
    } else {
        return array[N / 2];
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: ./benchmark <num_threads> <cache_block_size> <reps for confidence>\n");
        return 1;
    }

    int numThreads = atoi(argv[1]);
    int bytes = atoi(argv[2]);
    int reps = atoi(argv[3]);
    //int sameX = atoi(argv[3]);

    uint64_t results[reps];

    for(int i=0;i<reps;i++) {
        total_cycles = 0;
        bench(numThreads, bytes);
        results[i] = total_cycles;
    }

    uint64_t median_cycles = find_median(results, reps);

    printf("Median Elapsed cycles to flush %lu bytes by %d threads: %lu\n", bytes, numThreads, median_cycles);

    return 0;
}
