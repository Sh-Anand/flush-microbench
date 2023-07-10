#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <math.h>
#include <sched.h>

#define WARMUP_BYTES 4096
#define WARMUP_N 1000
#define BURN 1000000

#define FREQ 2000000000
#define MUL 20

int bytes_per_thread;
volatile uint64_t total_cycles = 0;

static inline uint64_t read_hw_counter() {
    uint64_t t;
    asm volatile("mrs %0, cntvct_el0" : "=r" (t));
    return t;
}

static inline uint64_t read_hw_counter_freq() {
    uint64_t f;
    asm volatile("mrs %0, cntfrq_el0" : "=r" (f));
    return f;
}

void* threadFuncclflush(void* arg) {

    int my_cpu = *(int *)arg;
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(my_cpu, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);

    int burn = 0;
    while(burn++ < BURN);

    int bytes = bytes_per_thread;

    void *x = malloc(bytes);

    // dirty each line
    for (int i = 0; i < bytes; i+=8) {
        *((uint64_t *) (x+i)) = i; 
    }

    asm volatile ("dsb sy" ::: "memory");
    // //flush
    uint64_t start = read_hw_counter();
    
    for (int i = 0; i < bytes; i += 64) {
        asm volatile (
        "dc civac, %[address]\n\t"
        :
        : [address] "r" (x+i)
        : "memory"
        );
    }
    asm volatile ("dsb sy" ::: "memory");
    
    uint64_t end = read_hw_counter();
    uint64_t elapsed = end - start;
    total_cycles += elapsed*MUL;

    free(x);

    return NULL;
}

void* threadFuncclflushopt(void* arg) {

    int my_cpu = *(int *)arg;
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(my_cpu, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);

    int burn = 0;
    while(burn++ < BURN);

    int bytes = bytes_per_thread;

    void *x = malloc(bytes);

    // dirty each line
    for (int i = 0; i < bytes; i+=8) {
        *((uint64_t *) (x+i)) = i; 
    }

    asm volatile ("dsb sy" ::: "memory");
    // //flush
    uint64_t start = read_hw_counter();
    
    for (int i = 0; i < bytes; i += 64) {
        asm volatile (
        "dc cvac, %[address]\n\t"
        :
        : [address] "r" (x+i)
        : "memory"
        );
    }

    asm volatile ("dsb sy" ::: "memory");
    
    uint64_t end = read_hw_counter();
    uint64_t elapsed = end - start;
    total_cycles += elapsed*MUL;
    
    free(x);

    return NULL;
}



void warmup(int type) {
    void *x = malloc(WARMUP_BYTES);

    // dirty each line
    for (int i = 0; i < WARMUP_BYTES; i+=8) {
        *((uint64_t *) (x+i)) = i; 
    }

    asm volatile ("dsb sy" ::: "memory");
    // //flush

    for (int i = 0; i < WARMUP_BYTES; i += 64) {
        if(type == 0)
            asm volatile (
            "dc civac, %[address]\n\t"
            :
            : [address] "r" (x+i)
            : "memory"
            );
        else if(type == 1)
            asm volatile (
            "dc cvac, %[address]\n\t"
            :
            : [address] "r" (x+i)
            : "memory"
            );
    }

    asm volatile ("dsb sy" ::: "memory");
    

    free(x);

}


void bench(int numThreads, int bytes, void *(*func)(void *)) {
// Initialize thread identifiers
    pthread_t* threads = malloc(numThreads * sizeof(pthread_t));

    bytes_per_thread = bytes/numThreads;

    // Create threads
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);

    for (int i = 0; i < numThreads; i++) {
        int core = i%2;
        CPU_ZERO(&cpuset);
        CPU_SET(core, &cpuset);
        if (pthread_create(&threads[i], NULL, func, &core) != 0) {
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

double calculate_mean(uint64_t* array, int N) {
    uint64_t sum = 0;

    for (int i = 0; i < N; i++) {
        sum += array[i];
    }

    return (double)sum / N;
}

double calculate_median(uint64_t* array, int N) {
    qsort(array, N, sizeof(uint64_t), compare_uint64_t);

    if (N % 2 == 0) {
        uint64_t mid1 = array[N / 2 - 1];
        uint64_t mid2 = array[N / 2];
        return (double)(mid1 + mid2) / 2;
    } else {
        return (double)array[N / 2];
    }
}

double calculate_standard_deviation(uint64_t* array, int N, double mean) {
    double sum_of_squares = 0;

    for (int i = 0; i < N; i++) {
        double diff = array[i] - mean;
        sum_of_squares += diff * diff;
    }

    double variance = sum_of_squares / N;
    return sqrt(variance);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: ./benchmark <reps for confidence>\n");
        return 1;
    }

    int flush_types = 2;
    char *flush_array[] = {"dccivac", "dccvac"};
    void *(*flush_func_array[])(void *) = { threadFuncclflush, threadFuncclflushopt };

    int reps = atoi(argv[1]);
    //int sameX = atoi(argv[3]);

    for(int type=0;type<flush_types;type++) {
        printf("%s\n", flush_array[type]);
        printf("threads, bytes, mean, median, stddev\n");

        for (int numThreads=1;numThreads<=8;numThreads*=2) {
            
            int starter_bytes = numThreads*64;
            for (int z=0;z<WARMUP_N;z++) {
                warmup(type);
            }

            for(int bytes=starter_bytes;bytes<=16384;bytes*=2) {

                uint64_t results[reps];

                for(int i=0;i<reps;i++) {
                    total_cycles = 0;
                    bench(numThreads, bytes, flush_func_array[type]);
                    results[i] = total_cycles;
                }

                double mean_cycles = calculate_mean(results, reps);
                double median_cycles = calculate_median(results, reps);
                double stddev_cycles = calculate_standard_deviation(results, reps, mean_cycles);

                printf("%d, %lu, %lf, %lf, %lf\n", numThreads, bytes, mean_cycles, median_cycles, stddev_cycles);

            }
        }
    }

    return 0;
}