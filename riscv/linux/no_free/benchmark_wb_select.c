#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <cflush.h>
#include <math.h>

int bytes_per_thread;

uint64_t results_raw[10][8];

void* threadFunc(void* arg) {
    int thread = *(int *)arg;
    int bytes = bytes_per_thread;

    for (int r = 0; r < 10; r++) {
    // dirty each line

        void *x = malloc(bytes);
        for (int i = 0; i < bytes; i+=64) {
            *((uint64_t *) (x+i)) = 1; 
        }

        // //flush
        uint64_t start = read_csr(cycle);
        
        for (int i = 0; i < bytes; i += 64) {
            CBO_CLEAN_FN(x+i);
        }
        asm volatile ("fence rw, rw");
        
        uint64_t end = read_csr(cycle);
        uint64_t elapsed = end - start;
        results_raw[r][thread] = elapsed;
    }

    return NULL;
}

void bench(int numThreads, int bytes) {
// Initialize thread identifiers
    pthread_t* threads = malloc(numThreads * sizeof(pthread_t));

    bytes_per_thread = bytes/numThreads;

    int thread_ids[numThreads];

    // Create threads
    for (int i = 0; i < numThreads; i++) {
        thread_ids[i] = i;
        if (pthread_create(&threads[i], NULL, threadFunc, &(thread_ids[i])) != 0) {
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
    if (argc != 4) {
        printf("Usage: ./benchmark <num_threads> <reps for confidence> <bytes>\n");
        return 1;
    }

    int numThreads = atoi(argv[1]);
    int reps = atoi(argv[2]);
    int bytes = atoi(argv[3]);

    printf("threads, bytes, mean, median, stddev\n");

    uint64_t results[10];

    for(int i=0;i<10;i++)
        results[i] = 0;

    bench(numThreads, bytes);

    for(int i=0;i<10;i++)
        for(int j=0;j<numThreads;j++)
            results[i] += results_raw[i][j];

    double mean_cycles = calculate_mean(results, 10);
    double median_cycles = calculate_median(results, 10);
    double stddev_cycles = calculate_standard_deviation(results, 10, mean_cycles);

    printf("%d, %d, %lf, %lf, %lf\n", numThreads, bytes, mean_cycles, median_cycles, stddev_cycles);

    return 0;
}
