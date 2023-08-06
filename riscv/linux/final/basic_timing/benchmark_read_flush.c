#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <cflush.h>
#include <math.h>

int bytes_per_thread;
volatile uint64_t total_cycles[8];

void* threadFunc(void* arg) {
    int bytes = bytes_per_thread;

    int threadid = *(int *)arg;

    void *x = malloc(bytes);

    // dirty each line
    for (int i = 0; i < bytes; i+=64) {
        *((uint64_t *) (x+i)) = i; 
    }

    // //flush
    uint64_t start = read_csr(cycle);
    
    for (int i = 0; i < bytes; i += 64) {
        CBO_FLUSH_FN(x+i);
    }

    asm volatile ("fence rw, rw");
    
    uint64_t end = read_csr(cycle);
    uint64_t elapsed = end - start;
    total_cycles[threadid] += elapsed;

    free(x);

    return NULL;
}

void bench(int numThreads, int bytes) {
// Initialize thread identifiers
    pthread_t* threads = malloc(numThreads * sizeof(pthread_t));

    bytes_per_thread = bytes/numThreads;

    int threadid[numThreads];
    // Create threads
    for (int i = 0; i < numThreads; i++) {
        threadid[i] = i;
        if (pthread_create(&threads[i], NULL, threadFunc, &threadid[i]) != 0) {
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
    if (argc != 3) {
        printf("Usage: ./benchmark <num_threads> <reps for confidence>\n");
        return 1;
    }

    int numThreads = atoi(argv[1]);
    int reps = atoi(argv[2]);
    //int sameX = atoi(argv[3]);

    printf("threads, bytes, mean, median, stddev, max, min, Q1, Q3\n");

    int starter_bytes = numThreads*64;

    for(int bytes=starter_bytes;bytes<=16384;bytes*=2) {

        uint64_t results[reps];

        for(int i=0;i<reps;i++) {
			for (int k = 0; k < numThreads; ++k) {
				total_cycles[k] = 0;
			}
            bench(numThreads, bytes);
            /* find max of each thread */
			results[i] = total_cycles[0];
			for (int j = 1; j < numThreads; ++j) {
				if (total_cycles[j] > total_cycles[j - 1]) {
					results[i] = total_cycles[j];
				}
			}
        }

    	qsort(results, reps, sizeof(uint64_t), compare_uint64_t);
        double mean_cycles = calculate_mean(results, reps);
        double median_cycles = calculate_median(results, reps);
        double stddev_cycles = calculate_standard_deviation(results, reps, mean_cycles);
		double Q1_cycles = calculate_median(results, reps/2);
		double Q3_cycles;		
		if (reps % 2  == 0) {
			Q3_cycles = calculate_median(&results[reps/2], reps/2);
		} else {
			Q3_cycles = calculate_median(&results[reps/2 + 1], reps/2);
		}
		uint64_t min_cycles = results[0];
		uint64_t max_cycles = results[reps - 1];

        printf("%d, %u, %lf, %lf, %lf, %lu, %lu, %lf, %lf\n", numThreads, bytes, mean_cycles, median_cycles, stddev_cycles, max_cycles, min_cycles, Q1_cycles, Q3_cycles);

    }

    return 0;
}
