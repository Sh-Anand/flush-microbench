// See LICENSE for license details.

#include "common.h"
#include "cflush.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "util.h"

#define U32 *(volatile unsigned int *)
#define MEMORY_MEM_ADDR 0x80000000UL

void thread_entry(int cid, int nc) {

    int buf[2048]={};
    if(cid==0) { //core 0
        for(unsigned int i=0;i<256;++i) {
            U32(MEMORY_MEM_ADDR+0x80000+i*4) = i*4;
        }
    } else if(cid==1){ //core 1
        for(volatile int cnt=0;cnt<=20;++cnt);
        for(unsigned int i=0;i<256;++i) {
            buf[i]=U32(MEMORY_MEM_ADDR+0x80000+i*4);
        }
    }

    FLUSH_D_REG(0x80080000UL);
    FLUSH_D_REG(0x80080040UL);
    FLUSH_D_REG(0x80080080UL);
  
    barrier(nc);
    exit(0);
}
