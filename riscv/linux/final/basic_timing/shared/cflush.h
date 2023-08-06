// Based on code from https://github.com/ucb-bar/hwacha-template/blob/master/tests/rocc.h

#ifndef C_FLUSH_H
#define C_FLUSH_H

#include <stdint.h>
#include <assert.h>

#define STR1(x) #x
#ifndef STR
#define STR(x) STR1(x)
#endif

#define CFLUSH_D_L1   0x10200f //0xCC000073 //Original opcode is 0xFC000073
#define CDISCARD_D_L1 0xCC200073
#define CBO_CLEAN	  0x10200f
#define CBO_FLUSH	  0x20200f

#define CUSTOMX(X, xd, xs1, xs2, rd, rs1, rs2) 		  \
  (X)                                   |             \
  (rd                 << (7))           |             \
  (xs2                << (7+5))         |             \
  (xs1                << (7+5+1))       |             \
  (xd                 << (7+5+2))       |             \
  (rs1                << (7+5+3))       |             \
  (rs2                << (7+5+3+5))


#define CFLUSH_D_L1_ALL() CUSTOMX(CFLUSH_D_L1,0,0,0,0,0,0)
#define CFLUSH_D_L1_REG(rs1) CUSTOMX(CFLUSH_D_L1,0,0,0,0,rs1,0)
#define CDISCARD_D_L1_REG(rs1) CUSTOMX(CDISCARD_D_L1,0,0,0,0,rs1,0)
#define CBO_CLEAN_REG(rs1) CUSTOMX(CBO_CLEAN,0,0,0,0,rs1,0)
#define CBO_FLUSH_REG(rs1) CUSTOMX(CBO_FLUSH,0,0,0,0,rs1,0)

#define FLUSH_D_ALL() { 								   \
    asm volatile(".word " STR(CFLUSH_D_L1_ALL()) "\n\t" :: \
				 : "memory");                              \
}


//rs1 is address
//rs1_n is the register number to use
#define CFLUSH_D_L1_INST(rs1, rs1_n) {                                          \
	register uint64_t rs1_ asm("x" #rs1_n) = (uint64_t) rs1;                    \
	asm volatile(".word " STR(CFLUSH_D_L1_REG(rs1_n)) "\n\t" ::[_rs1] "r"(rs1_) \
			     : "memory");                                                   \
}

#define FLUSH_D_REG(rs1) CFLUSH_D_L1_INST(rs1, 13)

//rs1 is address
//rs1_n is the register number to use
#define CBO_FLUSH_INST(rs1, rs1_n) {                                          \
	register uint64_t rs1_ asm("x" #rs1_n) = (uint64_t) rs1;                    \
	asm volatile(".word " STR(CBO_FLUSH_REG(rs1_n)) "\n\t" ::[_rs1] "r"(rs1_) \
			     : "memory");                                                   \
}

#define CBO_FLUSH_FN(rs1) CBO_FLUSH_INST(rs1, 13)

//rs1 is address
//rs1_n is the register number to use
#define CBO_CLEAN_INST(rs1, rs1_n) {                                          \
	register uint64_t rs1_ asm("x" #rs1_n) = (uint64_t) rs1;                    \
	asm volatile(".word " STR(CBO_CLEAN_REG(rs1_n)) "\n\t" ::[_rs1] "r"(rs1_) \
			     : "memory");                                                   \
}

#define CBO_CLEAN_FN(rs1) CBO_CLEAN_INST(rs1, 13)

#define read_csr(reg) ({ unsigned long __tmp; \
  asm volatile ("csrr %0, " #reg : "=r"(__tmp)); \
  __tmp; })
#define rdtime() read_csr(time)
#define rdcycle() read_csr(cycle)

#endif  // C_FLUSH_H
