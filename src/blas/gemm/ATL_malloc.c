/* TODO add copyright. */

#include "atlas_malloc.h"

#  include <stdint.h> /* uintptr_t*/

#define POINTER_SIZE sizeof(void*)
/* Kludge fix. There is bug in the HLC. sizeof(void*) returns 4 for
 * hsail64 target where pointer is 64 bits. */
/* #define POINTER_SIZE 8 */
#define MAX(x, y) (x) > (y) ? (x) : (y)

/* #define DEBUG_SIMPLE_MALLOC */
#ifdef DEBUG_SIMPLE_MALLOC
#include <stdio.h>
#include <assert.h>
#  define ASSERT(x) assert(x)
#  define DEBUG_MSG printf
#  define PING printf("%d: ping\n", __LINE__)
#else
#  define ASSERT(x)
#  define DEBUG_MSG (void)
#  define PING do {} while(0);
#endif


/* TODO fix for code style. */

/* TODO describtion for simple malloc
   - not thread safe.
   - Fragmentation issue.
   - User may overrun the book keeping.
   - Freeing non-allocated mem has undefined behavior.
   - stack behavior: mem is allocated after the last allocated block in
     address.
   - no (type)malloc(n) kind of alignment.
*/

/* TODO handle cases
   - blob is smaller than header.
   - allocating zero bytes. */

#ifndef DIRECTHSA
static char globalMem[ATL_MaxMalloc];
static MemBlob globalMemBlobData = { &globalMem, sizeof(globalMem), NULL };
MemBlob* globalMemBlob = &globalMemBlobData;
#endif

typedef struct MemBlockHeader MemBlockHeader;

struct MemBlockHeader {
  /* Block is free if non-zero. */
  unsigned free;
  /* Size of the reserved block in bytes. Always aligned to  */
  mem_size_t size;
  /* Pointer to previous reserved block. */
  MemBlockHeader* prev;
};



void Mjoin(init_mem_blob,PHSA_FN)(MemBlob* memBlob, void* blob,
                               mem_size_t size)
{
  memBlob->blob = blob;
  memBlob->size = size;
  memBlob->tail = NULL;
  /* DEBUG_MSG ("init blob. Start addr = %p\n", blob); */
}

/* Get block header for the pointer returned by simple_malloc.
   Return NULL if error. */
/* static MemBlockHeader* getHeader(MemBlob* memBlob, void* ptr) */
/* { */
/*   /\* Note: risky assumption that the ptr is valid. *\/ */
/*   return (MemBlockHeader*)ptr - sizeof(MemHeaderBlock); */
/* } */

/* Aligns pointer for block header. */

HSA_FUNCTION
static void* Mjoin(alignForHeader,PHSA_FN)(void* ptr) {
  /* Can not use c11's alignof() :(  */
  const uintptr_t align = MAX(POINTER_SIZE, sizeof(mem_size_t));
  if (((uintptr_t)ptr % align) != 0)
    ptr = (char*)ptr + (align - ((uintptr_t)ptr % align));
  return ptr;
}

/* Return pointer to next location for new block allocation.
   Pointer is aligned for the block header. */

HSA_FUNCTION
static void* Mjoin(nextBlockLoc,PHSA_FN)(MemBlockHeader* blkHeader) {
  ASSERT(blkHeader);
  void* next =
     (char*)blkHeader + sizeof(MemBlockHeader) + blkHeader->size;
  return Mjoin(alignForHeader,PHSA_FN)(next);
}

/* Book block reservation and return pointer to block data - the
   memory the user asked for.  */

HSA_FUNCTION
void* Mjoin(reserveBlock,PHSA_FN)(MemBlockHeader* newTail,
                               MemBlockHeader* prev,
                               mem_size_t size) {
  newTail->size = size;
  newTail->free = 0;
  newTail->prev = prev;
  /* DEBUG_MSG("init header at %p. (prev=%p)\n", */
  /*           (void*)newTail, (void*)prev); */
  return (char*)newTail + sizeof(MemBlockHeader);
}

/* Return bytes of available memory for the allocation. */

HSA_FUNCTION
static mem_size_t Mjoin(available_mem,PHSA_FN)(MemBlob* memBlob) {
  if (memBlob->tail == NULL) {
    if (sizeof(MemBlockHeader) > memBlob->size)
      return 0;
    return memBlob->size - sizeof(MemBlockHeader);
  }
  MemBlockHeader* tail = (MemBlockHeader*)memBlob->tail;
  uintptr_t next = (uintptr_t)Mjoin(nextBlockLoc,PHSA_FN)(tail) +
     sizeof(MemBlockHeader);
  uintptr_t end = (uintptr_t)((char*)memBlob->blob + memBlob->size);
  if (next > end)
    return 0;
  return end - next;
}

/* Moves tail pointer for freed memory. Return pointer to next topmost
   reserved block or NULL if all memory is freed. */

HSA_FUNCTION
static MemBlockHeader* Mjoin(cleanUp,PHSA_FN) (MemBlockHeader* tail)
{
  while (tail != NULL && tail->free) {
    /* DEBUG_MSG("cleaning up block (%p)\n", (void*)tail); */
    tail = tail->prev;
  }
  return tail;
}

HSA_FUNCTION
void* Mjoin(simple_malloc,PHSA_FN)(MemBlob* memBlob, mem_size_t size)
{
   if (memBlob == NULL) return NULL;

  mem_size_t freeMem = Mjoin(available_mem,PHSA_FN)(memBlob);
  if (size > freeMem) {
     /* DEBUG_MSG("out of memory (tried alloc %llub, %llub left).\n", */
     /*           size, freeMem); */
#ifndef DIRECTHSA
     static unsigned fail_count = 0;
     if ((fail_count % 1000) == 0)
        printf("simple malloc failed: "
               "out of memory (tried alloc %llub, %llub left).\n",
               size, freeMem);
     fail_count++;
#endif
    return NULL;
  }

  MemBlockHeader* tail = (MemBlockHeader*)memBlob->tail;
  MemBlockHeader* next = NULL;
  MemBlockHeader* prev = NULL;
  if (tail == NULL)
     memBlob->tail = tail = next =
        (MemBlockHeader*)Mjoin(alignForHeader,PHSA_FN)(memBlob->blob);
  /* DirectHSA specific check. Without this GCC emits __builtin_trap()
   * which fails during finalization.  */
  if (tail == NULL) return NULL;
  if (next == NULL) {
     next = (MemBlockHeader*)Mjoin(nextBlockLoc,PHSA_FN)(tail);
     prev = tail;
     memBlob->tail = next;
  }

  void* addr = Mjoin(reserveBlock,PHSA_FN)(next, prev, size);
  /* DEBUG_MSG("Allocated mem at %p of size %llub (%llub left).\n", */
  /*           addr, size, Mjoin(available_mem,PHSA_FN)(memBlob)); */
  return addr;
}


HSA_FUNCTION
void Mjoin(simple_free,PHSA_FN)(MemBlob* memBlob, void* ptr)
{
  if (ptr == NULL) /* Quick return for NULL. */
    return;

  MemBlockHeader* tail = (MemBlockHeader*)memBlob->tail;
  MemBlockHeader* found = tail;
  /* TODO may just assume the ptr is always valid. For free() it is
     undefined behavior to double free, free unallocated mem
     etc. anyway */
  /* Pointer to the assumed block header. */
  void* toSearch = (char*)ptr - sizeof(MemBlockHeader);

  while (found != NULL && found != toSearch) {
    found = found->prev;
  }

  if (found) {
    found->free = 1;
    /* DEBUG_MSG("Free mem at %p\n", ptr); */
  } else {
    /* DEBUG_MSG("Nothing to free at %p\n", ptr); */
  }

  /* TODO update prev pointer of the block after the freed one? */
  memBlob->tail = Mjoin(cleanUp,PHSA_FN)(tail);
}
