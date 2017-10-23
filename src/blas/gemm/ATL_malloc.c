/*
 *             Automatically Tuned Linear Algebra Software v3.10.3
 *                    (C) Copyright 2017 Parmance
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "atlas_malloc.h"

#  include <stdint.h> /* uintptr_t*/

#define POINTER_SIZE sizeof(void*)
#define MAX(x, y) (x) > (y) ? (x) : (y)

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

typedef struct mem_block_header {
  /* Block is free if non-zero. */
  unsigned free;
  /* Size of the reserved block in bytes. Always aligned to  */
  mem_size_t size;
  /* Pointer to previous reserved block. */
  struct mem_block_header* prev;
} mem_block_header;



void Mjoin(init_mem_blob,PHSA)(MemBlob* memBlob, void* blob,
                               mem_size_t size)
{
  memBlob->blob = blob;
  memBlob->size = size;
  memBlob->tail = NULL;
}

HSA_FUNCTION
static void* Mjoin(align_for_header,PHSA)(void* ptr) {
  /* Can not use c11's alignof() :(  */
  const uintptr_t align = MAX(POINTER_SIZE, sizeof(mem_size_t));
  if (((uintptr_t)ptr % align) != 0)
    ptr = (char*)ptr + (align - ((uintptr_t)ptr % align));
  return ptr;
}

/* Return pointer to next location for new block allocation.
   Pointer is aligned for the block header. */

HSA_FUNCTION
static void* Mjoin(next_block_loc,PHSA)(mem_block_header* blkHeader) {
  void* next =
     (char*)blkHeader + sizeof(mem_block_header) + blkHeader->size;
  return Mjoin(align_for_header,PHSA)(next);
}

/* Book block reservation and return pointer to block data - the
   memory the user asked for.  */

HSA_FUNCTION
void* Mjoin(reserve_block,PHSA)(mem_block_header* newTail,
                               mem_block_header* prev,
                               mem_size_t size) {
  newTail->size = size;
  newTail->free = 0;
  newTail->prev = prev;
  return (char*)newTail + sizeof(mem_block_header);
}

/* Return bytes of available memory for the allocation. */

HSA_FUNCTION
static mem_size_t Mjoin(available_mem,PHSA)(MemBlob* memBlob) {
  if (memBlob->tail == NULL)
  {
    if (sizeof(mem_block_header) > memBlob->size)
      return 0;
    return memBlob->size - sizeof(mem_block_header);
  }
  mem_block_header* tail = (mem_block_header*)memBlob->tail;
  uintptr_t next = (uintptr_t)Mjoin(next_block_loc,PHSA)(tail) +
     sizeof(mem_block_header);
  uintptr_t end = (uintptr_t)((char*)memBlob->blob + memBlob->size);
  if (next > end)
    return 0;
  return end - next;
}

/* Moves tail pointer for freed memory. Return pointer to next topmost
   reserved block or NULL if all memory is freed. */

HSA_FUNCTION
static mem_block_header* Mjoin(clean_up,PHSA) (mem_block_header* tail)
{
  while (tail != NULL && tail->free) tail = tail->prev;
  return tail;
}

HSA_FUNCTION
void* Mjoin(simple_malloc,PHSA)(MemBlob* memBlob, mem_size_t size)
{
   if (memBlob == NULL) return NULL;

  mem_size_t freeMem = Mjoin(available_mem,PHSA)(memBlob);
  if (size > freeMem)
  {
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

  mem_block_header* tail = (mem_block_header*)memBlob->tail;
  mem_block_header* next = NULL;
  mem_block_header* prev = NULL;
  if (tail == NULL)
     memBlob->tail = tail = next =
        (mem_block_header*)Mjoin(align_for_header,PHSA)(memBlob->blob);
  /* DirectHSA specific kludge fix. Without this GCC emits
   * __builtin_trap() which fails during BRIG finalization.
   */
  if (tail == NULL) return NULL;
  if (next == NULL)
  {
     next = (mem_block_header*)Mjoin(next_block_loc,PHSA)(tail);
     prev = tail;
     memBlob->tail = next;
  }

  void* addr = Mjoin(reserve_block,PHSA)(next, prev, size);
  return addr;
}

HSA_FUNCTION
void Mjoin(simple_free,PHSA)(MemBlob* memBlob, void* ptr)
{
  if (ptr == NULL) /* Quick return for NULL. */
    return;

  mem_block_header* tail = (mem_block_header*)memBlob->tail;
  mem_block_header* found = tail;
  /* TODO may just assume the ptr is always valid. For free() it is
     undefined behavior to double free, free unallocated mem
     etc. anyway */
  /* Pointer to the assumed block header. */
  void* toSearch = (char*)ptr - sizeof(mem_block_header);

  while (found != NULL && found != toSearch) found = found->prev;

  if (found) found->free = 1;

  /* TODO update prev pointer of the block after the freed one? */
  memBlob->tail = Mjoin(clean_up,PHSA)(tail);
}
