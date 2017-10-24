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

#ifndef ATLAS_MALLOC_H
#define ATLAS_MALLOC_H

/*
 * Copied from atlas_lvl3.h
 */
#ifndef ATL_MaxMalloc
   #ifdef ATL_MaxMalloc_MB
      #define ATL_MaxMalloc (((size_t)(ATL_MaxMalloc_MB))<<20)
   #else
      #define ATL_MaxMalloc 67108864
   #endif
#endif

typedef unsigned mem_size_t;

typedef struct {
  /*
   * Block from the memory is allocated.
   */
  void* blob;
  /*
   * Size of the memory block.
   */
  mem_size_t size;
  /* Pointer to last allocated block. Only NULL if no memory has been
   * allocated
   */
  void* tail;
} MemBlob;

/*
 * nclude here because of MemBlob declaration.
 */
#include "atlas_misc.h"

#ifndef DIRECTHSA
extern MemBlob* globalMemBlob;
#endif

#if defined(DIRECTHSA)
#undef PHSA
#define PHSA _hsa_function
#endif

HSA_FUNCTION
void Mjoin(init_mem_blob,PHSA)(MemBlob* memblob,
                               void* blob, mem_size_t size);
HSA_FUNCTION
void* Mjoin(simple_malloc,PHSA)(MemBlob* memblob, mem_size_t size);
HSA_FUNCTION
void Mjoin(simple_free,PHSA)(MemBlob* memblob, void* prt);

#endif
