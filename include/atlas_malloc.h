
#ifndef ATLAS_MALLOC_H
#define ATLAS_MALLOC_H


/* TODO: fix for code style. */

/* Copied from atlas_lvl3.h */
#ifndef ATL_MaxMalloc
   #ifdef ATL_MaxMalloc_MB
      #define ATL_MaxMalloc (((size_t)(ATL_MaxMalloc_MB))<<20)
   #else
      #define ATL_MaxMalloc 67108864
   #endif
#endif

typedef unsigned mem_size_t;

/* TODO make this opaque. */
typedef struct {
  /* Block from the memory is allocated. */
  void* blob;
  /* Size of the memory block. */
  mem_size_t size;
  /* Pointer to last allocated block. Only NULL if no memory has been
     allocated */
  void* tail;
} MemBlob;

/* Include here because of MemBlob declaration. */
#include "atlas_misc.h"

#ifndef DIRECTHSA
extern MemBlob* globalMemBlob;
#endif

/* Kludge. PHSA and PHSA_FN macros are not what is wanted because no
 * -DSREAL is not given.  */
#if defined(DIRECTHSA)
#  if defined(HSADECLS) /* DEVTEMP */
#    undef PHSA_FN
#    define PHSA_FN _hsa_function
#  else
#    undef PHSA_FN
#    define PHSA_FN _hsa
#  endif
#endif

HSA_FUNCTION
void Mjoin(init_mem_blob,PHSA_FN)(MemBlob* memblob,
                                  void* blob, mem_size_t size);
HSA_FUNCTION
void* Mjoin(simple_malloc,PHSA_FN)(MemBlob* memblob, mem_size_t size);
HSA_FUNCTION
void Mjoin(simple_free,PHSA_FN)(MemBlob* memblob, void* prt);

#endif
