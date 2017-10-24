#include <assert.h>

#include "atlas_malloc.h"

#define MEM_SIZE 128

int
main()
{
  char mem[MEM_SIZE];
  MemBlob blob;
  init_mem_blob(&blob, mem, sizeof(mem));

  char *a, *b, *c;

  /* Smoke test */
  a = ATL_Malloc(&blob, 12);
  ATL_Free(&blob, a);

  /* Out of memory (while nothing is allocated) */
  a = ATL_Malloc(&blob, MEM_SIZE + 12);
  assert(!a && "Allocated more mem than available.");

  /* Allocate zero bytes. */
  a = ATL_Malloc(&blob, 0);
  ATL_Free(&blob, a);

  /* Allocate multiple times, free to test clean up. */
  a = ATL_Malloc(&blob, 7);
  b = ATL_Malloc(&blob, 12);
  ATL_Free(&blob, a);
  ATL_Free(&blob, b);

  /* Same as above but little more complex */
  a = ATL_Malloc(&blob, 7);
  b = ATL_Malloc(&blob, 0);
  c = ATL_Malloc(&blob, 12);
  ATL_Free(&blob, b);
  ATL_Free(&blob, c);
  ATL_Free(&blob, a);

  /* Exhaust mem allocation. */
  int* pointers[MEM_SIZE];
  for (unsigned i = 0; i < MEM_SIZE; i++) {
    pointers[i] = 0;
  }
  unsigned primes[] = { 23, 19, 17, 13, 11, 7, 5, 3, 2 , 1, 0 };
  unsigned p = 0;
  for (unsigned i = 0; i < MEM_SIZE; i++) {
    pointers[i] = ATL_Malloc(&blob, primes[p]);
    if (!pointers[i])
      p++;
    if (!primes[p])
      break;
  }

  for (unsigned i = 0; i < MEM_SIZE; i++) {
    if (pointers[i])
      ATL_Free(&blob, pointers[i]);
  }

  return 0;
}
