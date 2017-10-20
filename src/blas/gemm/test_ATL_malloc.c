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
  a = simple_malloc(&blob, 12);
  simple_free(&blob, a);

  /* Out of memory (while nothing is allocated) */
  a = simple_malloc(&blob, MEM_SIZE + 12);
  assert(!a && "Allocated more mem than available.");

  /* Allocate zero bytes. */
  a = simple_malloc(&blob, 0);
  simple_free(&blob, a);

  /* Allocate multiple times, free to test clean up. */
  a = simple_malloc(&blob, 7);
  b = simple_malloc(&blob, 12);
  simple_free(&blob, a);
  simple_free(&blob, b);

  /* Same as above but little more complex */
  a = simple_malloc(&blob, 7);
  b = simple_malloc(&blob, 0);
  c = simple_malloc(&blob, 12);
  simple_free(&blob, b);
  simple_free(&blob, c);
  simple_free(&blob, a);

  /* Exhaust mem allocation. */
  int* pointers[MEM_SIZE];
  for (unsigned i = 0; i < MEM_SIZE; i++) {
    pointers[i] = 0;
  }
  unsigned primes[] = { 23, 19, 17, 13, 11, 7, 5, 3, 2 , 1, 0 };
  unsigned p = 0;
  for (unsigned i = 0; i < MEM_SIZE; i++) {
    pointers[i] = simple_malloc(&blob, primes[p]);
    if (!pointers[i])
      p++;
    if (!primes[p])
      break;
  }

  for (unsigned i = 0; i < MEM_SIZE; i++) {
    if (pointers[i])
      simple_free(&blob, pointers[i]);
  }

  return 0;
}
