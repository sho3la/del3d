// del3d - slot allocator acceptance test. Self-contained.
//
// The container is not interesting in itself; it exists to fix the order in
// which cells and vertices are handed out, which is part of del3d's output
// contract. This pins the four properties that order depends on:
//   * block sizes: 14 first, then 16 more each time (14, 30, 46, ...)
//   * a fresh block hands out its slots front to back
//   * destroy/create is LIFO, so a freed slot is the next one reused
//   * iteration is memory order, skipping free slots
#include "tds.h"
#include "compact_container.h"
#include <cstdio>
using namespace del3d::detail;
int main(){
  CompactContainer<int> cc;
  // A fresh block is consumed in ascending index order.
  for (int i=0;i<14;i++) { int k=cc.create(); if(k!=i){printf("FAIL alloc order %d!=%d\n",k,i);return 1;} }
  int k14 = cc.create();               // exhausts block 1, forces block 2
  if (k14 != 14) { printf("FAIL new block start %d\n", k14); return 1; }
  if (cc.capacity() != 14+30) { printf("FAIL capacity %zu\n", cc.capacity()); return 1; }
  // Freed slots come back in reverse order of freeing.
  cc.destroy(5); cc.destroy(9);
  int a=cc.create(), b=cc.create();
  if (a!=9 || b!=5) { printf("FAIL LIFO %d %d\n",a,b); return 1; }
  // Iteration is ascending index order and skips the hole left at 3.
  cc.destroy(3);
  int prev=-1, count=0; bool ok=true;
  cc.for_each([&](CompactContainer<int>::Index i){ if(i<=prev) ok=false; prev=i; ++count; });
  if(!ok || count!=14){ printf("FAIL iteration ok=%d count=%d\n",(int)ok,count); return 1; }
  Tds tds; (void)tds;   // the data structure must at least be constructible
  printf("compact_container + tds: OK (blocks 14/30, LIFO reuse, memory-order iteration)\n");
  return 0;
}
