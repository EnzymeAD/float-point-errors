#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "fcn.h"
#include "mesh.h"
#include "opt.h"

#ifndef MICROSEC
#  define MICROSEC 1000000
#endif

int main(int argc, char **argv)
{
  Mesh *mesh = NULL;

  double m, m1, m2;
  double s, s1, s2;
  double t1, t2;

  struct rusage r0, r1, r2;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <input_mesh> <output_mesh>\n", argv[0]);
    return -1;
  }

  getrusage(RUSAGE_SELF, &r0);
  if (readMesh(argv[1], &mesh)) {
    freeMesh(&mesh);
    return -1;
  }
  getrusage(RUSAGE_SELF, &r1);

  m1 = (double) r1.ru_utime.tv_usec;
  m2 = (double) r0.ru_utime.tv_usec;
  m = m1 - m2;
    
  s1 = (double) r1.ru_utime.tv_sec;
  s2 = (double) r0.ru_utime.tv_sec;
  s = s1 - s2;

  t1 = s + m / MICROSEC;

  m1 = (double) r1.ru_stime.tv_usec;
  m2 = (double) r0.ru_stime.tv_usec;
  m = m1 - m2;

  s1 = (double) r1.ru_stime.tv_sec;
  s2 = (double) r0.ru_stime.tv_sec;
  s = s1 - s2;

  t2 = s + m / MICROSEC;

  printf("Read: %5.4e System: %5.4e Total: %5.4e\n", t1, t2, t1+t2);

  getrusage(RUSAGE_SELF, &r1);
  optMesh(mesh, 1.0e-6, 2); 
  getrusage(RUSAGE_SELF, &r2);

  m1 = (double) r2.ru_utime.tv_usec;
  m2 = (double) r1.ru_utime.tv_usec;
  m = m1 - m2;
    
  s1 = (double) r2.ru_utime.tv_sec;
  s2 = (double) r1.ru_utime.tv_sec;
  s = s1 - s2;

  t1 = s + m / MICROSEC;

  m1 = (double) r2.ru_stime.tv_usec;
  m2 = (double) r1.ru_stime.tv_usec;
  m = m1 - m2;

  s1 = (double) r2.ru_stime.tv_sec;
  s2 = (double) r1.ru_stime.tv_sec;
  s = s1 - s2;

  t2 = s + m / MICROSEC;

  printf("Optimize: User: %5.4e System: %5.4e Total: %5.4e\n", t1, t2, t1+t2);

  getrusage(RUSAGE_SELF, &r1);
  writeMesh(argv[2], mesh);
  getrusage(RUSAGE_SELF, &r2);

  m1 = (double) r2.ru_utime.tv_usec;
  m2 = (double) r1.ru_utime.tv_usec;
  m = m1 - m2;
    
  s1 = (double) r2.ru_utime.tv_sec;
  s2 = (double) r1.ru_utime.tv_sec;
  s = s1 - s2;

  t1 = s + m / MICROSEC;

  m1 = (double) r2.ru_stime.tv_usec;
  m2 = (double) r1.ru_stime.tv_usec;
  m = m1 - m2;

  s1 = (double) r2.ru_stime.tv_sec;
  s2 = (double) r1.ru_stime.tv_sec;
  s = s1 - s2;

  t2 = s + m / MICROSEC;

  printf("Write: User: %5.4e System: %5.4e Total: %5.4e\n", t1, t2, t1+t2);

  freeMesh(&mesh);
  
  getrusage(RUSAGE_SELF, &r2);
  m1 = (double) r2.ru_utime.tv_usec;
  m2 = (double) r0.ru_utime.tv_usec;
  m = m1 - m2;
    
  s1 = (double) r2.ru_utime.tv_sec;
  s2 = (double) r0.ru_utime.tv_sec;
  s = s1 - s2;

  t1 = s + m / MICROSEC;

  m1 = (double) r2.ru_stime.tv_usec;
  m2 = (double) r0.ru_stime.tv_usec;
  m = m1 - m2;

  s1 = (double) r2.ru_stime.tv_sec;
  s2 = (double) r0.ru_stime.tv_sec;
  s = s1 - s2;

  t2 = s + m / MICROSEC;

  printf("Total: User: %5.4e System: %5.4e Total: %5.4e\n", t1, t2, t1+t2);
  return 0;
}

