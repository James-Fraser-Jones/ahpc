#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MASTER  0
#define NSPEEDS 9
#define NY 28
#define NX 3
#define heightH 4

typedef struct
{
  float speeds[NSPEEDS];
} t_speed;

int main(){

  printf("Grid size = %lu = %lu * %d * %d * %d\n", (sizeof(t_speed) * (heightH * NX)),
  sizeof(float), NSPEEDS, heightH, NX);

  t_speed* cells = (t_speed*)malloc(sizeof(t_speed) * (heightH * NX));

  return 0;
}

/*

[][][][][] ... grid as 1-d

[][][][] ... a row of the grid

[t_speed] ... a cell of the grid

[fffffffff] ... a t_speed of a cell

[bbbb,bbbb,bbbb,bbbb...] a set of floats of a t_speed

All memory is completely contiguous so there is no gaps between any of the bytes whatsoever.

//the size of t_speed is only just large enough to store NSPEEDS number of floats
//conclusion: every byte of memory used in the t_speed struct is used to represent
//one of the 9 floats that it's storing

//malloc allocates contiguous memory to store all the structs so it is literally a
//series of 4 * NSPEEDS * heightH * NX bytes with every byte representing one of the
//NSPEEDS * heightH * NX floats, every float representing one of the heightH * NX t_speed
//structs and every t_speed struct representing one of the NX cells in the grid
*/
