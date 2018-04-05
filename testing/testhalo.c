#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define MASTER  0
#define NSPEEDS 9
#define NY 28
#define NX 3
#define MAXITERS 20

//local compile and execute commands:
//mpicc -std=c99 -Wall -O3 testhalo.c -lm -o testhalo
//mpiexec -n 4 ./testhalo

typedef struct
{
  int    nx;            /* no. of cells in x-direction */
  int    ny;            /* no. of cells in y-direction */
  int    maxIters;	/* no. of iterations */
  int    reynolds_dim;  /* dimension for Reynolds number */
  float density;       /* density per link */
  float accel;         /* density redistribution */
  float omega;         /* relaxation parameter */
  int heightH; //ny for local process (including halos)
} t_param;

typedef struct
{
  float speeds[NSPEEDS];
} t_speed;

void init(t_param* params, t_speed** cells_ptr, int myRank, int numProcs, float** av_vels_ptr){
  params->ny = NY;
  params->nx = NX;
  params->maxIters = MAXITERS;
  params->heightH = NY/numProcs + 2;
  *cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params->heightH * params->nx)); /* grid containing fluid densities */

  for (int jj = 0; jj < params->heightH; jj++){
    for (int ii = 0; ii < params->nx; ii++){
      for (int kk = 0; kk < NSPEEDS; kk++){
        (*cells_ptr)[ii + jj*params->nx].speeds[kk] = 0.0f;
      }
    }
  }

  for (int jj = 1; jj < params->heightH-1; jj++){
    for (int ii = 0; ii < params->nx; ii++){
      for (int kk = 0; kk < NSPEEDS; kk++){
        (*cells_ptr)[ii + jj*params->nx].speeds[kk] = myRank*(params->heightH-2) + jj + kk + ii;
      }
    }
  }

  *av_vels_ptr = (float*)malloc(sizeof(float) * params->maxIters);
}

void print(t_param* params, t_speed* cells, int myRank){
  printf("Rank: %d\n", myRank);
  for (int jj = params->heightH-1; jj >= 0; jj--){
    for (int ii = 0; ii < params->nx; ii++){
      printf("(");
      for (int kk = 0; kk < NSPEEDS; kk++){
        float speed = cells[ii + jj*params->nx].speeds[kk];
        printf("%.2g, ", speed);
      }
      printf(") ");
    }
    printf("\n");
  }
  printf("\n\n");
}

int halo_exchange(const t_param params, t_speed* cells, int up, int down, MPI_Status status,
                  t_speed* botHalo, t_speed* botRow, t_speed* topRow, t_speed* topHalo){

  MPI_Sendrecv(topRow, NSPEEDS*params.nx, MPI_FLOAT, up, 0,
        botHalo, NSPEEDS*params.nx, MPI_FLOAT, down, 0, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(botRow, NSPEEDS*params.nx, MPI_FLOAT, down, 0,
        topHalo, NSPEEDS*params.nx, MPI_FLOAT, up, 0, MPI_COMM_WORLD, &status);

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[]){

  t_param params;
  t_speed* cells;

  int myRank;
  int numProcs;
  int up;
  int down;
  MPI_Status status;     /* struct used by MPI_Recv */
  t_speed* botHalo;
  t_speed* topHalo;
  t_speed* botRow;
  t_speed* topRow;
  float* av_vels = NULL;

  char sendbuf[BUFSIZ]; //these are just used for the printing messages
  char recvbuf[BUFSIZ];

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

  up = (myRank + 1) % numProcs;
  down = (myRank == 0) ? (myRank + numProcs - 1) : (myRank - 1);

  init(&params, &cells, myRank, numProcs, &av_vels);

  botHalo = &cells[0];
  botRow = &cells[params.nx];
  topRow = &cells[(params.heightH-2)*params.nx];
  topHalo = &cells[(params.heightH-1)*params.nx];

  ///////////////////////////////////////////////////////////////////////////////////////////

  /*
  halo_exchange(params, cells, up, down, status, botHalo, botRow, topRow, topHalo);
  //*/

  /* organised print
  if (myRank != (numProcs-1)){
    MPI_Recv(recvbuf, BUFSIZ, MPI_CHAR, up, 0, MPI_COMM_WORLD, &status);
  }
  print(&params, cells, myRank);
  sprintf(sendbuf, "Ok you can go now.");
  MPI_Ssend(sendbuf,strlen(sendbuf)+1, MPI_CHAR, down, 0, MPI_COMM_WORLD);
  if (myRank == (numProcs-1)){
    MPI_Recv(recvbuf, BUFSIZ, MPI_CHAR, up, 0, MPI_COMM_WORLD, &status);
  }
  //*/

  //*
  //set avvel numbers
  for (int i = 0; i < params.maxIters; i++){
    av_vels[i] = i * (myRank+1);
  }
  //*/

  if (myRank == MASTER){
    int height = params.ny / numProcs;
    int wrap = params.ny % numProcs;
    t_speed* all_cells_ptr = NULL; //pointer to start of array of cells
    t_speed* cell_insert = NULL; //pointer to current part of array where you want to read in data
    float* temp_vels = NULL;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*
    //allocate and clear the memory
    all_cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params.ny * params.nx));
    for (int i = 0; i < (params.ny * params.nx); i++){
      for (int j = 0; j < NSPEEDS; j++){
        all_cells_ptr[i].speeds[j] = 0.0f;
      }
    }

    //add master's portion of the grid
    cell_insert = all_cells_ptr;
    for (int i = 0; i < ((params.heightH-2) * params.nx); i++){
      cell_insert[i] = botRow[i];
    }
    cell_insert = &cell_insert[(params.heightH-2) * params.nx];

    //add every other rank's portion of the grid
    for (int rank = 1; rank < numProcs; rank++){
      int curHeight = height;
      if (rank < wrap){
        curHeight += 1;
      }
      MPI_Recv(cell_insert, (curHeight * params.nx * NSPEEDS), MPI_FLOAT, rank, 0, MPI_COMM_WORLD, &status);
      cell_insert = &cell_insert[curHeight * params.nx];
    }
    //*/
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*
    //allocate and clear memory for temp_vels
    temp_vels = malloc(sizeof(float) * params.maxIters);
    for (int i = 0; i < params.maxIters; i++){
      temp_vels[i] = 0.0f;
    }

    //receive av_vels from every rank and add them to master rank
    for (int rank = 1; rank < numProcs; rank++){
      MPI_Recv(temp_vels, params.maxIters, MPI_FLOAT, rank, 0, MPI_COMM_WORLD, &status);
      for (int i = 0; i < params.maxIters; i++){
        av_vels[i] += temp_vels[i];
      }
    }

    //free up memory again
    free(temp_vels);
    temp_vels = NULL;
    //*/
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*
    //print off grid
    for (int jj = params.ny-1; jj >= 0; jj--){
      for (int ii = 0; ii < params.nx; ii++){
        printf("(");
        for (int kk = 0; kk < NSPEEDS; kk++){
          float speed = all_cells_ptr[ii + jj*params.nx].speeds[kk];
          printf("%.2g, ", speed);
        }
        printf(") ");
      }
      printf("\n");
    }
    printf("\n\n");

    //print off av_vels
    for (int i = 0; i < params.maxIters; i++){
      printf("%.2g, ", av_vels[i]);
    }
    //*/
  }
  else{
    //*
    //send grid over
    MPI_Send(botRow, ((params.heightH-2) * params.nx * NSPEEDS), MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
    //send av_vels over
    MPI_Send(av_vels, params.maxIters, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
    //*/
  }

  MPI_Finalize();
  return 0;
}

/* More elegant halo exchange, 2 (large) messages sent and received per process per exchange
MPI_Sendrecv(topRow, NSPEEDS*NX, MPI_FLOAT, up, tag,
      botHalo, NSPEEDS*NX, MPI_FLOAT, down, tag, MPI_COMM_WORLD, &status);
MPI_Sendrecv(botRow, NSPEEDS*NX, MPI_FLOAT, down, tag,
      topHalo, NSPEEDS*NX, MPI_FLOAT, up, tag, MPI_COMM_WORLD, &status);
//*/

/////////////////////////////////////////////////////////////

/* Brute force halo exchange, 2 * NX * NSPEEDS messages sent and received per process per exchange
for (int i = 0; i < params.nx; i++){
  for (int j = 0; j < NSPEEDS; j++){
    //send top row receive bottom halo
    sprintf(sendbuf, "%f", cells[i+(params.heightH-2)*params.nx].speeds[j]); //create message
    MPI_Sendrecv(sendbuf, strlen(sendbuf)+1, MPI_CHAR, up, tag, //send the message
          recvbuf, BUFSIZ, MPI_CHAR, down, tag, MPI_COMM_WORLD, &status);
    cells[i].speeds[j] = atoi(recvbuf); //use received message
    //send bottom row receive top halo
    sprintf(sendbuf, "%f", cells[i+params.nx].speeds[j]); //create message
    MPI_Sendrecv(sendbuf, strlen(sendbuf)+1, MPI_CHAR, down, tag, //send the message
          recvbuf, BUFSIZ, MPI_CHAR, up, tag, MPI_COMM_WORLD, &status);
    cells[i+(params.heightH-1)*params.nx].speeds[j] = atoi(recvbuf); //use received message
  }
}
//*/

//*/

/*
//pack messages
position = 0;
MPI_Pack(&cells[params.nx].speeds[0], NSPEEDS*params.nx, MPI_FLOAT,
      bottomRow, BUFSIZ, &position, MPI_COMM_WORLD);


MPI_Sendrecv(bottomRow, strlen(bottomRow)+1, MPI_CHAR, down, tag, //send the message
      topRow, BUFSIZ, MPI_CHAR, up, tag, MPI_COMM_WORLD, &status);
//*/

/*
position = 0;
MPI_Unpack(topRow, BUFSIZ, &position, &topHalo[0], NSPEEDS*params.nx, MPI_FLOAT, MPI_COMM_WORLD);

if (myRank == MASTER){
  printf("\n");
  for (int i = 0; i < NSPEEDS*params.nx; i++){
    printf("%.2g, ", topHalo[i]);
  }
  printf("\n");
}

for (int i = 0; i < params.nx; i++){
  for (int j = 0; j < NSPEEDS; j++){
    cells[params.nx*(params.heightH-1)+i].speeds[j] = topHalo[i*NSPEEDS + j];
  }
}
*/

/*
MPI_Pack(&cells[params.nx*(params.heightH-2)].speeds[0], NSPEEDS*params.nx, MPI_FLOAT, topRow, BUFSIZ,
    &position, MPI_COMM_WORLD);

MPI_Sendrecv(bottomRow, strlen(sendbuf)+1, MPI_CHAR, down, tag, //send the message
      topHalo, BUFSIZ, MPI_CHAR, up, tag, MPI_COMM_WORLD, &status);

position = 0;
MPI_Unpack(topRow, BUFSIZ, &position, &bottomHalo[0], NSPEEDS*params.nx, MPI_FLOAT, MPI_COMM_WORLD);
*/

/*
if (myRank == MASTER){
  printf("\n");
  for (int i = 0; i < NSPEEDS*params.nx; i++){
    printf("%.2g, ", topHalo[i]);
  }
  printf("\n");
}

if (myRank == MASTER){
  printf("\n");
  for (int i = 0; i < NSPEEDS*params.nx; i++){
    printf("%.2g, ", bottomHalo[i]);
  }
  printf("\n");
}
*/

/*
int block_lengths[1];      // num of each type in a 'block' of the derived type
MPI_Aint displacements[1]; // associated memory displacements for each block
MPI_Datatype typelist[1];  // the actual intrinsic types comprising our bundle
MPI_Datatype halo;       // the new type we'll make
block_lengths[0] = NSPEEDS;
typelist[0] = MPI_FLOAT;
displacements[0] = 0;
MPI_Type_struct(1, block_lengths, displacements, typelist, &halo);
MPI_Type_commit(&halo);
//*/

/*
float *topRow, *topHalo, *botRow, *botHalo;
topHalo = &cells[(params.heightH-1)*params.nx].speeds[0];
topRow = &cells[(params.heightH-2)*params.nx].speeds[0];
botRow = &cells[params.nx].speeds[0];
botHalo = &cells[0].speeds[0];
*/
