
main(){
  //create a buffer for the reduction on both the host and the GPU
  float * reduction_buffer = malloc(sizeof(float));
  #pragma omp target enter data map(alloc: reduction_buffer[0:1])
  {}

  //main loop
  for (int itr = 0; itr < num_iterations; itr++){
    // Reset the reduction value to 0.0
    reduction_buffer[0] = 0.0f;
    #pragma omp target update to(reduction_buffer[0:1])
    {}

    // Execute vecadd on the target device
#pragma omp target teams distribute parallel for reduction(+:reduction_buffer[0])
    for (int i = 0; i < N; i++){
      reduction_buffer[0] += a[i] * b[i];
    }

    // Get the value from the reduction
    #pragma omp target update from(reduction_buffer[0:1])
    {}

    if (fabs(reduction_buffer[0] - (6.f * (float)N)) > 0.00001f) {
      printf("Incorrect answer at iteration %d %lf\n", itr, reduction_buffer[0]);
      correct_results = 0;
    }
  }

  //release the reduction buffer on the host and the device
  #pragma omp target exit data map(release: reduction_buffer[0:1])
  {}
  free(reduction_buffer);

  return 0;
}
