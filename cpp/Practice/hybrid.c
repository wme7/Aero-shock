/*
   This is a hello world program utilizing both MPI and OpenMP.

   In order to coordinate output, all output is handled by the master
   process.  Within the master process, first, each thread says hello.
   Once this is completed, the master thread waits for MPI sends from
   each of the other processes.  The first piece of data is how many
   threads the process has.  This is sent by the master thread of the
   remote process.  Then, each thread will send a thread ID, process
   rank, and processor name to the master process.  This will then be
   formatted and sent to standard output as a hello from the sending
   thread.
*/

// Include the MPI header <mpi.h> and the OpenMP header <omp.h>
// The MPI header should be included before stdio.h.

#include <mpi.h>
#include <omp.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

   int rank;                           // Rank ID of the current process
   int nproc;                          // Total number of processes
   int nthreads;                       // Total number of threads
   int threadID;                       // ID of the current thread
   int namelen;                        // Length of the processor name
   int required=MPI_THREAD_SERIALIZED; // Required level of MPI threading support
   /* Each thread will call MPI routines, but these calls will be coordinated
      to occur only one at a time within a process.
   */
   int provided;                       // Provided level of MPI threading support
   char name[MPI_MAX_PROCESSOR_NAME];  // Name of the processor
   int dThread;                        // Display thread ID
   int dRank;                          // Display rank ID
   int dNamelen;                       // Length of display name
   char dName[MPI_MAX_PROCESSOR_NAME]; // Display processor name
   int sNthreads;                      // nthreads from sender
   MPI_Status stat;                    // Status from MPI calls
   int r;                              // Rank loop counter
   int t;                              // Thread loop counter

   // Initialize MPI with threading

   MPI_Init_thread(&argc, &argv, required, &provided);

   // Determine the MPI rank, number of processes, and processor name

   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Get_processor_name(name,&namelen);

   // Check the threading support level

   if (provided < required) {
   
      // Insufficient support, degrade to 1 thread and warn the user
      
      if (rank == 0) {
         printf("Warning:  This MPI implementation provides insufficient");
         printf(" threading support.\n");
      }
      omp_set_num_threads(1);
   }

   // The multithreaded section where all threads will say hello
      
   #pragma omp parallel default(shared) private(threadID)
   {

   // All processes should get the total number of threads, each
   // threads needs to know its own ID.
      
   threadID=omp_get_thread_num();   // Get the thread ID
   nthreads=omp_get_num_threads();  // Get the total number of threads

   // Time to say hello, the master process performs all output.
   // Within the master process, each thread will handle its own
   // output, the master thread will handle output from all threads
   // of all other processes.

   if (rank == 0) {
   
      // The master process outputs from its own threads
      // This section is done by every OpenMP thread, but only one at a time.
      // This requires MPI_THREAD_SERIALIZED.
      
      #pragma omp critical
      {
         printf("Hello from thread %d of %d in rank %d of %d on %s.\n",
            threadID, nthreads, rank, nproc, name);
      }  // End of #pragma omp critical
      
      #pragma omp barrier
      
      // Now, receive data from each of the other processes and
      // give an appropriate greeting.  Only the master thread
      // should do this.  Since only the master thread is calling
      // MPI, this is an example of MPI_THREAD_FUNNELED.
      
      #pragma omp master
      {
      
      for (r=1;r<nproc;r++) {
      
         // Get the number of threads in the sender
         
         MPI_Recv(&sNthreads, 1, MPI_INT, r, 10*r, MPI_COMM_WORLD, &stat);
         
         for (t=0;t<sNthreads;t++) {
         
            // For each thread, get the rank ID, thread ID, and name
            
            MPI_Recv(&dRank, 1, MPI_INT, r, 10*r+1, MPI_COMM_WORLD, &stat);
            MPI_Recv(&dThread, 1, MPI_INT, r, 10*r+2, MPI_COMM_WORLD, &stat);
            MPI_Recv(&dNamelen, 1, MPI_INT, r, 1000*r+10*dThread, MPI_COMM_WORLD, &stat);
            MPI_Recv(dName, dNamelen+1, MPI_CHAR, r, 1000*r+10*dThread+1, MPI_COMM_WORLD, &stat);
            
            printf("Hello from thread %d of %d in rank %d of %d on %s.\n",
               dThread, sNthreads, dRank, nproc, dName);
         }
      }
      
      }  // End of #pragma omp master
      
   } else { // All other processes will send their data to the master
   
      // Only the master sends the number of threads.  MPI_THREAD_FUNNELED
   
      #pragma omp master
      {
      MPI_Send(&nthreads, 1, MPI_INT, 0, 10*rank, MPI_COMM_WORLD);
      }  // End of #pragma omp master
      
      #pragma omp critical
      {
      
      // Each thread will send its own data, but there is no
      // particular order required, so a critical section works
      // exactly as needed.  As such, this requires MPI_THREAD_SERIALIZED
      
      MPI_Send(&rank, 1, MPI_INT, 0, 10*rank+1, MPI_COMM_WORLD);
      MPI_Send(&threadID, 1, MPI_INT, 0, 10*rank+2, MPI_COMM_WORLD);
      MPI_Send(&namelen, 1, MPI_INT, 0, 1000*rank+10*threadID, MPI_COMM_WORLD);
      MPI_Send(name, namelen+1, MPI_CHAR, 0, 1000*rank+10*threadID+1, MPI_COMM_WORLD);
      
      }  // End of #pragma omp critical

   }
   
   }  // End of #pragma omp parallel
   
   // Close out MPI and the program

   MPI_Finalize();

   return 0;
}