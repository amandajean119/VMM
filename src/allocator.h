/*
 * allocator.h
 *
 *  Created on: Apr 3, 2013
 *  Author:
 */

#include <stdint.h>

#ifndef ALLOCATOR_H_
#define ALLOCATOR_H_

typedef uint64_t timeUs;

struct vcore;
/*
 * Physical core data structure
 */

typedef struct pcore{

  int id;
  int speedKhz;			// pcpu speed in Khz
  float maxUtilization;		// Maximum utilization (Range 0 - 100)
  float utilization;		// Achieved utilization
  int nrVCores;
  struct vcore * vcores;
} pcore;

/*
 * Virtual core data structure
 */

typedef struct vcore{

  int id;
  timeUs slice;	   // Size of the slice in usec
  timeUs period;   // Size of the period in usec
  int speedKhz;	   // vcpu speed in Khz
  pcore * pcore;   // Physical core allocated to this vcore
} vcore;


/*
 * Virtual Machine data structure
 */

typedef struct vm {

  int nrVcores;	   // Number of virtual cores in this virtual machine
  int tdf;         // Time dilation factor for this vm
  vcore * vcores;  // Pointer to the array of virtual cores of this vm
} vm;


/*
 * Physical Machine data structure
 */

typedef struct machine {

  int id;
  int nrPcores;
  pcore * pcores;       // Array of physical cores in this physical machine
  vm * vm;		        // Pointer to the vm
  float fitness;	    // Fitness function for this physical machine
    			        // configuration
} machine;


/*
 * Allocator configuration data structure
 */

typedef struct allocatorConfig{

  timeUs minSlice;       // Minimum slice allowed for a vcore
  timeUs maxSlice;       // Maximum slice allowed for a vcore
  timeUs minPeriod;      // Minimum period allowed for a vcore
  timeUs maxPeriod;      // Maximum period allowed for a vcore
  int minTdf;		 // Minimum time dilation factor allowed for a vm
  int maxTdf;		 // Maximum time dilation factor allowed for a vm
  int minCPU;            // Minimum CPU speed
  int maxCPU;   	 // Maximum CPU speed
  int minU;		 // Minimum utilization factor
  int maxU;		 // Maximum utilization factor

} allocConfig;


/*
 * General Data structure for the Genetic algorithm
 */

typedef struct ga {

  int populationSize;
  machine * population;  // Array with the different configurations
			 // for a physical machine. Each element of the
			 // array represents an individual of the
			 // population
  allocConfig * aConfig;	// Allocator Configuration
  float bestFitness;
  int bestMachineIndex;
  float avgFitness;
} ga;


#endif /* ALLOCATOR_H_ */
