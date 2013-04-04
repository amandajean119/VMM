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

/*
 * General Data structure for the Genetic algorithm
 */

struct ga {

	struct machine * population; // Array with the different configurations
								 // for a physical machine. Each element of the
								 // array represents an individual of the
								 // population
    struct allocatorConfig * aConfig;	// Allocator Configuration
};


/*
 * Allocator configuration data structure
 */

struct allocatorConfig{

	 timeUs minSlice;       // Minimum slice allowed for a vcore
	 timeUs maxSlice;       // Maximum slice allowed for a vcore
	 timeUs minPeriod;      // Minimum period allowed for a vcore
	 timeUs maxPeriod;      // Maximum period allowed for a vcore
	 int minTdf;			 // Minimum time dilation factor allowed for a vm
	 int maxTdf;		     // Maximum time dilation factor allowed for a vm


};


/*
 * Physical Machine data structure
 */

struct machine {

	struct pcore * pcores;  // Array of physical cores in this physical machine
    struct vm * vms;		// Pointer to the arrray of vms of this machine
    int fitness;			// Fitness function for this physical machine
    						// configuration
};


/*
 * Virtual Machine data structure
 */

struct vm {

	int nrVcores;			// Number of virtual cores in this virtual machine
	int tdf;                // Time dilation factor for this vm
	struct vcore * vcores;  // Pointer to the array of virtual cores of this
							// vm
};

/*
 * Virtual core data structure
 */

struct vcore{

	timeUs slice;			// Size of the slice in usec
	timeUs period;			// Size of the period in usec
	int speedKhz;			// vcpu speed in Khz
	struct phCore * pcore;  // Physical core allocated to this vcore
};

/*
 * Physical core data structure
 */

struct pcore{

	int speedKhz;			// pcpu speed in Khz
    int maxUtilization;		// Maximum utilization (Range 0 - 100)
    int utilization;		// Achieved utilization
};



#endif /* ALLOCATOR_H_ */
