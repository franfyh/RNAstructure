/*
 * A program that folds two strands of nucleic acids in a duplex.
 * These strands of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef DUPLEXFOLD_H
#define DUPLEXFOLD_H

#include "../RNA_class/HybridRNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class DuplexFold {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	DuplexFold();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool parse( int argc, char** argv );

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	void run();

 private:
	// Private variables.

	// Description of the calculation type.
	string calcType;

	// Input and output file names.
	string seqFile1;         // The first input sequence file.
	string seqFile2;         // The second input sequence file.
	string ctFile;           // The output ct file.

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;

	// The maximum internal bulge loop size.
	int maxLoop;
	
	// The maximum number of structures to generate.
	int maxStructures;

	// The maximum percent energy difference.
	double percent;

	// The temperature at which calculation occurs.
	double temperature;
	
	// The window size for calculation.
	int windowSize;
};

#endif /* DUPLEXFOLD_H */
