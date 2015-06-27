#include "../phmm/phmm_interface.h"
#include "../RNA_class/TwoRNA.h"
#include "../src/structure.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
phmm_interface::phmm_interface() {

	//  Initialize the maximum likelihood variable
	bool ML = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool phmm_interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "phmm" );
	parser->addParameterDescription( "seq 1", "The name of a file containing the first input sequence." );
	parser->addParameterDescription( "seq 2", "The name of a file containing the second input sequence." );
	parser->addParameterDescription( "out file", "The name of a file containing the output sequence." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the Maximum Likelihood option.
	vector<string> alignmentOptions;
	alignmentOptions.push_back( "-M" );
	alignmentOptions.push_back( "-m" );
	alignmentOptions.push_back( "--ML" );
	parser->addOptionFlagsNoParameters( alignmentOptions, "Specify that program should output a maximum likelihood alignment. Default is to output pairwise probabilities." );


	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seq1 = parser->getParameter( 1 );
		seq2 = parser->getParameter( 2 );
		outfile = parser->getParameter( 3 );
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the DNA option.
	if( !parser->isError() ) { isRNA = !parser->contains( dnaOptions ); }

	// Get the maximum likelihood option.
	if( !parser->isError() ) { ML = parser->contains( alignmentOptions ); }



	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void phmm_interface::run() {
	std::cout << "Running phmm!\n";
	structure ct1, ct2;
	int l1;
	int l2;
	int error = 0;
	if(!ct1.openseq(seq1.c_str())) {cerr << "ERROR: Could not open sequence file "<<seq1<<"\n"; error = 1;}
	if(!ct2.openseq(seq2.c_str())) {cerr << "ERROR: Could not open sequence file "<<seq2<<"\n"; error = 1;}

	if(error == 0){
		t_phmm_aln* phmm_aln = create_phmm_aln(&ct1,&ct2);						//t_phmm_aln class contains the two sequences and the hidden markov model
		l1 = phmm_aln->l1();										//get the length of the two sequences, used to write the output file
		l2 = phmm_aln->l2();
		if(!ML) {write_probability_array(phmm_aln->compute_posterior_probs(), outfile.c_str(),l1,l2);} 	//pairwise probability calculation (default behavior)
		else {write_ML_alignment(phmm_aln->compute_ML_alignment(), outfile.c_str(),l1,l2,seq1.c_str(),seq2.c_str());}		//maximum likelihood alignment calculation
		std::cout << "All done!\n";
	}
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	phmm_interface* runner = new phmm_interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
