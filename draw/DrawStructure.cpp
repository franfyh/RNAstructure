/*
 * A program that draws a structure and writes image output.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "DrawStructure.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DrawStructure::DrawStructure() {

	// Initialize the calculation type description.
	calcType = "Structure drawing";

	// Initialize the circular structure boolean to false.
	circular = false;

	// Initialize the encircling nucleotides boolean to true.
	encircle = true;

	// Initialize the linear structure boolean to false.
	flat = false;

	// Initialize the SVG boolean to false.
	isSVG = false;

	// Initialize the levorotatory flag to false.
	levorotatory = false;

	// Initialize the output structure number to -1, which signifies that all structures should be output.
	number = -1;

	// Initialize the text annotation flag to false.
	textAnnotation = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool DrawStructure::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "draw" );
	parser->addParameterDescription( "ct file", "The name of a file containing CT data for the structure to be drawn." );
	parser->addParameterDescription( "output file", "The name of an image file to which output will be written. Usually, this is a Postscript image file, although the user can specify that it be an SVG image file instead. To specify SVG images, the \"--svg\" flag must be specified in conjunction with a structure number flag." );

	// Add the circular option.
	vector<string> circleOptions;
	circleOptions.push_back( "-c" );
	circleOptions.push_back( "-C" );
	circleOptions.push_back( "--circle" );
	parser->addOptionFlagsNoParameters( circleOptions, "Specify that the structure should be drawn with its backbone stretched around a circle. Note that pseudoknotted structures will be drawn circularized even if this option is not specified. Default is to show a collapsed structure." );

	// Add the flat option.
	vector<string> flatOptions;
	flatOptions.push_back( "-f" );
	flatOptions.push_back( "-F" );
	flatOptions.push_back( "--flat" );
	parser->addOptionFlagsNoParameters( flatOptions, "Specify that the structure should be drawn with its backbone stretched in a straight line. Default is to show a collapsed structure." );

	// Add the levorotatory option.
	vector<string> levorotatoryOptions;
	levorotatoryOptions.push_back( "-l" );
	levorotatoryOptions.push_back( "-L" );
	levorotatoryOptions.push_back( "--levorotatory" );
	parser->addOptionFlagsNoParameters( levorotatoryOptions, "Specify that the drawn structure is rendered counterclockwise. Default is to render drawn structures clockwise." );

	// Add the number option.
	vector<string> numberOptions;
	numberOptions.push_back( "-n" );
	numberOptions.push_back( "-N" );
	numberOptions.push_back( "--number" );
	parser->addOptionFlagsWithParameters( numberOptions, "Specify the index of a particular structure in the predicted CT to be compared with the accepted CT, one-indexed. Default is -1, which signifies all structures output to one file." );

	// Add the probability option.
	vector<string> probabilityOptions;
	probabilityOptions.push_back( "-p" );
	probabilityOptions.push_back( "-P" );
	probabilityOptions.push_back( "--probability" );
	parser->addOptionFlagsWithParameters( probabilityOptions, "Specify the name of the partition function file from which base pairing probability data will be read for annotation. This file should describe pairing data for the predicted structure, not the accepted structure. Default is no probability annotation file used." );

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-s" );
	shapeOptions.push_back( "-S" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify the name of the file from which SHAPE data will be read for annotation. Default is no SHAPE annotation file used." );

	// Add the SVG option.
	vector<string> svgOptions;
	svgOptions.push_back( "--svg" );
	parser->addOptionFlagsNoParameters( svgOptions, "Specify that the output file should be an SVG image file, rather than a Postscript image file. Note that only one SVG image can be written into a particular file, so the structure number flag must also be specified when writing an SVG document." );

	// Add the text annotation option.
	vector<string> textOptions;
	textOptions.push_back( "-t" );
	textOptions.push_back( "-T" );
	textOptions.push_back( "--text" );
	parser->addOptionFlagsWithParameters( textOptions, "Specify the name of the text file from which base pairing probability data will be read for annotation. This file should describe pairing data for the predicted structure, not the accepted structure. Default is no probability annotation file used." );

	// Add the uncircled option.
	vector<string> uncircledOptions;
	uncircledOptions.push_back( "-u" );
	uncircledOptions.push_back( "-U" );
	uncircledOptions.push_back( "--uncircled" );
	parser->addOptionFlagsNoParameters( uncircledOptions, "Specify that no circles should surround nucleotides when drawing. Default is to surround nucleotides with circles." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		inputFile = parser->getParameter( 1 );
		outputFile = parser->getParameter( 2 );
	}

	// Get the circular option.
	if( !parser->isError() ) { circular = parser->contains( circleOptions ); }

	// Get the flat option.
	if( !parser->isError() ) { flat = parser->contains( flatOptions ); }

	// Get the levorotatory option.
	if( !parser->isError() ) { levorotatory = parser->contains( levorotatoryOptions ); }

	// Get the number option.
	if( !parser->isError() ) {
		parser->setOptionInteger( numberOptions, number );
		bool badNumber =
			( number != -1 ) &&
			( number < 0 );
		if( badNumber ) { parser->setError( "structure number" ); }
	}

	// Get the probability annotation file option.
	// This can be read either as a partition function save file or a text file.
	if( !parser->isError() ) {
		probabilityFile = parser->getOptionString( probabilityOptions, true );
		if( probabilityFile == "" ) {
			probabilityFile = parser->getOptionString( textOptions, true );
			if( probabilityFile != "" ) { textAnnotation = true; }
		}
	}

	// Get the SHAPE annotation file option.
	if( !parser->isError() ) { SHAPEFile = parser->getOptionString( shapeOptions, true ); }

	// Get the SVG option.
	if( !parser->isError() ) { isSVG = parser->contains( svgOptions ); }

	// Get the uncircled option.
	if( !parser->isError() ) { encircle = !parser->contains( uncircledOptions ); }

	// If the user requests levorotatory and linear, show an error, because linear structures should only be shown one way.
	if( levorotatory && flat ) {
		parser->setErrorSpecialized( "Linear structures cannot be rendered counterclockwise." );
	}

	// If in SVG mode and no structure number is specified, show an error message.
	if( ( !parser->isError() ) && parser->contains( svgOptions ) && ( !parser->contains( numberOptions ) ) ) {
		cerr << "No structure number specified with SVG structure." << endl
		     << "Please specify the structure to draw as an SVG image file." << endl;
		parser->setError();
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void DrawStructure::run() {

	// Create a variable that tracks the result of calculations.
	// Throughout, the calculation continues if result = "".
	string result = "";

	// Create variables to determine the range of structures that will be drawn.
	// Initialize them so the range of structures will only be the specific structure given, if any.
	int start = number;
	int end = number;

	// Create a temporary strand for the purpose of determining the number and range of structures to be drawn.
	int structures = 0;
	if( result == "" ) {

		// Show a message that the information gathering phase has begun.
		cout << "Determining structure information..." << flush;

		// Create the temporary strand to get the number of structures, then delete it.
		RNA* strand = new RNA( inputFile.c_str(), 1 );
		ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
		result = checker->returnError();
		structures = strand->GetStructureNumber();
		delete strand;
		delete checker;

		// Use the number of structures to determine if the user gave an appropriate specific structure value, if necessary.
		// If the structure value is invalid, set the error string to show this.
		if( result == "" ) {
			if( number != -1 && !( number >= 1 && number <= structures ) ) {
				result = "Structure number out of range.\n";
			}
		}

		// If the user did not give a specific structure, set the range of structures to all possible structures.
		if( number == -1 ) {
			start = 1;
			end = structures;
		}

		// If no error occurred, print a message saying that the information gathering phase is done.
		// Otherwise, print the error message.
		if( result == "" ) { cout << "done." << endl; }
		else { cerr << endl << result; }
	}

	// If no errors have occurred in structure preparation, draw the structure(s).
	if( result == "" ) {

		// Print a message saying that the drawing phase has started.
		cout << "Drawing..." << endl;

		// If the output file already exists, delete the file so a new file with that name can be generated.
		ifstream test( outputFile.c_str() );
		bool exists = test.good();
		test.close();
		if( exists ) { remove( outputFile.c_str() ); }

		// For each structure in the range, get its data and draw it appropriately.
		for( int i = start; i <= end; i++ ) {

			// Print a message saying that a particular structure has started drawing.
			cout << "    Drawing structure " << i << "..." << flush;

			// Create the drawing handler and set whether it should circle nucleotides.
			StructureImageHandler* handler = new StructureImageHandler();
			handler->setNucleotidesCircled( encircle );

			// Read in the structure.
			// Note that if a pseudoknotted structure is fed to the radial routine, the radial routine automatically calls the circular routine.
			if( result == "" ) {
				if( circular ) { result = handler->readCircular( inputFile, i ); }
				else if( flat ) { result = handler->readLinear( inputFile, i ); }
				else { handler->readRadial( inputFile, i ); }
			}

			// Read in annotation, if necessary.
			if( result == "" ) {
				if( probabilityFile != "" ) {
					result = handler->addAnnotationProbability( probabilityFile, textAnnotation );
				} else if( SHAPEFile != "" ) {
					result = handler->addAnnotationSHAPE( SHAPEFile );
				}
			}

			// Flip the structure, if asked.
			if( ( result == "" ) && ( levorotatory == true ) ) { handler->flipHorizontally(); }

			// Write the structure image, as either SVG or Postscript.
			if( result == "" ) {
				if( isSVG ) { handler->writeSVG( outputFile ); }
				else {
					bool append = ( number == -1 );
					if( i == start ) {
						ifstream test( outputFile.c_str() );
						if( test ) { append = false; }
					}
					handler->writePostscript( outputFile, append );
				}
			}

			// Delete the drawing handler.
			delete handler;

			// Print a message saying that a particular structure has finished drawing, if no error occurred.
			// If an error did occur, break out of the loop.
			if( result == "" ) { cout << "done." << endl; }
			else break;
		}

		// If no error occurred, print a message saying that the drawing phase is done.
		// Otherwise, print the error message.
		if( result == "" ) { cout << "done." << endl; }
		else { cerr << endl << result << endl; }
	}

	// Print confirmation of run finishing.
	if( result == "" ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	DrawStructure* runner = new DrawStructure();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
