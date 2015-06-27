#ifndef CONSERVED_TRAINING_H
#define CONSERVED_TRAINING_H

#include <dlib/optimization.h>

#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>

#include "../src/defines.h"
#include "../RNA_class/RNA.h"
#include "../src/phmm/utils/xmath/matrix/matrix.h"
#include "../src/configfile.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

using namespace std;
//using namespace dlib;
typedef dlib::matrix<double,0,1> column_vector;

class conserved_training_interface {
public:
    conserved_training_interface();

    bool parse(int argc, char** argv);
    friend double conserved_accuracy_wrapper(const column_vector& m);
    friend const column_vector conserved_accuracy_derivative_wrapper(const column_vector& m);
    vector<double> parameters;
    double conserved_testing_accuracy(const column_vector& m);

private:

    vector<string> training_files;
    vector<string> testing_files;
    string conf_file;
//    int processors;


    double alpha;
    int training;
    double conserved_accuracy(const column_vector& m);
    const column_vector conserved_accuracy_derivative(const column_vector& m);
};

#endif
