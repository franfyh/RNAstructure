
#include "conserved_training.h"

conserved_training_interface::conserved_training_interface() {

//    processors = 1;
    parameters.assign(8,0);
    training = 0;
    alpha = 0;
}

bool conserved_training_interface::parse(int argc, char** argv){

    string type = "conserved_training";

    ParseCommandLine* parser = new ParseCommandLine(type);
    parser->addParameterDescription("configuration file", "The name of a file containing configuration data." );

    parser->parseLine(argc,argv);

    ConfigFile file(parser->getParameter(1) );
    conf_file = parser->getParameter(1);
    if(file.isValid() == false){
        delete parser;
        return false;
    }

    bool isReadableGroups =
        file.contains("TrainingFiles");
    if(isReadableGroups){
        string fileData = file.getOption<string>( "TrainingFiles" );
		unsigned int fileLast = fileData.length() - 1;
		if( fileData[0] == '{' && fileData[fileLast] == '}' ) {
			fileData = fileData.erase( 0, 1 );
			fileData = fileData.erase( fileLast - 1, 1 );
			fileLast = fileData.length() - 1;
			if( fileData[fileLast] == ';' ) { fileData = fileData.erase( fileLast, 1 ); }
			stringstream fileStr( fileData );
			string fileFile;
			while( fileStr.good() ) {
				getline( fileStr, fileFile, ';' );
				if( fileFile != "" ) { training_files.push_back( fileFile ); }
			}
		}  else {
			if( fileData[0] != '{' ) { parser->setErrorSpecialized( "File group has no start bracket." ); }
			else { parser->setErrorSpecialized( "File group has no end bracket." ); }
			delete parser;
			return false;
		}
                
                if(training_files.size() % 3 != 0){
			parser->setErrorSpecialized( "Number of sequence files does not divide three." );
			delete parser;
			return false;
		}
                

    }

    isReadableGroups =
        file.contains("TestingFiles");
   
    if(isReadableGroups){
        string fileData = file.getOption<string>( "TestingFiles" );
		unsigned int fileLast = fileData.length() - 1;
		if( fileData[0] == '{' && fileData[fileLast] == '}' ) {
			fileData = fileData.erase( 0, 1 );
			fileData = fileData.erase( fileLast - 1, 1 );
			fileLast = fileData.length() - 1;
			if( fileData[fileLast] == ';' ) { fileData = fileData.erase( fileLast, 1 ); }
			stringstream fileStr( fileData );
			string fileFile;
			while( fileStr.good() ) {
				getline( fileStr, fileFile, ';' );
				if( fileFile != "" ) { testing_files.push_back( fileFile ); }
			}
		}  else {
			if( fileData[0] != '{' ) { parser->setErrorSpecialized( "File group has no start bracket." ); }
			else { parser->setErrorSpecialized( "File group has no end bracket." ); }
			delete parser;
			return false;
		}
                
                if(testing_files.size() % 3 != 0){
			parser->setErrorSpecialized( "Number of sequence files does not divide three." );
			delete parser;
			return false;
		}
                

    }


    for(int i = 0; i < 8; ++i)
    {
        stringstream numStream( stringstream::in | stringstream::out );
        numStream << i;
        string num = numStream.str();

        string parName = "Parameter" + num;
        if(file.contains( parName )) {
            parameters[i] = file.getOption<double>(parName);
        }
    }  

    if( !parser->isError() ) {
        if( file.contains( "Alpha" ) ) {
            alpha = file.getOption<double>( "Alpha" );

            if( alpha < 0.0 ) { parser->setError( "Alpha" );
            }
            cout << "Alpha: "<<alpha<<"\n";
        }
    }

                
    //running~~~

    bool noError = ( parser->isError() == false );
    delete parser;
    return noError;

}

double conserved_training_interface::conserved_accuracy(const column_vector& m)
{
    vector<double> log_accuracy(training,0.0);
    double sum_log_accuracy = 0.0;
    for(int i=0;i<7;++i){
//        cerr <<"i "<<i<<" "<<0.5*alpha*m(i)*m(i)<<"\n";
//        cerr <<exp(m(i) )<<" "<<alpha<<"\n";
    
    }
    for(int i = 0; i< training_files.size()/3; ++i) {
        vector<string> temp_files;
        temp_files.push_back(training_files[i*3]);
        temp_files.push_back(training_files[i*3+1]);
        temp_files.push_back(training_files[i*3+2]);
        cerr << temp_files[0] << "\n";
        ifstream myfile(temp_files[0].c_str() );
        string temp;
        string sequence1;
        string sequence2;
        getline(myfile,temp);
        getline(myfile,sequence1);
        getline(myfile,temp);
        getline(myfile,sequence2);
//        sequence1 = sequence1.substr(0,sequence1.size()-2);
//        sequence2 = sequence2.substr(0,sequence2.size()-2);
        string target_sequence = "";
        map<int,int> alignment_map;
        int current_position1 = 0;
        int current_position2 = 0;
        for(int j = 0;j<sequence1.size();++j) {
            if(sequence2[j]!='-') ++current_position2;
            if(sequence1[j]!='-') {
                target_sequence.push_back(sequence1[j]);
                ++current_position1;
                if (sequence2[j]!='-'){alignment_map[current_position1] = current_position2;}
                else{alignment_map[current_position1] = 0;}
                    
                }
                
        }
        
        t_matrix* basepairing_extrinsic_info = new t_matrix(current_position1, current_position1, true);
        vector<double>* singlestranded_extrinsic_info = new vector<double>(current_position1+1,0);
        RNA* structure1 = new RNA(temp_files[1].c_str(),1,true);
        RNA* structure2 = new RNA(temp_files[2].c_str(),1,true);
        RNA* partition_function_target = new RNA(target_sequence.c_str() );
        for(int k =1;k<= current_position1;++k){

            for(int l=k+1;l<= current_position1;++l){
                if(alignment_map[k]==0 && alignment_map[l]==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(2) );}
                else if(alignment_map[k]!=0 && alignment_map[l]==0 && structure2->GetPair(alignment_map[k],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(6) );}
                else if(alignment_map[k]==0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[l],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(6) );}
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)==0 && structure2->GetPair(alignment_map[l],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(5) );}
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)!=0 && structure2->GetPair(alignment_map[l],1)!=0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(3) );}
                else {
//                    cerr << "something wrong with base pair "<<k<<" "<<l<<" in sequence pair "<<temp_files[0]<<"\n";
//                    exit(EXIT_FAILURE);
                    basepairing_extrinsic_info->x(k,l) = 0;

                }


           }
            
            if(alignment_map[k]==0){
                singlestranded_extrinsic_info->at(k) = exp(m(0) );
            }
            else if(alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)!=0){
                singlestranded_extrinsic_info->at(k) = exp(m(4) );
        }
        else if (alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)==0){
            singlestranded_extrinsic_info->at(k) = exp(m(1) );
        }
        else {
            cerr << "something wrong with nucleotide "<<k<<" in sequence pair "<<temp_files[0]<<"\n";
            exit(EXIT_FAILURE);}

            

        }
        double max_extrinsic_info = 0;
        for(int k =1;k<= current_position1;++k){
            if (singlestranded_extrinsic_info->at(k) > max_extrinsic_info)
                max_extrinsic_info = singlestranded_extrinsic_info->at(k);
        }
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                if(basepairing_extrinsic_info->x(k,l) > max_extrinsic_info*max_extrinsic_info)
                    max_extrinsic_info = sqrt(basepairing_extrinsic_info->x(k,l));
            }
        }
        
        for(int k =1;k<= current_position1;++k){
            singlestranded_extrinsic_info->at(k) /= max_extrinsic_info;
//            cerr << "k "<<singlestranded_extrinsic_info->at(k)<<"\n";
        }
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                basepairing_extrinsic_info->x(k,l) /= max_extrinsic_info*max_extrinsic_info;
                //              cerr <<"k l "<<basepairing_extrinsic_info->x(k,l)<<"\n";
            }
        }
        
        



        partition_function_target->SetSinglestranded( singlestranded_extrinsic_info);
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                partition_function_target->SetExtrinsic(k,l,basepairing_extrinsic_info->x(k,l) );

            }
        }
//        cerr << "target ";
        partition_function_target->PartitionFunction();
        double energy_partition_function = partition_function_target->GetEnsembleEnergy();

        RNA* equilibrium_constant_target = new RNA(target_sequence.c_str() );
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                if (structure1->GetPair(k,1)!=l){
                    basepairing_extrinsic_info->x(k,l) = 0;
                }
                equilibrium_constant_target->SetExtrinsic(k,l,basepairing_extrinsic_info->x(k,l) );

            }
            if (structure1 ->GetPair(k,1)!=0){
                singlestranded_extrinsic_info->at(k) = 0;
            }
        }
        
        equilibrium_constant_target->SetSinglestranded( singlestranded_extrinsic_info);
//        cerr << "constant ";
        equilibrium_constant_target->PartitionFunction();
        double energy_equilibrium_constant = equilibrium_constant_target->GetEnsembleEnergy();
        
       sum_log_accuracy += energy_partition_function/(-1*RKC*310.15);
//        cerr <<"accuracy "<<log_accuracy[i]<<"\n";
        sum_log_accuracy -= energy_equilibrium_constant/(-1*RKC*310.15);
        //       cerr <<"accuracy "<<log_accuracy[i]<<"\n";
        cerr << "partition "<<energy_partition_function<<"\n";
        cerr << "eq "<<energy_equilibrium_constant<<"\n";

        delete structure1;
        delete structure2;
        delete basepairing_extrinsic_info;
        delete singlestranded_extrinsic_info;
        delete partition_function_target;
        delete equilibrium_constant_target;
    
   }
    sum_log_accuracy /= training_files.size();

    for(int i=0;i<7;++i){
//        cerr <<"i "<<i<<" "<<0.5*alpha*m(i)*m(i)<<"\n";
//        cerr <<m(i)<<" "<<alpha<<"\n";
        sum_log_accuracy += 0.5*alpha*m(i)*m(i);

    }
//    cerr <<"final "<<sum_log_accuracy<<"\n";
 
    return sum_log_accuracy;
}

const column_vector conserved_training_interface::conserved_accuracy_derivative(const column_vector& m){
    column_vector res(7);
    for(int p=0;p<7;++p)
        res(p) = 0; 


    for(int i = 0; i< training_files.size()/3; ++i) {
          vector<string> temp_files;
        temp_files.push_back(training_files[i*3]);
        temp_files.push_back(training_files[i*3+1]);
        temp_files.push_back(training_files[i*3+2]);
          ifstream myfile(temp_files[0].c_str() );
        string temp;
        string sequence1;
        string sequence2;
        getline(myfile,temp);
        getline(myfile,sequence1);
        getline(myfile,temp);
        getline(myfile,sequence2);
//      sequence1 = sequence1.substr(0,sequence1.size()-2);
//      sequence2 = sequence2.substr(0,sequence2.size()-2);
//        cerr <<"sequence1 "<<sequence1<<"\n";
//        cerr <<"sequence2 "<<sequence2<<"\n";
        string target_sequence = "";
        map<int,int> alignment_map;
        int current_position1 = 0;
        int current_position2 = 0;
        for(int j = 0;j<sequence1.size();++j) {
            if(sequence2[j]!='-') ++current_position2;
            if(sequence1[j]!='-') {
                target_sequence.push_back(sequence1[j]);
                ++current_position1;
                if (sequence2[j]!='-'){alignment_map[current_position1] = current_position2;}
                else{alignment_map[current_position1] = 0;}
                    
                }
                
        }
        
        t_matrix* basepairing_extrinsic_info = new t_matrix(current_position1, current_position1, true);
        vector<double>* singlestranded_extrinsic_info = new vector<double>(current_position1+1,0);
        RNA* structure1 = new RNA(temp_files[1].c_str(),1,true);
        RNA* structure2 = new RNA(temp_files[2].c_str(),1,true);
        RNA* partition_function_target = new RNA(target_sequence.c_str() );
        for(int k =1;k<= current_position1;++k){

            for(int l=k+1;l<= current_position1;++l){
                if(alignment_map[k]==0 && alignment_map[l]==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(2) );}
                else if(alignment_map[k]!=0 && alignment_map[l]==0 && structure2->GetPair(alignment_map[k],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(6) );}
                else if(alignment_map[k]==0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[l],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(6) );}
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)==0 && structure2->GetPair(alignment_map[l],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(5) );}
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)!=0 && structure2->GetPair(alignment_map[l],1)!=0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(3) );}
                else {
//                    cerr << "something wrong with base pair "<<k<<" "<<l<<" in sequence pair "<<temp_files[0]<<"\n";
//                    exit(EXIT_FAILURE);
                    basepairing_extrinsic_info->x(k,l) = 0;

}


           }
            
            if(alignment_map[k]==0){
                singlestranded_extrinsic_info->at(k) = exp(m(0) );
            }
            else if(alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)!=0){
                singlestranded_extrinsic_info->at(k) = exp(m(4) );
        }
        else if (alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)==0){
            singlestranded_extrinsic_info->at(k) = exp(m(1) );
        }
        else {
            cerr << "something wrong with nucleotide "<<k<<" in sequence pair "<<temp_files[0]<<"\n";
            exit(EXIT_FAILURE);}

            

        }
                

         double max_extrinsic_info = 0;
        for(int k =1;k<= current_position1;++k){
            if (singlestranded_extrinsic_info->at(k) > max_extrinsic_info)
                max_extrinsic_info = singlestranded_extrinsic_info->at(k);
        }
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                if(basepairing_extrinsic_info->x(k,l) > max_extrinsic_info*max_extrinsic_info)
                    max_extrinsic_info = sqrt(basepairing_extrinsic_info->x(k,l));
            }
        }
        
        for(int k =1;k<= current_position1;++k){
            singlestranded_extrinsic_info->at(k) /= max_extrinsic_info;
        }
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                basepairing_extrinsic_info->x(k,l) /= max_extrinsic_info*max_extrinsic_info;
            }
        }
        


        partition_function_target->SetSinglestranded( singlestranded_extrinsic_info);
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                partition_function_target->SetExtrinsic(k,l,basepairing_extrinsic_info->x(k,l) );

            }
        }
        partition_function_target->PartitionFunction();

        for(int k =1;k<= current_position1;++k){
            double singlestranded_probability = 1;
            for(int l=1; l<=k-1;++l) 
                singlestranded_probability -= partition_function_target->GetPairProbability(l,k);
            
            for(int l=k+1;l<= current_position1;++l){
                singlestranded_probability -= partition_function_target->GetPairProbability(k,l);
                if(alignment_map[k]==0 && alignment_map[l]==0){
                    res(2) += partition_function_target->GetPairProbability(k,l);
                    res(2) -= (structure1->GetPair(k) == l);
                }
                else if(alignment_map[k]!=0 && alignment_map[l]==0 && structure2->GetPair(alignment_map[k],1)==0){
                    res(6) += partition_function_target->GetPairProbability(k,l);
                    res(6) -= (structure1->GetPair(k) == l);
                }
                else if(alignment_map[k]==0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[l],1)==0){
                    res(6) += partition_function_target->GetPairProbability(k,l);
                    res(6) -= (structure1->GetPair(k) == l);
                }
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)==0 && structure2->GetPair(alignment_map[l],1)==0){
                    res(5) += partition_function_target->GetPairProbability(k,l);
                    res(5) -= (structure1->GetPair(k) == l);
                }
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)!=0 && structure2->GetPair(alignment_map[l],1)!=0){
                    res(3) += partition_function_target->GetPairProbability(k,l);
                    res(3) -= (structure1->GetPair(k) == l);
               
                }
                else {
//                    cerr << "something wrong with base pair "<<k<<" "<<l<<" in sequence pair "<<temp_files[0]<<"\n";
//                    exit(EXIT_FAILURE);
                }



            }
            
            if(alignment_map[k]==0){
                res(0) += singlestranded_probability;
                res(0) -= (structure1->GetPair(k) == 0);

            }
            else if(alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)!=0){
                res(4) += singlestranded_probability;
                res(4) -= (structure1->GetPair(k) == 0);
        }
        else if (alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)==0){
            res(1) += singlestranded_probability;
            res(1) -= (structure1->GetPair(k) == 0);
        }
        else {
            cerr << "something wrong with nucleotide "<<k<<" in sequence pair "<<temp_files[0]<<"\n";
            exit(EXIT_FAILURE);}


        }
        delete basepairing_extrinsic_info;
        delete singlestranded_extrinsic_info;
        delete structure1;
        delete structure2;
        delete partition_function_target;
      


      }

    for(int p=0;p<7;++p){
        res(p) /= training_files.size();
        res(p) += alpha*m(p); 
    }

 
      return res;

}



double conserved_training_interface::conserved_testing_accuracy(const column_vector& m)
{
    vector<double> log_accuracy(training,0.0);
    double sum_log_accuracy = 0.0;
    for(int i = training; i< testing_files.size()/3; ++i) {
        vector<string> temp_files;
        temp_files.push_back(testing_files[i*3]);
        temp_files.push_back(testing_files[i*3+1]);
        temp_files.push_back(testing_files[i*3+2]);
        ifstream myfile(temp_files[0].c_str() );
//        cerr << temp_files[0] <<"\n";
        string temp;
        string sequence1;
        string sequence2;
        getline(myfile,temp);
        getline(myfile,sequence1);
        getline(myfile,temp);
        getline(myfile,sequence2);
//        sequence1 = sequence1.substr(0,sequence1.size()-2);
//        sequence2 = sequence2.substr(0,sequence2.size()-2);
//        cerr <<"sequence1 "<<sequence1<<"\n";
//        cerr <<"sequence2 "<<sequence2<<"\n";
        string target_sequence = "";
        map<int,int> alignment_map;
        int current_position1 = 0;
        int current_position2 = 0;
        for(int j = 0;j<sequence1.size();++j) {
            if(sequence2[j]!='-') ++current_position2;
            if(sequence1[j]!='-') {
                target_sequence.push_back(sequence1[j]);
                ++current_position1;
                if (sequence2[j]!='-'){alignment_map[current_position1] = current_position2;}
                else{alignment_map[current_position1] = 0;}
                    
                }
                
        }
//        cerr <<"target sequence "<<target_sequence<<"\n";
        t_matrix* basepairing_extrinsic_info = new t_matrix(current_position1, current_position1, true);
        vector<double>* singlestranded_extrinsic_info = new vector<double>(current_position1+1,0);
        RNA* structure1 = new RNA(temp_files[1].c_str(),1,true);
        RNA* structure2 = new RNA(temp_files[2].c_str(),1,true);
        RNA* partition_function_target = new RNA(target_sequence.c_str() );
        for(int k =1;k<= current_position1;++k){

            for(int l=k+1;l<= current_position1;++l){
                if(alignment_map[k]==0 && alignment_map[l]==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(2) );}
                else if(alignment_map[k]!=0 && alignment_map[l]==0 && structure2->GetPair(alignment_map[k],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(6) );}
                else if(alignment_map[k]==0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[l],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(6) );}
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)==0 && structure2->GetPair(alignment_map[l],1)==0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(5) );}
                else if(alignment_map[k]!=0 && alignment_map[l]!=0 && structure2->GetPair(alignment_map[k],1)!=0 && structure2->GetPair(alignment_map[l],1)!=0){
                    basepairing_extrinsic_info->x(k,l) = exp(m(3) );}
                else {
                    //      cerr << "something wrong with base pair "<<k<<" "<<l<<" in sequence pair "<<temp_files[0]<<"\n";
//                    exit(EXIT_FAILURE);
                    basepairing_extrinsic_info->x(k,l) = 0;

                }

                // partition_function_target->SetExtrinsic(k,l,basepairing_extrinsic_info->x(k,l) );

           }
            
            if(alignment_map[k]==0){
                singlestranded_extrinsic_info->at(k) = exp(m(0) );
            }
            else if(alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)!=0){
                singlestranded_extrinsic_info->at(k) = exp(m(4) );
        }
        else if (alignment_map[k]!=0 && structure2->GetPair(alignment_map[k],1)==0){
            singlestranded_extrinsic_info->at(k) = exp(m(1) );
        }
        else {
            //      cerr << "something wrong with nucleotide "<<k<<" in sequence pair "<<temp_files[0]<<"\n";
            exit(EXIT_FAILURE);}

            

        }
       

        
/*
         double max_extrinsic_info = 0;
        for(int k =1;k<= current_position1;++k){
            if (singlestranded_extrinsic_info->at(k) > max_extrinsic_info)
                max_extrinsic_info = singlestranded_extrinsic_info->at(k);
        }
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                if(basepairing_extrinsic_info->x(k,l) > max_extrinsic_info*max_extrinsic_info)
                    max_extrinsic_info = sqrt(basepairing_extrinsic_info->x(k,l));
            }
        }
        
        for(int k =1;k<= current_position1;++k){
            singlestranded_extrinsic_info->at(k) /= max_extrinsic_info;
        }
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                basepairing_extrinsic_info->x(k,l) /= max_extrinsic_info*max_extrinsic_info;
            }
        }
*/      




        
        partition_function_target->SetSinglestranded( singlestranded_extrinsic_info);

        
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                partition_function_target->SetExtrinsic(k,l,basepairing_extrinsic_info->x(k,l) );

            }
        }
        partition_function_target->PartitionFunction();
        double energy_partition_function = partition_function_target->GetEnsembleEnergy();

        RNA* equilibrium_constant_target = new RNA(target_sequence.c_str() );
        for(int k =1;k<= current_position1;++k){
            for(int l=k+1;l<= current_position1;++l){
                if (structure1->GetPair(k,1)!=l){
                    basepairing_extrinsic_info->x(k,l) = 0;
                }
                else {
//                    cerr<<i<<" " <<k<<" "<<l<<" "<<basepairing_extrinsic_info->x(k,l)<<" "<<target_sequence[k-1]<<" "<<target_sequence[l-1]<<"\n";
                }
                equilibrium_constant_target->SetExtrinsic(k,l,basepairing_extrinsic_info->x(k,l) );

            }
            if (structure1 ->GetPair(k,1)!=0){
                singlestranded_extrinsic_info->at(k) = 0;
            }
            else {
//                cerr<<i<<" "<<k<<" "<<singlestranded_extrinsic_info->at(k)<<"\n";
            }
        }
        
        equilibrium_constant_target->SetSinglestranded( singlestranded_extrinsic_info);
//        cerr << "testing "<<i<<"\n";
        equilibrium_constant_target->PartitionFunction();
        double energy_equilibrium_constant = equilibrium_constant_target->GetEnsembleEnergy();
        
       sum_log_accuracy += energy_partition_function/(-1*RKC*310.15);
//        cerr <<"accuracy "<<log_accuracy[i]<<"\n";
        sum_log_accuracy -= energy_equilibrium_constant/(-1*RKC*310.15);
        //       cerr <<"accuracy "<<log_accuracy[i]<<"\n";
        //     cerr << "partition "<<energy_partition_function<<"\n";
        //         cerr << "eq "<<energy_equilibrium_constant<<"\n";

        delete structure1;
        delete structure2;
        delete basepairing_extrinsic_info;
        delete singlestranded_extrinsic_info;
        delete partition_function_target;
        delete equilibrium_constant_target;
    
   }
    

//    for(int i=0;i<7;++i){
//        cerr <<"i "<<i<<" "<<0.5*alpha*m(i)*m(i)<<"\n";
//        cerr <<m(i)<<" "<<alpha<<"\n";
//        sum_log_accuracy += 0.5*alpha*m(i)*m(i);

        //  }
//    cerr <<"final "<<sum_log_accuracy<<"\n";
    return sum_log_accuracy/testing_files.size();
} 

conserved_training_interface* runner = new conserved_training_interface();

double conserved_accuracy_wrapper(const column_vector& m){
    return runner->conserved_accuracy(m);
}

const column_vector conserved_accuracy_derivative_wrapper(const column_vector& m){
    return runner->conserved_accuracy_derivative(m);
}

int main( int argc, char* argv[] ){


    bool parseable = runner->parse(argc,argv);
    column_vector starting_point(7);  
    if (parseable == true){
        starting_point(0) = log(runner->parameters[1]);
        starting_point(1) = log(runner->parameters[2]);
        starting_point(2) = log(runner->parameters[3]);
        
        starting_point(3) = log(runner->parameters[4]);
        starting_point(4) = log(runner->parameters[5]);
        starting_point(5) = log(runner->parameters[6]);
        starting_point(6) = log(runner->parameters[7]);
    
        dlib::find_min(dlib::lbfgs_search_strategy(30), 
                       dlib::objective_delta_stop_strategy(1e-2).be_verbose(),
                       conserved_accuracy_wrapper,
                       conserved_accuracy_derivative_wrapper,
                       starting_point,
                       -1);
        
        cout << "Parameter1: "<<exp(starting_point(0) )<<"\n";
        cout << "Parameter2: "<<exp(starting_point(1) )<<"\n";
        cout << "Parameter3: "<<exp(starting_point(2) )<<"\n";
        cout << "Parameter4: "<<exp(starting_point(3) )<<"\n";
        cout << "Parameter5: "<<exp(starting_point(4) )<<"\n";
        cout << "Parameter6: "<<exp(starting_point(5) )<<"\n";
        cout << "Parameter7: "<<exp(starting_point(6) )<<"\n";
/*
        starting_point(0) = log(0.986149);
        starting_point(1) = log(3.77393);
        starting_point(2) = log(1.00976);
        starting_point(3) = log(1.97015);
        starting_point(4) = log(0.257632);
        starting_point(5) = log(0.516168);
        starting_point(6) = log(0.994545);

        starting_point(0) = log(1);
        starting_point(1) = log(1);
        starting_point(2) = log(1);
        starting_point(3) = log(1);
        starting_point(4) = log(1);
        starting_point(5) = log(1);
        starting_point(6) = log(1);
*/

        double testing_accuracy = runner->conserved_testing_accuracy(starting_point);
        cout << "Testing Accuracy: "<<testing_accuracy<<"\n";
        

    }
    
    delete runner;
    return 0;

}
