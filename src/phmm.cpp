#include "phmm.h"
using std::string;

t_phmm_aln* create_phmm_aln(structure *ct1, structure *ct2){

	//put seq1 nucleotides into a vector
    	vector<char>* seq1_nucs = new vector<char>();
    	for(int i = 1; i <= ct1->GetSequenceLength(); i++){
		seq1_nucs->push_back(ct1->nucs[i]);
    	} // i loop.
	
	// Sequence2 nucleotides.
	vector<char>* seq2_nucs = new vector<char>();
	for(int j = 1; j <= ct2->GetSequenceLength(); j++){
		seq2_nucs->push_back(ct2->nucs[j]);
	} // i loop.
	
	// Allocate the structures.
	t_structure* str1 = new t_structure("seq1", seq1_nucs);
	t_structure* str2 = new t_structure("seq2", seq2_nucs);
	
	// Allocate the structures.
	t_phmm_aln* phmm_aln = new t_phmm_aln(str1, str2);
	
	return phmm_aln;
}

void write_probability_array(t_pp_result* pp_result, const char* outputfile,int l1,int l2){
/*
	for(int i = 0; i < l1;i++){
		std::cout <<"\n";
		for(int j = 0; j< l2;j++)	
			std::cout << pp_result->aln_probs[i][j] <<"\n";
	}
*/
	ofstream out(outputfile);
	for(int j = 0; j < l2;j++)
		out << "\t" << j+1;
	for(int i = 0; i < l1;i++){
		out <<"\n" << i+1;
		for(int j = 0; j< l2;j++)	
			out <<"\t" << pp_result->aln_probs[i][j];
	}
	out.close();

}

void write_ML_alignment(t_ML_result* ML_result, const char* outputfile, int l1, int l2, const char* seq1, const char* seq2){
	ofstream out(outputfile);
	out << "Maximum likelihood alignment between " << seq1 << " and " << seq2 <<"\n\n" ;
	for(int i = 0;i<ML_result->seq1_aln_line->size();i++)
		out << ML_result->seq1_aln_line->at(i);
	out << "\n";
	for(int i = 0;i<ML_result->seq2_aln_line->size();i++)
		out << ML_result->seq2_aln_line->at(i);
	out.close();

}
