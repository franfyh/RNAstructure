
/* 	RNA Secondary Structure Prediction, Using the Algorithm of Zuker
	C++ version by David Mathews, copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006

	Addition programming by Zhi John Lu, 2005, 2006, to revise internal loop discovery to an O(N^3) algorithm.

	Programmed for Isis Pharmaceuticals, the Turner Lab, and the Mathews Lab
	Department of Biochemistry and Biophysics, University of Rochester Medical Center

	Revised on the basis of current research:
	Mathews, Sabina, Zuker, & Turner.  1999.  JMB: 288:911-940.
	Mathews, Childs, Disney, Schroeder, Zuker, & Turner.  2004. PNAS: 101:7287-7292. 
 */

#include "structure.h"
#include "algorithm.h"
#ifdef _WINDOWS_GUI
#include "../RNAstructure_windows_interface/platform.h"
#else 
#include "platform.h" 
#endif //_WINDOWS
#include <cstring>
#include <cmath>
#include <cstdlib>
using namespace std;




//#define maxfil 250    //maximum length of file names

#define maxtloop 200 //maximum tetraloops allowed (info read from tloop)
//#define ctheaderlength 125 //maximum length of string containing info on sequence

#define numlen 8  //maximum digits in a number
#define maxnopair 600 //maximum number of bases that can be forced single
#undef debugmode //a flag to turn off debugging features
//#define debugmode //a flag to turn on debugging features
//flags for debugging
#undef timer
//#define timer //flag to indicate the code execution should be timed

//********************************functions:

/*	Function efn2

	Calculates the free energy of each structure in a structure called structure.

	Structures cannot have pseudoknots

	Structnum indicates which structure to calculate the free energy
	the default, 0, indicates "all" structures

	simplemb = true means that the energy function for mb loops is linear,
	so identical to the dynamic programming algorithms.
	The default is false, so logarithmic in the number of unpaired nucleotides.

	This neew version of efn2 is now stacking aware -- if ct->stacking==true, it
	will use the stacking information stored in ct, rather than inferring stacking
	by free energy minimization.

 */
#ifndef INSTRUMENTED//Execute if pre-compiler flag INSTRUMENTED is not defined. Instrumented is only defined when compiling NAPSS.
//If INSTRUMENTED is defined, functions 'dynamic' and 'fill' will be overloaded

void efn2(datatable *data,structure *ct, int structnum, bool simplemb, const char *outputfilename)
{
	int i,j,k,open,null,stz,count,sum,ip,jp;
	stackstruct stack;
	forceclass fce(ct->GetSequenceLength());
	int start, stop;
	ofstream out;
	char temp[maxfil];
	int tempsum;
	int energy;

	/*	stack = a place to keep track of where efn2 is located in a structure
		inter = indicates whether there is an intermolecular interaction involved in
		a multi-branch loop
	 */

	stack.sp = 0;  //set stack counter

	if (ct->intermolecular) {//this indicates an intermolecular folding
		for (i=0;i<3;i++) {
			forceinterefn(ct->inter[i],ct,&fce);
		}
	}
	if (structnum!=0) {
		start = structnum;
		stop = structnum;

	}
	else {
		start = 1;
		stop = ct->GetNumberofStructures();
	}

	if (outputfilename!=NULL) { //open a file to contain the thermodynamic details
		out.open(outputfilename);
	}

	for (count=start;count<=stop;count++){//one structure at a time

		if (outputfilename!=NULL) { //open an output file for debugging


			out << "Thermodynamic details for structure # "<<count<<"\n";
		}
//<<<<<<< algorithm.cpp
		
//		
//=======

		energy = ergexterior(count, ct, data);
		ct->SetEnergy(count,energy);
//>>>>>>> 1.40

		if (outputfilename!=NULL) {
			sprintf(temp,DIGITS,(float (energy)/conversionfactor));
			//gcvt((float (ct->energy[count]))/conversionfactor,6,temp);
			out << "\tExterior loop = "<<temp<<"\n";	
		}

		i=1;
		while (i<ct->GetSequenceLength()) { 	
			if (ct->GetPair(i,count)!=0) {
				push(&stack,i,ct->GetPair(i,count),1,0);
				i = ct->GetPair(i,count);
			}
			i++;
		}

subroutine://loop starts here (I chose the goto statement to make the code more readable)
		pull(&stack,&i,&j,&open,&null,&stz);//take a substructure off the stack

		while (stz!=1){
			while (ct->GetPair(i,count)==j) { //are i and j paired?
				if (outputfilename!=NULL) {
					tempsum=0;
				}
//<<<<<<< algorithm.cpp
			
//				while (ct->GetPair(i+1,count)==j-1) {//are i,j and i+1,j-1 stacked?
//					energy += erg1(i,j,i+1,j-1,ct,data);
//=======

				while (ct->GetPair(i+1,count)==j-1) {//are i,j and i+1,j-1 stacked?
					ct->SetEnergy(count,ct->GetEnergy(count) + erg1(i,j,i+1,j-1,ct,data));
//>>>>>>> 1.40

					if (outputfilename!=NULL) {
						tempsum = tempsum + erg1(i,j,i+1,j-1,ct,data); 
						sprintf(temp,DIGITS,(float (erg1(i,j,i+1,j-1,ct,data)))/conversionfactor);
						//gcvt((float (erg1(i,j,i+1,j-1,ct,data)))/conversionfactor,6,temp);
						out << "\t\tStack = "<<temp<<"  for stack of "<<
							i<<"-"<<j<<"\n";
					}
					i++;
					j--;
				}

				if (outputfilename!=NULL) {
					sprintf(temp,DIGITS,(float (tempsum))/conversionfactor);
					//gcvt((float (tempsum))/conversionfactor,6,temp);
					out << "\tHelix total = "<<temp<<"\n";	
				}

				sum = 0;
				k = i + 1;
				// now efn2 is past the paired region, so define
				//		 the intervening non-paired segment

				while (k<j) {
					if (ct->GetPair(k,count)>k)	{
						sum++;
						ip = k;
						k = ct->GetPair(k,count) + 1;
						jp = k-1;
					}
					else if (ct->GetPair(k,count)==0) k++;
				}

				if (sum==0) {//hairpin loop
					
					ct->SetEnergy(count,ct->GetEnergy(count) + erg3(i,j,ct,data,fce.f(i,j)));
					if (outputfilename!=NULL) {
						sprintf(temp,DIGITS,(float (erg3(i,j,ct,data,fce.f(i,j))))/conversionfactor);
						//gcvt((float (erg3(i,j,ct,data,fce.f(i,j))))/conversionfactor,6,temp);
						out << "\tHairpin = "<<temp<<"  for closure of "<<
							i<<"-"<<j<<"\n";	
					}
					goto subroutine;
				}
				else if (sum==1) {//bulge/internal loop

					
					ct->SetEnergy(count,ct->GetEnergy(count) + erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j)));
					if (outputfilename!=NULL) {

						sprintf(temp,DIGITS,(float (erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))))/conversionfactor);
						//gcvt((float (erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))))/conversionfactor,6,temp);
						out << "\tInternal/bulge = "<<temp<<"  for closure of "<<
							i<<"-"<<j<<"\t" << "Size = " << (ip-i+j-jp-2);
						if (	(ip-i+j-jp-2) == 1) {
							//This is a bulge loop with a single bulge, so report the stack information
							sprintf(temp,DIGITS,((float)data->stack[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]])/conversionfactor);
							//gcvt(((float)data->stack[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]])/conversionfactor,6,temp);
							out << " BP stack across bulge = " << temp;
						}
						out << "\n";	
					}
					i = ip;
					j = jp;
				}
				else {//multi-branch loop
//<<<<<<< algorithm.cpp
				
//					energy += ergmulti(count, i, ct, data, simplemb);
				
				
//=======

					ct->SetEnergy(count,ct->GetEnergy(count)+ 
						ergmulti(count, i, ct, data, simplemb));


//>>>>>>> 1.40
					if (outputfilename!=NULL) {
						sprintf(temp,DIGITS,(float (ergmulti(count, i, ct, data, simplemb)))/conversionfactor);
						//gcvt((float (ergmulti(count, i, ct, data, simplemb)))/conversionfactor,6,temp);
						out << "\tMultiloop = "<<temp<<"  for closure of "<<
							i<<"-"<<j<<"\n";	
					}

					//put the exiting helices on the stack:
					sum++;//total helixes = sum + 1
					i++;
					for (k=1;k<sum;k++) {
						while (ct->GetPair(i,count)==0) i++;
						push (&stack,i,ct->GetPair(i,count),1,0);
						i = ct->GetPair(i,count)+1;
					}

					goto subroutine;

				}
			}    
		}
//<<<<<<< algorithm.cpp
	
		//Strore the energy in the underlying structure class:
//		ct->SetEnergy(count,energy);

		//If a report file is being written, write the total energy to it:
//=======

//>>>>>>> 1.40
		if (outputfilename!=NULL) {
			sprintf(temp,DIGITS,(float (ct->GetEnergy(count)))/conversionfactor);
			out << "\n\n\tTotal energy = "<<temp<<"\n\n\n\n";
		}
	}

	if (outputfilename!=NULL) {
		out.close();
	}   
	return;
}

//This function adds the pair between i and j to structure ct.  
//This pair is added to the last structure structured in ct.
//This function is careful to register the pair for the real nucleotide indexes.  This is important because the pair could have been discovered in an exterior fragment.
void registerbasepair(structure *ct, short i, short j) {
	
	//i and j are paired -- put that information in structure ct
	if (j<=ct->GetSequenceLength()) {
		ct->SetPair(i,j,ct->GetNumberofStructures());
	}
	else if (i>ct->GetSequenceLength()) {

		i = i - ct->GetSequenceLength();
		j = j - ct->GetSequenceLength();
		ct->SetPair(i,j,ct->GetNumberofStructures());
	}
	else {
		j = j-ct->GetSequenceLength();
		ct->SetPair(i,j,ct->GetNumberofStructures());
	}
	return;
}


//This function, given a pair between ii and ji, will determine base pairs
//in the lowest free energy structure that contain that pair.

//Pairs are stored in structure ct at the last structure number recorded in ct.
//This returns an error code, where zero is no error and non-zero indicates a traceback error.

//If quickstructure=true, then the function knows that the fragment to traceback is from 1 to N, and the W5 array should be traced.
#ifdef DYNALIGN_II
int trace(structure *ct, datatable *data, int ii, int ji,
          arrayclass *v, arrayclass *w, arrayclass *wmb, arrayclass *w2, arrayclass *wmb2, 
          bool *lfce, forceclass *fce, integersize *w3, integersize *w5,bool *mod,arrayclass *we,integersize energy,short open,short pair, bool quickstructure)
#else
int trace(structure *ct, datatable *data, int ii, int ji,
		arrayclass *v, arrayclass *w, arrayclass *wmb, arrayclass *w2, arrayclass *wmb2, 
		bool *lfce, forceclass *fce, integersize *w3, integersize *w5,bool *mod,bool quickstructure)
#endif
{
	stackclass *stack;
	register short number;
	short i,j,k,a,b,d;
#ifndef DYNALIGN_II
        short open, pair;
	integersize energy;
#else
#endif

	bool found;
	int tracebackerror=0;
	register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},{0,1,0,1,0,0},{0,0,0,0,0,0}};
//<<<<<<< algorithm.cpp
//	                        {0,1,0,1,0,0},{0,0,0,0,0,0}};
	
//	number = ct->GetSequenceLength();
//=======
		

	number = ct->GetSequenceLength();
//>>>>>>> 1.40
	stack = new stackclass();

	
#ifdef DYNALIGN_II
	stack->push(ii,ji,open,energy,pair);
#else
	//initialize the stack:
	if (quickstructure) {
		//Tracebeck from W5:
		stack->push(1,ct->GetSequenceLength(),1,w5[ct->GetSequenceLength()],0);

		ct->AddStructure();
		//ct->numofstructures=1;

	}

	else stack->push(ii,ji,0,v->f(ii,ji),1); //normal traceback from a pair
#endif
	
	//zero all basepairs:
	//Move this to calling function...
	//if (ji<=(number)) {
	//	for (k=ii;k<=ji;k++) ct->basepr[ct->numofstructures][k] = 0;
	//}
	//else {
	//	for (k=1;k<=(ji-(number));k++) ct->basepr[ct->numofstructures][k]=0;
	//	for (k=ii;k<=(number);k++) ct->basepr[ct->numofstructures][k]=0;
	//}


	//Add forced pairs to the underlying structure ct
	for (k=0;k<ct->GetNumberofPairs();k++) {
		ct->SetPair(ct->GetPair5(k),ct->GetPair3(k),ct->GetNumberofStructures());
	}


	while (stack->pull(&i,&j,&open,&energy,&pair)) {

		found = true;
		//keep taking pairs off the stack until the stack is empty:
		if (pair==0) {
			found = false;

			//a bifurcated segment has been taken off the stack:
			if (open) {
				//this is an exterior loop
				//so either i==1 or j==n
#ifndef DYNALIGN_II
				//remove nucleotides:
				if (j==number) {
					while(energy==w3[i+1]+ ct->SHAPEss_give_value(i)&&i<number) {
						energy-=ct->SHAPEss_give_value(i);
						i++;				
						if (i==number) break;
					}
				}
				if (i==1) {
					//i == 1
					while(energy==w5[j-1]+ ct->SHAPEss_give_value(j)&&j!=1) {
						energy-=ct->SHAPEss_give_value(j);
						j--;
					}
				}

				if (i!=number&&j!=1) {

					if (energy==v->f(i+1,j) + penalty(i+1,j,ct,data) +
							erg4(j,i+1,i,2,ct,data,lfce[i])) {
						i++;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[j]]) {
						if (energy==v->f(i+2,j-1) + penalty(i+1,j,ct,data) +
								erg4(j,i+1,i,2,ct,data,lfce[i])+erg1(i+1,j,i+2,j-1,ct,data)) {

							stack->push(i+2,j-1,0,v->f(i+2,j-1),1);
							registerbasepair(ct,i+1,j);
							found = true;
						}
					}
					if (!found&&(energy==v->f(i,j-1) + penalty(i,j-1,ct,data) +
								erg4(j-1,i,j,1,ct,data,lfce[j]))) {
						j--;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if ((!found)&&(mod[i]||mod[j-1])&&inc[ct->numseq[i]][ct->numseq[j-1]]) {
						if ((energy==v->f(i+1,j-2) + penalty(i,j-1,ct,data) +
									erg4(j-1,i,j,1,ct,data,lfce[j])
									+erg1(i,j-1,i+1,j-2,ct,data))) {

							registerbasepair(ct,i,j-1);
							stack->push(i+1,j-2,0,v->f(i+1,j-2),1);
							found = true;
						}
					}

					if (!found&&(energy==v->f(i+1,j-1) + penalty(i+1,j-1,ct,data) +
								data->tstack[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]]
								[ct->numseq[i]]+checknp(lfce[i],lfce[j])+ ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(j))) {
						j--;
						i++;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if (!found&&(mod[i+1]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-1]]) {

						if (energy==v->f(i+2,j-2) + penalty(i+1,j-1,ct,data) +
								data->tstack[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]]
								[ct->numseq[i]]+checknp(lfce[i],lfce[j])+erg1(i+1,j-1,i+2,j-2,ct,data)+ ct->SHAPEss_give_value(i)+ 
								ct->SHAPEss_give_value(j)) {

							registerbasepair(ct,i+1,j-1);
							stack->push(i+2,j-2,0,v->f(i+2,j-2),1);
							found = true;
						}
					}

					if (!found&&(energy==v->f(i,j)+penalty(i,j,ct,data))) {
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if (!found&&(mod[i]||mod[j])&&inc[ct->numseq[i]][ct->numseq[j]]) {
						if (energy==v->f(i+1,j-1)+penalty(i+1,j-1,ct,data)+erg1(i,j,i+1,j-1,ct,data)) {
							registerbasepair(ct,i,j);
							stack->push(i+1,j-1,0,v->f(i+1,j-1),1);
							found = true;
						}
					}
				}


#else

									//start by removing nucleotides:
				while (energy==we->f(i+1,j)) {
				  i++;
				  energy = we->f(i,j);
				  if (i==j) break;
				  
				}
				while (energy==we->f(i,j-1)) {
				  j--;
				  energy = we->f(i,j);
				  if (i==j) break;
				  
				}

					//now check for a stem:
					
				if (energy==v->f(i+1,j)+erg4(j,i+1,i,2,ct,data,lfce[i])+
				    penalty(i+1,j,ct,data)) {
				  
				  i++;
				  stack->push(i,j,0,v->f(i,j),1);
				  found = true;
				}
				
				else if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[j]]) {

				  if (energy==v->f(i+2,j-1)+erg4(j,i+1,i,2,ct,data,lfce[i])+
				      penalty(i+1,j,ct,data)+erg1(i+1,j,i+2,j-1,ct,data)) {

				    registerbasepair(ct,i+1,j);
				    stack->push(i+2,j-1,0,v->f(i+2,j-1),1);
				    found = true;
				  }
				}

				if (!found&&(energy==v->f(i,j-1)+erg4(j-1,i,j,1,ct,data,lfce[j])+
					     penalty(i,j-1,ct,data))) {
				  
				  j--;
				  stack->push(i,j,0,v->f(i,j),1);
				  found = true;
				}
				    else if (!found&&(mod[i]||mod[j-1])&&inc[ct->numseq[i]][ct->numseq[j-1]]) {
				      
				      if ((energy==v->f(i+1,j-2)+erg4(j-1,i,j,1,ct,data,lfce[j])+
					   penalty(i,j-1,ct,data)+erg1(i,j-1,i+1,j-2,ct,data))) {
					
					registerbasepair(ct,i,j-1);
					stack->push(i+1,j-2,0,v->f(i+1,j-2),1);
					found = true;
				      }
				    }

				    if (!found&&energy==v->f(i+1,j-1)+
					data->tstack[ct->numseq[j-1]][ct->numseq[i+1]]//tstkm?
					[ct->numseq[j]][ct->numseq[i]] + checknp(lfce[i],lfce[j])+
					penalty(i+1,j-1,ct,data)) {
				      
				      i++;
				      j--;
				      stack->push(i,j,0,v->f(i,j),1);
				      found = true;
				    }
				    else if (!found&&(mod[i+1]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-1]]) {
				      if (energy==v->f(i+2,j-2)+
					  data->tstack[ct->numseq[j-1]][ct->numseq[i+1]]
					  [ct->numseq[j]][ct->numseq[i]] + checknp(lfce[i],lfce[j])+
					  penalty(i+1,j-1,ct,data)
					  +erg1(i+1,j-1,i+2,j-2,ct,data)) {

					registerbasepair(ct,i+1,j-1);
					stack->push(i+2,j-2,0,v->f(i+2,j-2),1);
					found = true;
				      }
				    }

				    if (!found&&energy==v->f(i,j)+penalty(i,j,ct,data)) {
				      stack->push(i,j,0,v->f(i,j),1);
				      found = true;
				    }
				    else if (!found&&(mod[i]||mod[j])&&inc[ct->numseq[i]][ct->numseq[j]]) {
				      
				      if (energy==v->f(i+1,j-1)+penalty(i,j,ct,data)+erg1(i,j,i+1,j-1,ct,data)) {
					
					registerbasepair(ct,i,j);
					stack->push(i+1,j-1,0,v->f(i+1,j-1),1);
					found = true;
				      }
				    }

#endif//yinghan1


			}

			else {
				//(open==0), a multibranch loop

				if (ct->intermolecular) {
					//intermolcular folding 

					//start by removing nucleotides:
					while (energy==w2->f(i+1,j)) {
						i++;
						energy = w2->f(i,j);
						if (i==j) break;
					}
					while (energy==w2->f(i,j-1)) {
						j--;
						energy = w2->f(i,j);
						if (i==j) break;
					}

					//following commented out by DHM (5/17/08) -- this should be handled by WMB2
					//if (fce->f(i,i)&INTER) {

					//	if (energy==w2->f(i+1,j)-INFINITE_ENERGY+data->init) {
					//		i++;
					//		energy = w2->f(i,j);

					//	}
					//}
					//if (fce->f(j,j)&INTER) {

					//	if (energy==w2->f(i,j-1)-INFINITE_ENERGY+data->init) {
					//		j--;
					//		energy = w2->f(i,j);

					//	}
					//}

					while (energy==w2->f(i+1,j)) {
						i++;
						energy = w2->f(i,j);
						if (i==j) break;
					}
					while (energy==w2->f(i,j-1)) {
						j--;
						energy = w2->f(i,j);
						if (i==j) break;
					}

					//now check for a stem:					
					if (energy==v->f(i+1,j)+erg4(j,i+1,i,2,ct,data,lfce[i])+
							penalty(i+1,j,ct,data)) {

						i++;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[j]]) {

						if (energy==v->f(i+2,j-1)+erg4(j,i+1,i,2,ct,data,lfce[i])+
								penalty(i+1,j,ct,data)+erg1(i+1,j,i+2,j-1,ct,data)) {

							registerbasepair(ct,i+1,j);
							stack->push(i+2,j-1,0,v->f(i+2,j-1),1);
							found = true;
						}
					}
					if (!found&&(energy==v->f(i,j-1)+erg4(j-1,i,j,1,ct,data,lfce[j])+
								penalty(i,j-1,ct,data))) {

						j--;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}
					else if (!found&&(mod[i]||mod[j-1])&&inc[ct->numseq[i]][ct->numseq[j-1]]) { 

						if (!found&&(energy==v->f(i+1,j-2)+erg4(j-1,i,j,1,ct,data,lfce[j])+
									penalty(i,j-1,ct,data)+erg1(i,j-1,i+1,j-2,ct,data))) {

							registerbasepair(ct,i,j);
							stack->push(i+1,j-2,0,v->f(i+1,j-2),1);
							found = true;
						}
					}

					if (!found&&(energy==v->f(i+1,j-1)+
								data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
								[ct->numseq[j]][ct->numseq[i]] + ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(j) + checknp(lfce[i],lfce[j])+
								penalty(i+1,j-1,ct,data))) {

						i++;
						j--;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if (!found&&(mod[i+1]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-1]]) {
						if (!found&&(energy==v->f(i+2,j-2)+
									data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
									[ct->numseq[j]][ct->numseq[i]] + ct->SHAPEss_give_value(i) + ct->SHAPEss_give_value(j) + checknp(lfce[i],lfce[j])+
									penalty(i+1,j-1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data))) {

							registerbasepair(ct,i+1,j-1);
							stack->push(i+2,j-2,0,v->f(i+2,j-2),1);
							found = true;
						}
					}

					if (!found&&(energy==v->f(i,j)+penalty(i,j,ct,data))) {
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if (!found&&(mod[i]||mod[j])&&inc[ct->numseq[i]][ct->numseq[j]]) {
						if (!found&&(energy==v->f(i+1,j-1)+penalty(i,j,ct,data)
									+erg1(i,j,i+1,j-1,ct,data))) {

							registerbasepair(ct,i,j);
							stack->push(i+1,j-1,0,v->f(i+1,j-1),1);
							found = true;
						}
					}
				}

				if (!found) {

					//start by removing nucleotides:
					while (energy==w->f(i+1,j)+data->eparam[6]+ ct->SHAPEss_give_value(i)) {
						i++;
						energy = w->f(i,j);
						if (i==j) break;

					}
					while (energy==w->f(i,j-1)+data->eparam[6]+ ct->SHAPEss_give_value(j)&&j!=ct->GetSequenceLength()+1) {
						j--;
						energy = w->f(i,j);
						if (i==j) break;

					}

					//now check for a stem:

					if (energy==v->f(i+1,j)+erg4(j,i+1,i,2,ct,data,lfce[i])+
							penalty(i+1,j,ct,data)+data->eparam[10]+data->eparam[6]) {

						i++;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}

					else if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[j]]) {

						if (energy==v->f(i+2,j-1)+erg4(j,i+1,i,2,ct,data,lfce[i])+
								penalty(i+1,j,ct,data)+data->eparam[10]+data->eparam[6]
								+erg1(i+1,j,i+2,j-1,ct,data)) {

							registerbasepair(ct,i+1,j);
							stack->push(i+2,j-1,0,v->f(i+2,j-1),1);
							found = true;
						}
					}

					if (!found&&(energy==v->f(i,j-1)+erg4(j-1,i,j,1,ct,data,lfce[j])+
								penalty(i,j-1,ct,data)+data->eparam[10]+data->eparam[6])) {

						j--;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}
					else if (!found&&(mod[i]||mod[j-1])&&inc[ct->numseq[i]][ct->numseq[j-1]]) {

						if ((energy==v->f(i+1,j-2)+erg4(j-1,i,j,1,ct,data,lfce[j])+
									penalty(i,j-1,ct,data)+data->eparam[10]+data->eparam[6]
									+erg1(i,j-1,i+1,j-2,ct,data))) {

							registerbasepair(ct,i,j-1);
							stack->push(i+1,j-2,0,v->f(i+1,j-2),1);
							found = true;
						}
					}

					if (!found&&energy==v->f(i+1,j-1)+
							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
							[ct->numseq[j]][ct->numseq[i]] + checknp(lfce[i],lfce[j])+
							penalty(i+1,j-1,ct,data)+data->eparam[10]+2*data->eparam[6]
							+ ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(j)) {

						i++;
						j--;
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}
					else if (!found&&(mod[i+1]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-1]]) {
						if (energy==v->f(i+2,j-2)+
								data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
								[ct->numseq[j]][ct->numseq[i]] + checknp(lfce[i],lfce[j])+
								penalty(i+1,j-1,ct,data)+data->eparam[10]+2*data->eparam[6]
								+erg1(i+1,j-1,i+2,j-2,ct,data)+ ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(j)) {

							registerbasepair(ct,i+1,j-1);
							stack->push(i+2,j-2,0,v->f(i+2,j-2),1);
							found = true;
						}
					}

					if (!found&&energy==v->f(i,j)+penalty(i,j,ct,data)+data->eparam[10]) {
						stack->push(i,j,0,v->f(i,j),1);
						found = true;
					}
					else if (!found&&(mod[i]||mod[j])&&inc[ct->numseq[i]][ct->numseq[j]]) {

						if (energy==v->f(i+1,j-1)+penalty(i,j,ct,data)+data->eparam[10]
								+erg1(i,j,i+1,j-1,ct,data)) {

							registerbasepair(ct,i,j);
							stack->push(i+1,j-1,0,v->f(i+1,j-1),1);
							found = true;
						}
					}
				}
			}
		}

		else if (pair==2) {
			found = false;
			while (energy==wmb->f(i+1,j)+data->eparam[6]+ ct->SHAPEss_give_value(i)) {
				i++;
				energy= wmb->f(i,j);
				if (i==j) break;
			}
			while (energy==wmb->f(i,j-1)+data->eparam[6]+ ct->SHAPEss_give_value(j)) {
				j--;
				energy = wmb->f(i,j);
				if (i==j) break;
			}

			if (ct->intermolecular) {
				while (energy==w2->f(i+1,j)) {
					i++;
					energy = w2->f(i,j);
					if (i==j) break;

				}
				while (energy==w2->f(i,j-1)) {
					j--;
					energy = w2->f(i,j);
					if (i==j) break;
				}
				//Commented out by DHM (5/17/08) - Should be handled by wmb2
				//if (fce->f(i,i)&INTER) {

				//	if (energy==w2->f(i+1,j)-INFINITE_ENERGY+data->init) {
				//		i++;
				//		energy = w2->f(i,j);

				//	}
				//}
				//if (fce->f(j,j)&INTER) {

				//	if (energy==w2->f(i,j-1)-INFINITE_ENERGY+data->init) {
				//		j--;
				//		energy = w2->f(i,j);

				//	}
				//}

				while (energy==w2->f(i+1,j)) {
					i++;
					energy = w2->f(i,j);
					if (i==j) break;

				}
				while (energy==w2->f(i,j-1)) {
					j--;
					energy = w2->f(i,j);
					if (i==j) break;
				}
			}
		}

		if(pair==2||!found) {
			//the structure must bifurcate

			//found = false;

			if (open) {

#ifndef DYNALIGN_II
				//exterior loop:
				if (i==1) {
					k=1;
					while(k<=j-minloop&&!found) {
						if (energy==w5[k]+v->f(k+1,j)+penalty(k+1,j,ct,data)) {
							stack->push(1,k,1,w5[k],0);
							stack->push(k+1,j,0,v->f(k+1,j),1);
							found = true;
						}
						else if (mod[k+1]||mod[j]&&inc[ct->numseq[k+1]][ct->numseq[j]]) {

							if (energy==w5[k]+v->f(k+2,j-1)+penalty(k+1,j,ct,data)
									+erg1(k+1,j,k+2,j-1,ct,data)) {

								registerbasepair(ct,k+1,j);
								stack->push(1,k,1,w5[k],0);
								stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
								found = true;
							}
						}
						if (!found&&(energy==w5[k]+v->f(k+2,j)+penalty(k+2,j,ct,data)+
									erg4(j,k+2,k+1,2,ct,data,lfce[k+1]))) {
							stack->push(1,k,1,w5[k],0);
							stack->push(k+2,j,0,v->f(k+2,j),1);
							found = true;
						}
						else if (!found&&(mod[k+2]||mod[j])&&inc[ct->numseq[k+2]][ct->numseq[j]]) {


							if (energy==w5[k]+v->f(k+3,j-1)+penalty(k+2,j,ct,data)+
									erg4(j,k+2,k+1,2,ct,data,lfce[k+1])
									+erg1(k+2,j,k+3,j-2,ct,data)) {

								registerbasepair(ct,k+2,j);
								stack->push(1,k,1,w5[k],0);
								stack->push(k+3,j-1,0,v->f(k+3,j-1),1);
								found = true;
							}
						}

						if (!found&&(energy==w5[k]+v->f(k+1,j-1)+penalty(k+1,j-1,ct,data)+
									erg4(j-1,k+1,j,1,ct,data,lfce[j]))) {

							stack->push(1,k,1,w5[k],0);
							stack->push(k+1,j-1,0,v->f(k+1,j-1),1);
							found = true;
						}
						else if (!found&&(mod[k+1]||mod[j-1])&&inc[ct->numseq[k+1]][ct->numseq[j-1]]) {

							if (!found&&(energy==w5[k]+v->f(k+2,j-2)+penalty(k+1,j-1,ct,data)+
										erg4(j-1,k+1,j,1,ct,data,lfce[j])
										+erg1(k+1,j-1,k+2,j-2,ct,data))) {

								registerbasepair(ct,k+1,j-1);
								stack->push(1,k,1,w5[k],0);
								stack->push(k+2,j-2,0,v->f(k+2,j-2),1);
								found = true;

							}
						}

						if (!found&&(energy==w5[k]+v->f(k+2,j-1)+penalty(k+2,j-1,ct,data)+
									data->tstack[ct->numseq[j-1]][ct->numseq[k+2]]
									[ct->numseq[j]][ct->numseq[k+1]]+checknp(lfce[j],lfce[k+1])
									+ ct->SHAPEss_give_value(j)+ ct->SHAPEss_give_value(k+1))) {

							stack->push(1,k,1,w5[k],0);
							stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
							found = true;
						}

						else if (!found&&(mod[k+2]||mod[j-1])&&inc[ct->numseq[k+2]][ct->numseq[j-1]]) {

							if ((energy==w5[k]+v->f(k+3,j-2)+penalty(k+2,j-1,ct,data)+
										data->tstack[ct->numseq[j-1]][ct->numseq[k+2]]
										[ct->numseq[j]][ct->numseq[k+1]]+checknp(lfce[j],lfce[k+1])
										+erg1(k+2,j-1,k+3,j-2,ct,data)
										+ ct->SHAPEss_give_value(j)+ ct->SHAPEss_give_value(k+1))) {

								registerbasepair(ct,k+2,j-1);
								stack->push(1,k,1,w5[k],0);
								stack->push(k+3,j-2,0,v->f(k+3,j-2),1);
								found = true;
							}
						}

						if (!found) {
							//check for coaxial stacking:
							a = k+minloop+1;
							while (a<j-minloop&&!found) {
								if (energy==w5[k-1]+v->f(k,a)+v->f(a+1,j)+
										ergcoaxflushbases(k,a,a+1,j,ct,data)+
										penalty(k,a,ct,data)+penalty(a+1,j,ct,data)) {

									if (k>1) stack->push(1,k-1,1,w5[k-1],0);
									stack->push(k,a,0,v->f(k,a),1);
									stack->push(a+1,j,0,v->f(a+1,j),1);
									found = true;
								}

								else if (mod[k]||mod[a]||mod[a+1]||mod[j]) {
									if ((mod[k]||mod[a])&&(mod[a+1]||mod[j])&&inc[ct->numseq[k]][ct->numseq[a]]
											&&inc[ct->numseq[a+1]][ct->numseq[j]]) {


										if (energy==w5[k-1]+v->f(k+1,a-1)+v->f(a+2,j-1)+
												ergcoaxflushbases(k,a,a+1,j,ct,data)+
												penalty(k,a,ct,data)+penalty(a+1,j,ct,data)
												+erg1(k,a,k+1,a-1,ct,data)
												+erg1(a+1,j,a+2,j-1,ct,data)) {

											registerbasepair(ct,k,a);
											registerbasepair(ct,a+1,j);
											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k+1,a-1,0,v->f(k+1,a-1),1);
											stack->push(a+2,j-1,0,v->f(a+2,j-1),1);
											found = true;
										}
									}
									if (!found&&(mod[k]||mod[a])&&inc[ct->numseq[k]][ct->numseq[a]]) {

										if (energy==w5[k-1]+v->f(k+1,a-1)+v->f(a+1,j)+
												ergcoaxflushbases(k,a,a+1,j,ct,data)+
												penalty(k,a,ct,data)+penalty(a+1,j,ct,data)
												+erg1(k,a,k+1,a-1,ct,data)) {

											registerbasepair(ct,k,a);

											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k+1,a-1,0,v->f(k+1,a-1),1);
											stack->push(a+1,j,0,v->f(a+1,j),1);
											found = true;
										}
									}

									if (!found&&(mod[a+1]||mod[j])&&inc[ct->numseq[a+1]][ct->numseq[j]]) {

										if (energy==w5[k-1]+v->f(k,a)+v->f(a+2,j-1)+
												ergcoaxflushbases(k,a,a+1,j,ct,data)+
												penalty(k,a,ct,data)+penalty(a+1,j,ct,data)

												+erg1(a+1,j,a+2,j-1,ct,data)) {

											registerbasepair(ct,a+1,j);
											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k,a,0,v->f(k,a),1);
											stack->push(a+2,j-1,0,v->f(a+2,j-1),1);
											found = true;
										}
									}
								}

								if (!found&&(energy==w5[k-1]+v->f(k,a)+v->f(a+2,j-1)+
											ergcoaxinterbases2(k,a,a+2,j-1,ct,data)+
											penalty(k,a,ct,data)+penalty(a+2,j-1,ct,data))) {

									if (k>1) stack->push(1,k-1,1,w5[k-1],0);
									stack->push(k,a,0,v->f(k,a),1);
									stack->push(a+2,j-1,0,v->f(a+2,j-1),1);
									found = true;
								}

								else if (!found&&(mod[k]||mod[a]||mod[a+2]||mod[j-1])) {

									if ((mod[k]||mod[a])&&(mod[a+2]||mod[j-1])&&inc[ct->numseq[k]][ct->numseq[a]]
											&&inc[ct->numseq[a+2]][ct->numseq[j-1]]) {

										if ((energy==w5[k-1]+v->f(k+1,a-1)+v->f(a+3,j-2)+
													ergcoaxinterbases2(k,a,a+2,j-1,ct,data)+
													penalty(k,a,ct,data)+penalty(a+2,j-1,ct,data)
													+erg1(k,a,k+1,a-1,ct,data)
													+erg1(a+2,j-1,a+3,j-2,ct,data))) {

											registerbasepair(ct,k,a);
											registerbasepair(ct,a+2,j-1);
											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k+1,a-1,0,v->f(k+1,a-1),1);
											stack->push(a+3,j-2,0,v->f(a+3,j-2),1);
											found = true;
										}
									}

									if (!found&&(mod[k]||mod[a])&&inc[ct->numseq[k]][ct->numseq[a]]) {

										if ((energy==w5[k-1]+v->f(k+1,a-1)+v->f(a+2,j-1)+
													ergcoaxinterbases2(k,a,a+2,j-1,ct,data)+
													penalty(k,a,ct,data)+penalty(a+2,j-1,ct,data)
													+erg1(k,a,k+1,a-1,ct,data))) {

											registerbasepair(ct,k,a);

											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k+1,a-1,0,v->f(k+1,a-1),1);
											stack->push(a+2,j-1,0,v->f(a+2,j-1),1);
											found = true;
										}
									}

									if (!found&&(mod[a+2]||mod[j-1])&&inc[ct->numseq[a+2]][ct->numseq[j-1]]) {

										if ((energy==w5[k-1]+v->f(k,a)+v->f(a+3,j-2)+
													ergcoaxinterbases2(k,a,a+2,j-1,ct,data)+
													penalty(k,a,ct,data)+penalty(a+2,j-1,ct,data)
													+erg1(a+2,j-1,a+3,j-2,ct,data))) {

											registerbasepair(ct,a+2,j-1);
											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k,a,0,v->f(k,a),1);
											stack->push(a+3,j-2,0,v->f(a+3,j-2),1);
											found = true;
										}
									}
								}

								if (!found&&(energy==w5[k-1]+v->f(k+1,a)+v->f(a+2,j)+
											ergcoaxinterbases1(k+1,a,a+2,j,ct,data)+
											penalty(k+1,a,ct,data)+penalty(a+2,j,ct,data))) {

									if (k>1) stack->push(1,k-1,1,w5[k-1],0);
									stack->push(k+1,a,0,v->f(k+1,a),1);
									stack->push(a+2,j,0,v->f(a+2,j),1);
									found = true;
								}

								else if (!found&&(mod[k+1]||mod[a]||mod[a+2]||mod[j])) {

									if ((mod[k+1]||mod[a])&&(mod[a+2]||mod[j])&&inc[ct->numseq[k+1]][ct->numseq[a]]
											&&inc[ct->numseq[a+2]][ct->numseq[j]]) {

										if ((energy==w5[k-1]+v->f(k+2,a-1)+v->f(a+3,j-1)+
													ergcoaxinterbases1(k+1,a,a+2,j,ct,data)+
													penalty(k+1,a,ct,data)+penalty(a+2,j,ct,data)
													+erg1(k+1,a,k+2,a-1,ct,data)
													+erg1(a+2,j,a+3,j-1,ct,data))) {

											registerbasepair(ct,k+1,a);
											registerbasepair(ct,a+2,j);
											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k+2,a-1,0,v->f(k+2,a-1),1);
											stack->push(a+3,j-1,0,v->f(a+3,j-1),1);
											found = true;
										}
									}

									if (!found&&(mod[k+1]||mod[a])&&inc[ct->numseq[k+1]][ct->numseq[a]]) {

										if ((energy==w5[k-1]+v->f(k+2,a-1)+v->f(a+2,j)+
													ergcoaxinterbases1(k+1,a,a+2,j,ct,data)+
													penalty(k+1,a,ct,data)+penalty(a+2,j,ct,data)
													+erg1(k+1,a,k+2,a-1,ct,data))) {

											registerbasepair(ct,k+1,a);

											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k+2,a-1,0,v->f(k+2,a-1),1);
											stack->push(a+2,j,0,v->f(a+2,j),1);
											found = true;
										}
									}

									if (!found&&(mod[a+2]||mod[j])&&inc[ct->numseq[a+2]][ct->numseq[j]]) {

										if ((energy==w5[k-1]+v->f(k+1,a)+v->f(a+3,j-1)+
													ergcoaxinterbases1(k+1,a,a+2,j,ct,data)+
													penalty(k+1,a,ct,data)+penalty(a+2,j,ct,data)
													+erg1(a+2,j,a+3,j-1,ct,data))) {

											registerbasepair(ct,a+2,j);
											if (k>1) stack->push(1,k-1,1,w5[k-1],0);
											stack->push(k+1,a,0,v->f(k+1,a),1);
											stack->push(a+3,j-1,0,v->f(a+3,j-1),1);
											found = true;
										}
									}
								}
								a++;
							}
						}
						k++;
					}
				}
				else {
					//j==n
					k=i+minloop;
					while(k<=j&&!found) {
						if (energy==w3[k]+v->f(i,k-1)+penalty(i,k-1,ct,data)) {
							stack->push(k,number,1,w3[k],0);
							stack->push(i,k-1,0,v->f(i,k-1),1);
							found=true;
						}
						else if (mod[i]||mod[k-1]&&inc[ct->numseq[i]][ct->numseq[k-1]]) {

							if (energy==w3[k]+v->f(i+1,k-2)+penalty(i,k-1,ct,data)
									+erg1(i,k-1,i+1,k-2,ct,data)) {

								registerbasepair(ct,i,k-1);
								stack->push(k,number,1,w3[k],0);
								stack->push(i+1,k-2,0,v->f(i+1,k-2),1);
								found=true;
							}
						}
						if (!found&&(energy==w3[k]+v->f(i,k-2)+penalty(i,k-2,ct,data)+
									erg4(k-2,i,k-1,1,ct,data,lfce[k-1]))) {

							stack->push(k,number,1,w3[k],0);
							stack->push(i,k-2,0,v->f(i,k-2),1);
							found=true;

						}
						else if (!found&&(mod[i]||mod[k-2])&&inc[ct->numseq[i]][ct->numseq[k-2]]) {

							if ((energy==w3[k]+v->f(i+1,k-3)+penalty(i,k-2,ct,data)+
										erg4(k-2,i,k-1,1,ct,data,lfce[k-1])
										+erg1(i,k-2,i+1,k-3,ct,data))) {

								registerbasepair(ct,i,k-2);
								stack->push(k,number,1,w3[k],0);
								stack->push(i+1,k-3,0,v->f(i+1,k-3),1);
								found=true;
							}
						}

						if (!found&&(energy==w3[k]+v->f(i+1,k-1)+penalty(i+1,k-1,ct,data)+
									erg4(k-1,i+1,i,2,ct,data,lfce[i]))) {

							stack->push(k,number,1,w3[k],0);
							stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
							found=true;

						}
						else if (!found&&(mod[i+1]||mod[k-1])&&inc[ct->numseq[i+1]][ct->numseq[k-1]]) {

							if ((energy==w3[k]+v->f(i+2,k-2)+penalty(i+1,k-1,ct,data)+
										erg4(k-1,i+1,i,2,ct,data,lfce[i])
										+erg1(i+1,k-1,i+2,k-2,ct,data))) {

								registerbasepair(ct,i+1,k-1);
								stack->push(k,number,1,w3[k],0);
								stack->push(i+2,k-2,0,v->f(i+2,k-2),1);
								found=true;

							}

						}
						if (!found&&(energy==w3[k]+v->f(i+1,k-2)+penalty(i+1,k-2,ct,data)+
									data->tstack[ct->numseq[k-2]][ct->numseq[i+1]]
									[ct->numseq[k-1]][ct->numseq[i]]+
									checknp(lfce[k-1],lfce[i])+ ct->SHAPEss_give_value(k-1)+ ct->SHAPEss_give_value(i))) {

							stack->push(k,number,1,w3[k],0);
							stack->push(i+1,k-2,0,v->f(i+1,k-2),1);
							found=true;


						}

						else if (!found&&(mod[i+1]||mod[k-2])&&inc[ct->numseq[i+1]][ct->numseq[k-2]]) {
							if ((energy==w3[k]+v->f(i+2,k-3)+penalty(i+1,k-2,ct,data)+
										data->tstack[ct->numseq[k-2]][ct->numseq[i+1]]
										[ct->numseq[k-1]][ct->numseq[i]]+
										checknp(lfce[k-1],lfce[i])+erg1(i+1,k-2,i+2,k-3,ct,data)+ 
										ct->SHAPEss_give_value(k-1)+ ct->SHAPEss_give_value(i))) {

								registerbasepair(ct,i+1,k-2);
								stack->push(k,number,1,w3[k],0);
								stack->push(i+2,k-3,0,v->f(i+2,k-3),1);
								found=true;
							}
						}
						if (!found) {
							//check for coaxial stacking:
							a = i+minloop+1;
							while (a<k&&!found) {

								if (energy==w3[k+1]+v->f(i,a)+v->f(a+1,k)+
										ergcoaxflushbases(i,a,a+1,k,ct,data)+
										penalty(i,a,ct,data)+penalty(a+1,k,ct,data)) {

									if (k<number) stack->push(k+1,number,1,w3[k+1],0);
									stack->push(i,a,0,v->f(i,a),1);
									stack->push(a+1,k,0,v->f(a+1,k),1);
									found =true;

								}
								else if (mod[i]||mod[a]||mod[a+1]||mod[k]) {

									if ((mod[i]||mod[a])&&(mod[a+1]||mod[k])&&inc[ct->numseq[i]][ct->numseq[a]]
											&&inc[ct->numseq[a+1]][ct->numseq[k]]) {
										if (energy==w3[k+1]+v->f(i+1,a-1)+v->f(a+2,k-1)+
												ergcoaxflushbases(i,a,a+1,k,ct,data)+
												penalty(i,a,ct,data)+penalty(a+1,k,ct,data)
												+erg1(i,a,i+1,a-1,ct,data)
												+erg1(a+1,k,a+2,k-1,ct,data)) {

											registerbasepair(ct,i,a);
											registerbasepair(ct,a+1,k);
											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i+1,a-1,0,v->f(i+1,a-1),1);
											stack->push(a+2,k-1,0,v->f(a+2,k-1),1);
											found =true;

										}	
									}

									if ((mod[i]||mod[a])&&!found&&inc[ct->numseq[i]][ct->numseq[a]]) {
										if (energy==w3[k+1]+v->f(i+1,a-1)+v->f(a+1,k)+
												ergcoaxflushbases(i,a,a+1,k,ct,data)+
												penalty(i,a,ct,data)+penalty(a+1,k,ct,data)
												+erg1(i,a,i+1,a-1,ct,data)) {

											registerbasepair(ct,i,a);

											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i+1,a-1,0,v->f(i+1,a-1),1);
											stack->push(a+1,k,0,v->f(a+1,k),1);
											found =true;

										}	
									}

									if (!found&&(mod[a+1]||mod[k])&&inc[ct->numseq[a+1]][ct->numseq[k]]) {
										if (energy==w3[k+1]+v->f(i,a)+v->f(a+2,k-1)+
												ergcoaxflushbases(i,a,a+1,k,ct,data)+
												penalty(i,a,ct,data)+penalty(a+1,k,ct,data)
												+erg1(a+1,k,a+2,k-1,ct,data)) {


											registerbasepair(ct,a+1,k);
											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i,a,0,v->f(i,a),1);
											stack->push(a+2,k-1,0,v->f(a+2,k-1),1);
											found =true;

										}	

									}


								}
								if (!found&&(energy==w3[k+1]+v->f(i,a)+v->f(a+2,k-1)+
											ergcoaxinterbases2(i,a,a+2,k-1,ct,data)+
											penalty(i,a,ct,data)+penalty(a+2,k-1,ct,data))) {

									if (k<number) stack->push(k+1,number,1,w3[k+1],0);
									stack->push(i,a,0,v->f(i,a),1);
									stack->push(a+2,k-1,0,v->f(a+2,k-1),1);
									found=true;

								}
								else if (!found&&(mod[i]||mod[a]||mod[a+2]||mod[k-1])) {
									if ((mod[i]||mod[a])&&(mod[a+2]||mod[k-1])&&inc[ct->numseq[i]][ct->numseq[a]]
											&&inc[ct->numseq[a+2]][ct->numseq[k-1]]) {
										if ((energy==w3[k+1]+v->f(i+1,a-1)+v->f(a+3,k-2)+
													ergcoaxinterbases2(i,a,a+2,k-1,ct,data)+
													penalty(i,a,ct,data)+penalty(a+2,k-1,ct,data)
													+erg1(i,a,i+1,a-1,ct,data)
													+erg1(a+2,k-1,a+3,k-2,ct,data))) {

											registerbasepair(ct,i,a);
											registerbasepair(ct,a+2,k-1);
											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i+1,a-1,0,v->f(i+1,a-1),1);
											stack->push(a+3,k-2,0,v->f(a+3,k-2),1);
											found=true;

										}
									}

									if ((mod[i]||mod[a])&&!found&&inc[ct->numseq[i]][ct->numseq[a]]) {
										if ((energy==w3[k+1]+v->f(i+1,a-1)+v->f(a+2,k-1)+
													ergcoaxinterbases2(i,a,a+2,k-1,ct,data)+
													penalty(i,a,ct,data)+penalty(a+2,k-1,ct,data)
													+erg1(i,a,i+1,a-1,ct,data)
										    )) {

											registerbasepair(ct,i,a);

											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i+1,a-1,0,v->f(i+1,a-1),1);
											stack->push(a+2,k-1,0,v->f(a+2,k-1),1);
											found=true;

										}
									}

									if ((!found)&&(mod[a+2]||mod[k-1])&&inc[ct->numseq[a+2]][ct->numseq[k-1]]) {
										if ((energy==w3[k+1]+v->f(i,a)+v->f(a+3,k-2)+
													ergcoaxinterbases2(i,a,a+2,k-1,ct,data)+
													penalty(i,a,ct,data)+penalty(a+2,k-1,ct,data)
													+erg1(a+2,k-1,a+3,k-2,ct,data))) {


											registerbasepair(ct,a+2,k-1);
											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i,a,0,v->f(i,a),1);
											stack->push(a+3,k-2,0,v->f(a+3,k-2),1);
											found=true;

										}
									}
								}

								if (!found&&(energy==w3[k+1]+v->f(i+1,a)+v->f(a+2,k)+
											ergcoaxinterbases1(i+1,a,a+2,k,ct,data)+
											penalty(i+1,a,ct,data)+penalty(a+2,k,ct,data)
									    )) {

									if (k<number) stack->push(k+1,number,1,w3[k+1],0);
									stack->push(i+1,a,0,v->f(i+1,a),1);
									stack->push(a+2,k,0,v->f(a+2,k),1);
									found=true;

								}
								else if (!found&&(mod[i+1]||mod[a]||mod[a+2]||mod[k])) {
									if ((mod[i+1]||mod[a])&&(mod[a+2]||mod[k])
											&&inc[ct->numseq[i+1]][ct->numseq[a]]
											&&inc[ct->numseq[a+2]][ct->numseq[k]]) {

										if ((energy==w3[k+1]+v->f(i+2,a-1)+v->f(a+3,k-1)+
													ergcoaxinterbases1(i+1,a,a+2,k,ct,data)+
													penalty(i+1,a,ct,data)+penalty(a+2,k,ct,data)
													+erg1(i+1,a,i+2,a-1,ct,data)
													+erg1(a+2,k,a+3,k-1,ct,data))) {

											registerbasepair(ct,i+1,a);
											registerbasepair(ct,a+2,k);
											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i+2,a-1,0,v->f(i+2,a-1),1);
											stack->push(a+3,k-1,0,v->f(a+3,k-1),1);
											found=true;

										}	
									}

									if ((mod[i+1]||mod[a])&&(!found)&&inc[ct->numseq[i+1]][ct->numseq[a]]) {

										if ((energy==w3[k+1]+v->f(i+2,a-1)+v->f(a+2,k)+
													ergcoaxinterbases1(i+1,a,a+2,k,ct,data)+
													penalty(i+1,a,ct,data)+penalty(a+2,k,ct,data)
													+erg1(i+1,a,i+2,a-1,ct,data))) {

											registerbasepair(ct,i+1,a);

											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i+2,a-1,0,v->f(i+2,a-1),1);
											stack->push(a+2,k,0,v->f(a+2,k),1);
											found=true;

										}	
									}

									if ((!found)&&(mod[a+2]||mod[k])&&inc[ct->numseq[a+2]][ct->numseq[k]]) {

										if ((energy==w3[k+1]+v->f(i+1,a)+v->f(a+3,k-1)+
													ergcoaxinterbases1(i+1,a,a+2,k,ct,data)+
													penalty(i+1,a,ct,data)+penalty(a+2,k,ct,data)
													+erg1(a+2,k,a+3,k-1,ct,data))) {


											registerbasepair(ct,a+2,k);
											if (k<number) stack->push(k+1,number,1,w3[k+1],0);
											stack->push(i+1,a,0,v->f(i+1,a),1);
											stack->push(a+3,k-1,0,v->f(a+3,k-1),1);
											found=true;

										}	
									}
								}
								a++;
							}
						}
						k++;
					}
				}
#else
              			  k=i+1;
			  while(k<j&&!found) {
					if (energy==we->f(i,k)+we->f(k+1,j)) {
						stack->push(i,k,1,we->f(i,k),0);
						stack->push(k+1,j,1,we->f(k+1,j),0);
						found = true;
					}

					//also check for coaxial stacking:
					else if (energy==v->f(i,k)+v->f(k+1,j)+penalty(i,k,ct,data)+
					         penalty(k+1,j,ct,data)
					         +ergcoaxflushbases(i,k,k+1,j,ct,data)) {

					  stack->push(i,k,0,v->f(i,k),1);
					  stack->push(k+1,j,0,v->f(k+1,j),1);
					  found = true;
					}

					else if (mod[i]||mod[k]||mod[k+1]||mod[j]) {

					  if ((mod[i]||mod[k])&&(mod[k+1]||mod[j])&&inc[ct->numseq[i]][ct->numseq[k]]
					      &&inc[ct->numseq[k+1]][ct->numseq[j]]) {

					    if (energy==v->f(i+1,k-1)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
						penalty(k+1,j,ct,data)
						+ergcoaxflushbases(i,k,k+1,j,ct,data)
						+erg1(i,k,i+1,k-1,ct,data)
						+erg1(k+1,j,k+2,j-1,ct,data)) {

					      registerbasepair(ct,i,k);
					      registerbasepair(ct,k+1,j);
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					      found = true;
					    }		
					  }

					  if ((mod[i]||mod[k])&&(!found)&&inc[ct->numseq[i]][ct->numseq[k]]) {
					    
					    if (energy==v->f(i+1,k-1)+v->f(k+1,j)+penalty(i,k,ct,data)+
						penalty(k+1,j,ct,data)
						+ergcoaxflushbases(i,k,k+1,j,ct,data)
						+erg1(i,k,i+1,k-1,ct,data)) {

					      registerbasepair(ct,i,k);
					      
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+1,j,0,v->f(k+1,j),1);
					      found = true;
					    }		
					  }
					  if ((!found)&&(mod[k+1]||mod[j])&&inc[ct->numseq[k+1]][ct->numseq[j]]) {

					    if (energy==v->f(i,k)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
						penalty(k+1,j,ct,data)
						+ergcoaxflushbases(i,k,k+1,j,ct,data)
						+erg1(k+1,j,k+2,j-1,ct,data)) {

							
					      registerbasepair(ct,k+1,j);
					      stack->push(i,k,0,v->f(i,k),1);
					      stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					      found = true;
					    }		
					  }
					}
					if (!found&&(energy==v->f(i,k)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
					             penalty(k+2,j-1,ct,data)
					             +ergcoaxinterbases2(i,k,k+2,j-1,ct,data))) {

					  stack->push(i,k,0,v->f(i,k),1);
					  stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					  found = true;
					}
					else if (!found&&(mod[i]||mod[k]||mod[k+2]||mod[j-1])) {

					  if ((mod[i]||mod[k])&&(mod[k+2]||mod[j-1])&&inc[ct->numseq[i]][ct->numseq[k]]
					      &&inc[ct->numseq[k+2]][ct->numseq[j-1]]) {

					    if ((energy==v->f(i+1,k-1)+v->f(k+3,j-2)+penalty(i,k,ct,data)+
						 penalty(k+2,j-1,ct,data)
						 +ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
						 +erg1(i,k,i+1,k-1,ct,data)
						 +erg1(k+2,j-1,k+3,j-2,ct,data))) {

					      registerbasepair(ct,i,k);
					      registerbasepair(ct,k+2,j-1);
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+3,j-2,0,v->f(k+3,j-2),1);
					      found = true;
					    }	
					  }

					  if ((mod[i]||mod[k])&&(!found)&&inc[ct->numseq[i]][ct->numseq[k]]) {
					    
					    if ((energy==v->f(i+1,k-1)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
						 penalty(k+2,j-1,ct,data)
						 +ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
						 +erg1(i,k,i+1,k-1,ct,data))) {
					      
					      registerbasepair(ct,i,k);
								
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					      found = true;
					    }	
					  }

					  if ((!found)&&(mod[k+2]||mod[j-1])&&inc[ct->numseq[k+2]][ct->numseq[j-1]]) {
					    
					    if ((energy==v->f(i,k)+v->f(k+3,j-2)+penalty(i,k,ct,data)+
						 penalty(k+2,j-1,ct,data)+
						 ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
						 +erg1(k+2,j-1,k+3,j-2,ct,data))) {

								
					      registerbasepair(ct,k+2,j-1);
					      stack->push(i,k,0,v->f(i,k),1);
					      stack->push(k+3,j-2,0,v->f(k+3,j-2),1);
					      found = true;
					    }	
					  }
					}
					if (!found&&(energy==v->f(i+1,k)+v->f(k+2,j)+penalty(i+1,k,ct,data)+
					             penalty(k+2,j,ct,data)+
					             ergcoaxinterbases1(i+1,k,k+2,j,ct,data))) {

					  stack->push(i+1,k,0,v->f(i+1,k),1);
					  stack->push(k+2,j,0,v->f(k+2,j),1);
					  found = true;
					}

					else if (!found&&(mod[i+1]||mod[k]||mod[k+2]||mod[j])) {

					  if ((mod[i+1]||mod[k])&&(mod[k+2]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[k]]
					      &&inc[ct->numseq[k+2]][ct->numseq[j]]) {
					    if ((energy==v->f(i+2,k-1)+v->f(k+3,j-1)+penalty(i+1,k,ct,data)+
						 penalty(k+2,j,ct,data)+
						 ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
						 +erg1(i+1,k,i+2,k-1,ct,data)
						 +erg1(k+2,j,k+3,j-1,ct,data))) {

					      registerbasepair(ct,i+1,k);
					      registerbasepair(ct,k+2,j);
					      stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
					      stack->push(k+3,j-1,0,v->f(k+3,j-1),1);
					      found = true;
					    }
					  }

					  if ((mod[i+1]||mod[k])&&(!found)&&inc[ct->numseq[i+1]][ct->numseq[k]]) {
					    if ((energy==v->f(i+2,k-1)+v->f(k+2,j)+penalty(i+1,k,ct,data)+
						 penalty(k+2,j,ct,data)+
						 ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
						 +erg1(i+1,k,i+2,k-1,ct,data))) {

					      registerbasepair(ct,i+1,k);
								
					      stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
					      stack->push(k+2,j,0,v->f(k+2,j),1);
					      found = true;
					    }
					  }

					  if ((!found)&&(mod[k+2]||mod[j])&&inc[ct->numseq[k+2]][ct->numseq[j]]) {
					    if ((energy==v->f(i+1,k)+v->f(k+3,j-1)+penalty(i+1,k,ct,data)+
						 penalty(k+2,j,ct,data)+
						 ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
								
						 +erg1(k+2,j,k+3,j-1,ct,data))) {

					      
					      registerbasepair(ct,k+2,j);
					      stack->push(i+1,k,0,v->f(i+1,k),1);
					      stack->push(k+3,j-1,0,v->f(k+3,j-1),1);
					      found = true;
					    }
					  }
					}
					k++;
			  }

#endif                  


			}
			else {
				//open == 0, multiloop
				k=i+1;
				while(k<j&&!found) {
					if (energy==w->f(i,k)+w->f(k+1,j)) {
						stack->push(i,k,0,w->f(i,k),0);
						stack->push(k+1,j,0,w->f(k+1,j),0);
						found = true;
					}

					//also check for coaxial stacking:
					else if (energy==v->f(i,k)+v->f(k+1,j)+penalty(i,k,ct,data)+
							penalty(k+1,j,ct,data)+2*data->eparam[10]+
							ergcoaxflushbases(i,k,k+1,j,ct,data)) {

						stack->push(i,k,0,v->f(i,k),1);
						stack->push(k+1,j,0,v->f(k+1,j),1);
						found = true;
					}

					else if (mod[i]||mod[k]||mod[k+1]||mod[j]) {

						if ((mod[i]||mod[k])&&(mod[k+1]||mod[j])&&inc[ct->numseq[i]][ct->numseq[k]]
								&&inc[ct->numseq[k+1]][ct->numseq[j]]) {

							if (energy==v->f(i+1,k-1)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
									penalty(k+1,j,ct,data)+2*data->eparam[10]+
									ergcoaxflushbases(i,k,k+1,j,ct,data)
									+erg1(i,k,i+1,k-1,ct,data)
									+erg1(k+1,j,k+2,j-1,ct,data)) {

								registerbasepair(ct,i,k);
								registerbasepair(ct,k+1,j);
								stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
								stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
								found = true;
							}		
						}

						if ((mod[i]||mod[k])&&(!found)&&inc[ct->numseq[i]][ct->numseq[k]]) {

							if (energy==v->f(i+1,k-1)+v->f(k+1,j)+penalty(i,k,ct,data)+
									penalty(k+1,j,ct,data)+2*data->eparam[10]+
									ergcoaxflushbases(i,k,k+1,j,ct,data)
									+erg1(i,k,i+1,k-1,ct,data)) {

								registerbasepair(ct,i,k);

								stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
								stack->push(k+1,j,0,v->f(k+1,j),1);
								found = true;
							}		
						}
						if ((!found)&&(mod[k+1]||mod[j])&&inc[ct->numseq[k+1]][ct->numseq[j]]) {

							if (energy==v->f(i,k)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
									penalty(k+1,j,ct,data)+2*data->eparam[10]+
									ergcoaxflushbases(i,k,k+1,j,ct,data)
									+erg1(k+1,j,k+2,j-1,ct,data)) {


								registerbasepair(ct,k+1,j);
								stack->push(i,k,0,v->f(i,k),1);
								stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
								found = true;
							}		
						}
					}
					if (!found&&(energy==v->f(i,k)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
								penalty(k+2,j-1,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
								ergcoaxinterbases2(i,k,k+2,j-1,ct,data))) {

						stack->push(i,k,0,v->f(i,k),1);
						stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
						found = true;
					}
					else if (!found&&(mod[i]||mod[k]||mod[k+2]||mod[j-1])) {

						if ((mod[i]||mod[k])&&(mod[k+2]||mod[j-1])&&inc[ct->numseq[i]][ct->numseq[k]]
								&&inc[ct->numseq[k+2]][ct->numseq[j-1]]) {

							if ((energy==v->f(i+1,k-1)+v->f(k+3,j-2)+penalty(i,k,ct,data)+
										penalty(k+2,j-1,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
										ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
										+erg1(i,k,i+1,k-1,ct,data)
										+erg1(k+2,j-1,k+3,j-2,ct,data))) {

								registerbasepair(ct,i,k);
								registerbasepair(ct,k+2,j-1);
								stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
								stack->push(k+3,j-2,0,v->f(k+3,j-2),1);
								found = true;
							}	
						}

						if ((mod[i]||mod[k])&&(!found)&&inc[ct->numseq[i]][ct->numseq[k]]) {

							if ((energy==v->f(i+1,k-1)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
										penalty(k+2,j-1,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
										ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
										+erg1(i,k,i+1,k-1,ct,data))) {

								registerbasepair(ct,i,k);

								stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
								stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
								found = true;
							}	
						}

						if ((!found)&&(mod[k+2]||mod[j-1])&&inc[ct->numseq[k+2]][ct->numseq[j-1]]) {

							if ((energy==v->f(i,k)+v->f(k+3,j-2)+penalty(i,k,ct,data)+
										penalty(k+2,j-1,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
										ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
										+erg1(k+2,j-1,k+3,j-2,ct,data))) {


								registerbasepair(ct,k+2,j-1);
								stack->push(i,k,0,v->f(i,k),1);
								stack->push(k+3,j-2,0,v->f(k+3,j-2),1);
								found = true;
							}	
						}
					}
					if (!found&&(energy==v->f(i+1,k)+v->f(k+2,j)+penalty(i+1,k,ct,data)+
								penalty(k+2,j,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
								ergcoaxinterbases1(i+1,k,k+2,j,ct,data))) {

						stack->push(i+1,k,0,v->f(i+1,k),1);
						stack->push(k+2,j,0,v->f(k+2,j),1);
						found = true;
					}

					else if (!found&&(mod[i+1]||mod[k]||mod[k+2]||mod[j])) {

						if ((mod[i+1]||mod[k])&&(mod[k+2]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[k]]
								&&inc[ct->numseq[k+2]][ct->numseq[j]]) {
							if ((energy==v->f(i+2,k-1)+v->f(k+3,j-1)+penalty(i+1,k,ct,data)+
										penalty(k+2,j,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
										ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
										+erg1(i+1,k,i+2,k-1,ct,data)
										+erg1(k+2,j,k+3,j-1,ct,data))) {

								registerbasepair(ct,i+1,k);
								registerbasepair(ct,k+2,j);
								stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
								stack->push(k+3,j-1,0,v->f(k+3,j-1),1);
								found = true;
							}
						}

						if ((mod[i+1]||mod[k])&&(!found)&&inc[ct->numseq[i+1]][ct->numseq[k]]) {
							if ((energy==v->f(i+2,k-1)+v->f(k+2,j)+penalty(i+1,k,ct,data)+
										penalty(k+2,j,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
										ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
										+erg1(i+1,k,i+2,k-1,ct,data))) {

								registerbasepair(ct,i+1,k);

								stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
								stack->push(k+2,j,0,v->f(k+2,j),1);
								found = true;
							}
						}

						if ((!found)&&(mod[k+2]||mod[j])&&inc[ct->numseq[k+2]][ct->numseq[j]]) {
							if ((energy==v->f(i+1,k)+v->f(k+3,j-1)+penalty(i+1,k,ct,data)+
										penalty(k+2,j,ct,data)+2*data->eparam[10]+2*data->eparam[6]+
										ergcoaxinterbases1(i+1,k,k+2,j,ct,data)

										+erg1(k+2,j,k+3,j-1,ct,data))) {


								registerbasepair(ct,k+2,j);
								stack->push(i+1,k,0,v->f(i+1,k),1);
								stack->push(k+3,j-1,0,v->f(k+3,j-1),1);
								found = true;
							}
						}
					}
					k++;
				}

				if (ct->intermolecular) {
					k=i+1;
					while(k<j&&!found) {
						if (energy==w2->f(i,k)+w2->f(k+1,j)) {
							stack->push(i,k,0,w2->f(i,k),0);
							stack->push(k+1,j,0,w2->f(k+1,j),0);
							found = true;
						}

						k++;
					}


					if (i!=number) {
						if (!(fce->f(i,i)&INTER)) {
							if (energy==wmb2->f(i+1,j) ) {
								stack->push(i+1,j,0,w2->f(i+1,j),0); 
								found = true;
							}
						}
						else  if (energy==wmb2->f(i+1,j) + data->init - INFINITE_ENERGY) {
							stack->push(i+1,j,0,w2->f(i+1,j),0);
							found = true;
						}
					}
					if ((j!=number+1)&&(!found)) {
						if (!(fce->f(j,j)&INTER)) {
							if (energy==wmb2->f(i,j-1)) {
								stack->push(i,j-1,0,wmb2->f(i,j-1),0);
								found = true;
							}
						}
						else if (energy==wmb2->f(i,j-1) +data->init-INFINITE_ENERGY) {
							stack->push(i,j-1,0,wmb2->f(i,j-1),0);
							found = true;
						}
					}
				}
			}
			if (i==j) found = true;

			if (!found) {
				errmsg (100,2);
				tracebackerror=1;
			}
		}

		else if (pair==1) {
			//energy = v->f(i,j) So we have a basepair


			while (energy==v->f(i,j)) {
				//record the found base pair
				registerbasepair(ct,i,j);

				//check for stacked pair:
				if (energy==erg1(i,j,i+1,j-1,ct,data)+v->f(i+1,j-1)&&
						i!=number&&j!=number+1) {
					i++;
					j--;
					energy = v->f(i,j);
				}
				else break;
			}

			//now past the helical region:

			//check for hairpin loop:
			//if it is a hairpin loop, stop and return to the stack, otherwise
			//continue to define the loop closed by i and j
			if (energy!=erg3(i,j,ct,data,fce->f(i,j))) {

				found = false; //use the flag found to keep track as to whether a 
				//multiloop or exterior loop was found

				if (i==number) {
					//place the 5' fragment on the stack if big enough to have structure:
					if (energy==penalty(i,j,ct,data)+w5[j-number-1]) {
						if (j-number-1>minloop+1) {
							stack->push(1,j-number-1,1,w5[j-number-1],0);
							found = true;
						}
						else found=true;
					}
					if (!found&&energy==penalty(i,j,ct,data)+w5[j-number-2]+
							erg4(i,j,j-1,2,ct,data,lfce[j-1])) {
						if (j-number-2>minloop+1) {
							stack->push(1,j-number-2,1,w5[j-number-2],0);
							found = true;
						}
						else found=true;
					}

					//now consider stacking to the 3' fragment:
					//first consider coaxial stacking to the 5' fragment:
					a = j-number-minloop-2;
					while (a>0&&!found) {
						if (energy==w5[a-1]+penalty(i,j,ct,data)+
								penalty(a,j-number-1,ct,data)+v->f(a,j-number-1)+
								ergcoaxflushbases(a,j-number-1,j-number,i,ct,data)) {

							if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
							if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
							stack->push(a,j-number-1,0,v->f(a,j-number-1),1);
							found = true;
						}

						else if (mod[a]||mod[j-number-1]&&inc[ct->numseq[a]][ct->numseq[j-number-1]]) {

							if (energy==w5[a-1]+penalty(i,j,ct,data)+
									penalty(a,j-number-1,ct,data)+v->f(a+1,j-number-2)+
									ergcoaxflushbases(a,j-number-1,j-number,i,ct,data)
									+erg1(a,j-number-1,a+1,j-number-2,ct,data)) {

								registerbasepair(ct,a,j-number-1);
								if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
								if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
								stack->push(a+1,j-number-2,0,v->f(a+1,j-number-2),1);
								found = true;
							}

						}

						if (!found&&(energy==w5[a-1]+penalty(i,j,ct,data)+
									penalty(a+1,j-number-2,ct,data)+v->f(a+1,j-number-2)+
									ergcoaxinterbases1(a+1,j-number-2,j-number,i,ct,data))) {

							if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
							if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
							stack->push(a+1,j-number-2,0,v->f(a+1,j-number-2),1);
							found = true;

						}
						else if (!found&&(mod[a+1]||mod[j-number-2])&&inc[ct->numseq[i+1]][ct->numseq[j-number-2]]) {
							if ((energy==w5[a-1]+penalty(i,j,ct,data)+
										penalty(a+1,j-number-2,ct,data)+v->f(a+2,j-number-3)+
										ergcoaxinterbases1(a+1,j-number-2,j-number,i,ct,data)
										+erg1(a+1,j-number-2,a+2,j-number-3,ct,data))) {

								registerbasepair(ct,a+1,j-number-2);
								if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
								if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
								stack->push(a+2,j-number-3,0,v->f(a+2,j-number-3),1);
								found = true;

							}
						}
						if (!found&&i<number) {
							if (energy==w5[a-1]+penalty(i,j,ct,data)+
									penalty(a,j-number-2,ct,data)+v->f(a,j-number-2)+
									ergcoaxinterbases2(a,j-number-2,j-number,i,ct,data)) {

								if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
								if (i+2<number-minloop-1) stack->push(i+2,number,1,w3[i+2],0);
								stack->push(a,j-number-2,0,v->f(a,j-number-2),1);
								found = true;

							}
							else if (mod[a]||mod[j-number-2]&&inc[ct->numseq[a]][ct->numseq[j-number-2]]) {

								if (energy==w5[a-1]+penalty(i,j,ct,data)+
										penalty(a,j-number-2,ct,data)+v->f(a+1,j-number-3)+
										ergcoaxinterbases2(a,j-number-2,j-number,i,ct,data)
										+erg1(a,j-number-2,a+1,j-number-3,ct,data)) {

									registerbasepair(ct,a,j-number-2);
									if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
									if (i+2<number-minloop-1) stack->push(i+2,number,1,w3[i+2],0);
									stack->push(a+1,j-number-3,0,v->f(a+1,j-number-3),1);
									found = true;
								}
							}
						}
						a--;
					}
				}
				else {
					//this means that (j>i+2*minloop+2) 

					//check for multiloop closed by i and j

					//check for the case where the bifuraction is in a multibranch loop:

					if (energy==wmb->f(i+1,j-1)+data->eparam[10]+data->eparam[5]+penalty(i,j,ct,data)) {

						stack->push(i+1,j-1,0,wmb->f(i+1,j-1),2);
						found = true;

					}
					else if (energy==wmb->f(i+2,j-1)+penalty(i,j,ct,data)+
							data->eparam[10]+data->eparam[5]+data->eparam[6]+
							erg4(i,j,i+1,1,ct,data,lfce[i+1])) {

						stack->push(i+2,j-1,0,wmb->f(i+2,j-1),2);
						found = true;

					}
					else if (energy==wmb->f(i+1,j-2)+penalty(i,j,ct,data)+
							data->eparam[10]+data->eparam[5]+data->eparam[6]+
							erg4(i,j,j-1,2,ct,data,lfce[j-1])) {

						stack->push(i+1,j-2,0,wmb->f(i+1,j-2),2);
						found = true;

					}
					else if (energy==wmb->f(i+2,j-2)+penalty(i,j,ct,data)+
							data->eparam[10]+data->eparam[5]+2*data->eparam[6]+
							data->tstkm[ct->numseq[i]][ct->numseq[j]]
							[ct->numseq[i+1]][ct->numseq[j-1]]
							+checknp(lfce[i+1],lfce[j-1])+ ct->SHAPEss_give_value(i+1)+ ct->SHAPEss_give_value(j-1)) {

						stack->push(i+2,j-2,0,wmb->f(i+2,j-2),2);
						found = true;
					}
					k=i+1;
					while (k<j&&!found) {
						//if (k!=number) {

						if (k+1<j-1&&energy==penalty(i,j,ct,data)+
								penalty(i+1,k,ct,data)+v->f(i+1,k)+
								ergcoaxflushbases(j,i,i+1,k,ct,data)+w->f(k+1,j-1)+
								data->eparam[5]+2*data->eparam[10]) {

							stack->push(i+1,k,0,v->f(i+1,k),1);
							stack->push(k+1,j-1,0,w->f(k+1,j-1),0);
							found = true;
						}

						else if (k+1<j-1&&mod[i+1]||mod[k]&&inc[ct->numseq[i+1]][ct->numseq[k]]) {

							if (energy==penalty(i,j,ct,data)+
									penalty(i+1,k,ct,data)+v->f(i+2,k-1)+
									ergcoaxflushbases(j,i,i+1,k,ct,data)+w->f(k+1,j-1)+
									data->eparam[5]+2*data->eparam[10]
									+erg1(i+1,k,i+2,k-1,ct,data)) {

								registerbasepair(ct,i+1,k);
								stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
								stack->push(k+1,j-1,0,w->f(k+1,j-1),0);
								found = true;
							}
						}

						if (!found&&(k+1<j-2&&energy==penalty(i,j,ct,data)+
									penalty(i+2,k,ct,data)+v->f(i+2,k)+
									ergcoaxinterbases1(j,i,i+2,k,ct,data)+w->f(k+1,j-2)+
									data->eparam[5]+2*data->eparam[10]+2*data->eparam[6])) {

							stack->push(i+2,k,0,v->f(i+2,k),1);
							stack->push(k+1,j-2,0,w->f(k+1,j-2),0);
							found = true;
						}
						else if (!found&&(k+1<j-2&&(mod[i+2]||mod[k]))&&inc[ct->numseq[i+2]][ct->numseq[k]]) {

							if (energy==penalty(i,j,ct,data)+
									penalty(i+2,k,ct,data)+v->f(i+3,k-1)+
									ergcoaxinterbases1(j,i,i+2,k,ct,data)+w->f(k+1,j-2)+
									data->eparam[5]+2*data->eparam[10]+2*data->eparam[6]
									+erg1(i+2,k,i+3,k-1,ct,data)) {

								registerbasepair(ct,i+2,k);
								stack->push(i+3,k-1,0,v->f(i+3,k-1),1);
								stack->push(k+1,j-2,0,w->f(k+1,j-2),0);
								found = true;
							}
						}
						if (!found&&k+1<j-1&&energy==penalty(i,j,ct,data)+
								penalty(i+2,k-1,ct,data)+v->f(i+2,k-1)+
								ergcoaxinterbases2(j,i,i+2,k-1,ct,data)+w->f(k+1,j-1)+
								data->eparam[5]+2*data->eparam[10]+2*data->eparam[6]) {

							stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
							stack->push(k+1,j-1,0,w->f(k+1,j-1),0);
							found = true;


						}
						else if (!found&&(k+1<j-1)&&(mod[i+2]||mod[k-1])&&inc[ct->numseq[i+2]][ct->numseq[k-1]]) {

							if (energy==penalty(i,j,ct,data)+
									penalty(i+2,k-1,ct,data)+v->f(i+3,k-2)+
									ergcoaxinterbases2(j,i,i+2,k-1,ct,data)+w->f(k+1,j-1)+
									data->eparam[5]+2*data->eparam[10]+2*data->eparam[6]
									+erg1(i+2,k-1,i+3,k-2,ct,data)) {

								registerbasepair(ct,i+2,k-1);
								stack->push(i+3,k-2,0,v->f(i+3,k-2),1);
								stack->push(k+1,j-1,0,w->f(k+1,j-1),0);
								found = true;
							}
						}

						if (!found&&i+1<k-1&&energy==penalty(i,j,ct,data)+
								penalty(k,j-1,ct,data)+v->f(k,j-1)+
								ergcoaxflushbases(k,j-1,j,i,ct,data)+w->f(i+1,k-1)+
								data->eparam[5]+2*data->eparam[10]) {

							stack->push(k,j-1,0,v->f(k,j-1),1);
							stack->push(i+1,k-1,0,w->f(i+1,k-1),0);
							found = true;

						}
						else if (!found&&i+1<k-1&&(mod[k]||mod[j-1])&&inc[ct->numseq[k]][ct->numseq[j-1]]) {

							if (energy==penalty(i,j,ct,data)+
									penalty(k,j-1,ct,data)+v->f(k+1,j-2)+
									ergcoaxflushbases(k,j-1,j,i,ct,data)+w->f(i+1,k-1)+
									data->eparam[5]+2*data->eparam[10]
									+erg1(k,j-1,k+1,j-2,ct,data)) {

								registerbasepair(ct,k,j-1);
								stack->push(k+1,j-2,0,v->f(k+1,j-2),1);
								stack->push(i+1,k-1,0,w->f(i+1,k-1),0);
								found = true;
							}
						}

						if (!found&&i+2<k-1&&energy==penalty(i,j,ct,data)+
								penalty(k,j-2,ct,data)+v->f(k,j-2)+
								ergcoaxinterbases2(k,j-2,j,i,ct,data)+w->f(i+2,k-1)+
								data->eparam[5]+2*data->eparam[10]+2*data->eparam[6]
						   ){

							stack->push(k,j-2,0,v->f(k,j-2),1);
							stack->push(i+2,k-1,0,w->f(i+2,k-1),0);
							found = true;
						}

						else if (!found&&i+2<k-1&&(mod[k]||mod[j-2])&&inc[ct->numseq[k]][ct->numseq[j-2]]) {
							if (energy==penalty(i,j,ct,data)+
									penalty(k,j-2,ct,data)+v->f(k+1,j-3)+
									ergcoaxinterbases2(k,j-2,j,i,ct,data)+w->f(i+2,k-1)+
									data->eparam[5]+2*data->eparam[10]+2*data->eparam[6]
									+erg1(k,j-2,k+1,j-3,ct,data)){

								registerbasepair(ct,k,j-2);
								stack->push(k+1,j-3,0,v->f(k+1,j-3),1);
								stack->push(i+2,k-1,0,w->f(i+2,k-1),0);
								found = true;
							}
						}

						if (!found&&i+1<k-1&&energy==penalty(i,j,ct,data)+
								penalty(k+1,j-2,ct,data)+v->f(k+1,j-2)+
								ergcoaxinterbases1(k+1,j-2,j,i,ct,data)+w->f(i+1,k-1)+
								data->eparam[5]+2*data->eparam[10]+2*data->eparam[6]
						   ) {

							stack->push(k+1,j-2,0,v->f(k+1,j-2),1);
							stack->push(i+1,k-1,0,w->f(i+1,k-1),0);
							found = true;
						}
						else if (!found&&i+1<k-1&&(mod[k+1]||mod[j-2])&&inc[ct->numseq[k+1]][ct->numseq[j-2]]) {

							if (energy==penalty(i,j,ct,data)+
									penalty(k+1,j-2,ct,data)+v->f(k+2,j-3)+
									ergcoaxinterbases1(k+1,j-2,j,i,ct,data)+w->f(i+1,k-1)+
									data->eparam[5]+2*data->eparam[10]+2*data->eparam[6]
									+erg1(k+1,j-2,k+2,j-3,ct,data)) {

								registerbasepair(ct,k+1,j-2);
								stack->push(k+2,j-3,0,v->f(k+2,j-3),1);
								stack->push(i+1,k-1,0,w->f(i+1,k-1),0);
								found = true;
							}
						}

						//}
						if (k==number&&!found) {
							//else {
							//k==number
							if (energy==w3[i+1]+w5[j-number-1]+penalty(i,j,ct,data)) {


								if (i<number-minloop-2) stack->push(i+1,number,1,w3[i+1],0);
								if (j-number>minloop+3) stack->push(1,j-number-1,1,w5[j-number-1],0);
								found = true;
							}
							else if (energy==w3[i+2]+w5[j-number-1]+penalty(i,j,ct,data) +
									erg4(i,j,i+1,1,ct,data,lfce[i+1])) {

								if (i<number-minloop-3) stack->push(i+2,number,1,w3[i+2],0);
								if (j-number>minloop+3) stack->push(1,j-number-1,1,w5[j-number-1],0);
								found = true;

							}

							else if (j-number-2>-1) {

								if (energy==w3[i+1]+w5[j-number-2]+penalty(i,j,ct,data) +
										erg4(i,j,j-1,2,ct,data,lfce[j-1])) {

									if (i<number-minloop-2) stack->push(i+1,number,1,w3[i+1],0);
									if (j-number>minloop+4) stack->push(1,j-number-2,1,w5[j-number-2],0);
									found = true;

								}
								else if (energy==w3[i+2]+w5[j-number-2]+penalty(i,j,ct,data) +
										data->tstack[ct->numseq[i]][ct->numseq[j]]
										[ct->numseq[i+1]][ct->numseq[j-1]]
										+checknp(lfce[i+1],lfce[j-1])+ ct->SHAPEss_give_value(i+1)+ ct->SHAPEss_give_value(j-1)) {

									if (i<number-minloop-3) stack->push(i+2,number,1,w3[i+2],0);
									if (j-number>minloop+4) stack->push(1,j-number-2,1,w5[j-number-2],0);
									found = true;

								}
							}

							//also consider coaxial stacking:

							//first consider coaxial stacking to the 5' fragment:
							a = j-number-minloop-2;
							while (a>0&&!found) {
								if (energy==w5[a-1]+w3[i+1]+penalty(i,j,ct,data)+
										penalty(a,j-number-1,ct,data)+v->f(a,j-number-1)+
										ergcoaxflushbases(a,j-number-1,j-number,i,ct,data)) {

									if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
									if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
									stack->push(a,j-number-1,0,v->f(a,j-number-1),1);
									found = true;
								}

								else if (mod[a]||mod[j-number-1]&&inc[ct->numseq[a]][ct->numseq[j-number-1]]) {

									if (energy==w5[a-1]+w3[i+1]+penalty(i,j,ct,data)+
											penalty(a,j-number-1,ct,data)+v->f(a+1,j-number-2)+
											ergcoaxflushbases(a,j-number-1,j-number,i,ct,data)
											+erg1(a,j-number-1,a+1,j-number-2,ct,data)) {

										registerbasepair(ct,a,j-number-1);
										if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
										if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
										stack->push(a+1,j-number-2,0,v->f(a+1,j-number-2),1);
										found = true;
									}
								}

								if (!found&&energy==w5[a-1]+w3[i+1]+penalty(i,j,ct,data)+
										penalty(a+1,j-number-2,ct,data)+v->f(a+1,j-number-2)+
										ergcoaxinterbases1(a+1,j-number-2,j-number,i,ct,data)
								   ) {

									if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
									if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
									stack->push(a+1,j-number-2,0,v->f(a+1,j-number-2),1);
									found = true;
								}

								else if (!found&&(mod[a+1]||mod[j-number-2])&&inc[ct->numseq[a+1]][ct->numseq[j-number-2]]) {

									if (energy==w5[a-1]+w3[i+1]+penalty(i,j,ct,data)+
											penalty(a+1,j-number-2,ct,data)+v->f(a+2,j-number-3)+
											ergcoaxinterbases1(a+1,j-number-2,j-number,i,ct,data)
											+erg1(a+1,j-number-2,a+2,j-number-3,ct,data)) {

										registerbasepair(ct,a+1,j-number-2);
										if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
										if (i+1<number-minloop-1) stack->push(i+1,number,1,w3[i+1],0);
										stack->push(a+2,j-number-3,0,v->f(a+2,j-number-3),1);
										found = true;
									}
								}

								if (!found&&energy==w5[a-1]+w3[i+2]+penalty(i,j,ct,data)+
										penalty(a,j-number-2,ct,data)+v->f(a,j-number-2)+
										ergcoaxinterbases2(a,j-number-2,j-number,i,ct,data)) {

									if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
									if (i+2<number-minloop-1) stack->push(i+2,number,1,w3[i+2],0);
									stack->push(a,j-number-2,0,v->f(a,j-number-2),1);
									found = true;
								}

								else if (!found&&(mod[a]||mod[j-number-2])&&inc[ct->numseq[a]][ct->numseq[j-number-2]]) {

									if (energy==w5[a-1]+w3[i+2]+penalty(i,j,ct,data)+
											penalty(a,j-number-2,ct,data)+v->f(a+1,j-number-3)+
											ergcoaxinterbases2(a,j-number-2,j-number,i,ct,data)
											+erg1(a,j-number-2,a+1,j-number-3,ct,data)) {

										registerbasepair(ct,a,j-number-2);
										if (a-1>minloop+1) stack->push(1,a-1,1,w5[a-1],0);
										if (i+2<number-minloop-1) stack->push(i+2,number,1,w3[i+2],0);
										stack->push(a+1,j-number-3,0,v->f(a+1,j-number-3),1);
										found = true;
									}
								}
								a--;
							}

							//now consider stacking to the 3' fragment:
							a = i+minloop+2;
							while (a<=number&&!found) {
								if (energy==w5[j-number-1]+w3[a+1]+penalty(i,j,ct,data)+
										penalty(i+1,a,ct,data) +v->f(i+1,a)+
										ergcoaxflushbases(j-number,i,i+1,a,ct,data)) {

									if (j-number-1>minloop+1) stack->push(1,j-number-1,1,w5[j-number-1],0);
									if (a+1<number-minloop-1) stack->push(a+1,number,1,w3[a+1],0);

									stack->push(i+1,a,0,v->f(i+1,a),1);
									found = true;
								}

								else if (mod[i+1]||mod[a]&&inc[ct->numseq[i+1]][ct->numseq[a]]) {
									if (energy==w5[j-number-1]+w3[a+1]+penalty(i,j,ct,data)+
											penalty(i+1,a,ct,data) +v->f(i+2,a-1)+
											ergcoaxflushbases(j-number,i,i+1,a,ct,data)
											+erg1(i+1,a,i+2,a-1,ct,data)) {

										registerbasepair(ct,i+1,a);
										if (j-number-1>minloop+1) stack->push(1,j-number-1,1,w5[j-number-1],0);
										if (a+1<number-minloop-1) stack->push(a+1,number,1,w3[a+1],0);

										stack->push(i+2,a-1,0,v->f(i+2,a-1),1);
										found = true;
									}
								}

								if (j-number-2>-1) {
									if (!found&&energy==w5[j-number-2]+w3[a+1]+penalty(i,j,ct,data)+
											penalty(i+2,a,ct,data)+v->f(i+2,a)+
											ergcoaxinterbases1(j-number,i,i+2,a,ct,data)
									   ) {

										if (j-number-2>minloop+1) stack->push(1,j-number-2,1,w5[j-number-2],0);
										if (a+1<number-minloop-1) stack->push(a+1,number,1,w3[a+1],0);

										stack->push(i+2,a,0,v->f(i+2,a),1);
										found = true;
									}
									else if (!found&&(mod[i+2]||mod[a])&&inc[ct->numseq[i+2]][ct->numseq[a]]) {
										if (!found&&energy==w5[j-number-2]+w3[a+1]+penalty(i,j,ct,data)+
												penalty(i+2,a,ct,data)+v->f(i+3,a-1)+
												ergcoaxinterbases1(j-number,i,i+2,a,ct,data)
												+erg1(i+2,a,i+3,a-1,ct,data)) {

											registerbasepair(ct,i+2,a);
											if (j-number-2>minloop+1) stack->push(1,j-number-2,1,w5[j-number-2],0);
											if (a+1<number-minloop-1) stack->push(a+1,number,1,w3[a+1],0);

											stack->push(i+3,a-1,0,v->f(i+3,a-1),1);
											found = true;
										}
									}
								}

								if (!found&&energy==w5[j-number-1]+w3[a+1]+penalty(i,j,ct,data)+
										penalty(i+2,a-1,ct,data)+v->f(i+2,a-1)+
										ergcoaxinterbases2(j-number,i,i+2,a-1,ct,data)) {

									if (j-number-1>minloop+1) stack->push(1,j-number-1,1,w5[j-number-1],0);
									if (a+1<number-minloop-1) stack->push(a+1,number,1,w3[a+1],0);

									stack->push(i+2,a-1,0,v->f(i+2,a-1),1);
									found = true;

								}

								else if (!found&&(mod[i+2]||mod[a-1])&&inc[ct->numseq[i+2]][ct->numseq[a-1]]) {
									if (!found&&energy==w5[j-number-1]+w3[a+1]+penalty(i,j,ct,data)+
											penalty(i+2,a-1,ct,data)+v->f(i+3,a-2)+
											ergcoaxinterbases2(j-number,i,i+2,a-1,ct,data)
											+erg1(i+2,a-1,i+3,a-2,ct,data)) {

										registerbasepair(ct,i+2,a-1);
										if (j-number-1>minloop+1) stack->push(1,j-number-1,1,w5[j-number-1],0);
										if (a+1<number-minloop-1) stack->push(a+1,number,1,w3[a+1],0);

										stack->push(i+3,a-2,0,v->f(i+3,a-2),1);
										found = true;

									}
								}
								a++;
							}
						}
						if (found) break;
						k++;
						}
					}

					if (ct->intermolecular) {


						if (energy==wmb2->f(i+1,j-1)+penalty(i,j,ct,data)+INFINITE_ENERGY) {

							stack->push(i+1,j-1,0,wmb2->f(i+1,j-1),2);
							found = true;

						}
						else if (energy==wmb2->f(i+2,j-1)+penalty(i,j,ct,data)+
								erg4(i,j,i+1,1,ct,data,lfce[i+1])+INFINITE_ENERGY) {

							stack->push(i+2,j-1,0,wmb2->f(i+2,j-1),2);
							found = true;

						}
						else if (energy==wmb2->f(i+1,j-2)+penalty(i,j,ct,data)+
								erg4(i,j,j-1,2,ct,data,lfce[j-1])+INFINITE_ENERGY) {

							stack->push(i+1,j-2,0,wmb2->f(i+1,j-2),2);
							found = true;

						}
						else if (energy==wmb2->f(i+2,j-2)+penalty(i,j,ct,data)+
								INFINITE_ENERGY+
								data->tstkm[ct->numseq[i]][ct->numseq[j]]
								[ct->numseq[i+1]][ct->numseq[j-1]]
								+checknp(lfce[i+1],lfce[j-1])) {

							stack->push(i+2,j-2,0,wmb2->f(i+2,j-2),2);
							found = true;

						}
					}
					if (!found) {
						//neither a multiloop nor an exterior loop were found,
						//there must be a bulge/interior loop


						for (d=j-i-minloop;(d>=1)&&!found;d--) {
							for (a=i+1;(a<=j-1-d)&&!found;a++) {
								b = d+a;

								//if (abs(a-i+b-j)<=data->eparam[8]) { 

								if (energy==(erg2(i,j,a,b,ct,data,fce->f(i,a),fce->f(b,j))+
											v->f(a,b))) {
									i = a;
									j = b;
									stack->push(i,j,0,v->f(i,j),1);
									found = true;
								}
								//}
								if (mod[a]||mod[b]&&inc[ct->numseq[a]][ct->numseq[b]]) {
									if (energy==(erg2(i,j,a,b,ct,data,fce->f(i,a),fce->f(b,j))+
												v->f(a+1,b-1)+erg1(a,b,a+1,b-1,ct,data))) {
										i = a+1;
										j = b-1;
										registerbasepair(ct,a,b);
										stack->push(i,j,0,v->f(i,j),1);
										found = true;
									}
								}
							}
						}
					}

					if (!found) {//something went wrong					
						errmsg (100,2);
						tracebackerror=1;
					}
				}
			}
		}

		delete stack;
		return tracebackerror;
	}



#define maxsort 90000 //Starting point for the maximum number of basepairs within %cntrl8
//of the minimum free energy.


//Traceback predicts a set of low free energy structures using the mfold heuristic with the fill step information.
//cntrl8 is the maximum % difference in free energy of suboptimal structures if > 0
//	otherwise, cntrl8 is a maximum energy difference in kcal/mol*factor

//This returns an error code, where zero is no error and non-zero indicates a traceback error.
int traceback(structure *ct, datatable *data, arrayclass *v, arrayclass *w, arrayclass *wmb, arrayclass *w2,arrayclass *wmb2, integersize *w3, integersize *w5, forceclass *fce,
              bool *lfce,integersize vmin, int cntrl6, int cntrl8, int cntrl9, bool *mod) {

	bool flag,**mark;
	register int number;
	int ii,sort;
	int i;
	int iret,jret,numbp,count,count2,k1,k2,num,crit;
	int *heapi,*heapj;
	int cur,c,k,j,cntr;
	int ji;//(added during debugging)
	int tracebackerror=0;

	integersize *energy;

#if defined(debugmode)
	char filename[maxfil];
	char temp[20];
#endif

	//mark keeps track of which pairs have been formed by the suboptimal routine

	sort = maxsort;
	number = ct->GetSequenceLength();

	//Construct heapi and heapj composed of pairs ij such that a structure
	//	containing ij has energy less than a given percent (ie:cntrl8) of
	//	the minimum folding energy

	if (vmin>=0) { //no viable structure found
		//Add an empty structure to ct:
		ct->AddStructure();
		//Set the energy to zero, a random coil:
		ct->SetEnergy(1,0);
		return 0;
	}

	//dynamically allocate space for mark:
	mark = new bool *[number + 1];
	for (i=0;i<=(number);i++)
		mark[i] = new bool [number + 1];

	//This is the traceback portion of the dynamic algorithm

	flag = true;
	

	for (count=1;count<=(number);count++) {
		for (count2=1;count2<=(number);count2++) {
			mark[count][count2]=false;
		}
	}

	if (cntrl8> 0) crit= (short) (abs(vmin)*(float (cntrl8)/100.0));
	else crit = -cntrl8;

	crit = crit + vmin;

	energy = new integersize [sort+1];
	heapi = new int [sort+1];
	heapj = new int [sort+1];

	num = 0;
	i = 1;
	j = 2;
	while (i<(number)) {
		if (num==sort) {
			//allocate more space for the heap
			delete[] heapi;
			delete[] heapj;
			delete[] energy;
			sort = 10*sort;
			heapi = new int [sort+1];
			heapj = new int [sort+1];
			energy = new integersize [sort+1];
			i = 1;
			j = 2;
			num = 0;

		}

		if ((v->f(i,j)+v->f(j,i+number))<=crit) {

			num++;
			heapi[num]=i;
			heapj[num]=j;
			energy[num] = (v->f(i,j)+v->f(j,i+number));
			j = j+cntrl9+1;
			if (j>number) {
				i++;
				j=i+1;
			}
		}
		else if (mod[i]||mod[j]) {
			//add i-j pair to heap if it is adjacent to unpaired nucs in one direction
			if (v->f(i,j)<INFINITE_ENERGY) {
				cur = v->f(i,j)+v->f(j+1,i+number-1)+erg1(j,i+number,j+1,i+number-1,ct,data);
				if (cur<=crit) {
					num++;
					heapi[num]=i;
					heapj[num]=j;
					energy[num] = crit;
					j = j+cntrl9+1;
					if (j>number) {
						i++;
						j=i+1;
					}

				}
				else {
					j++;
					if (j>number) {
						i++;
						j=i+1;
					}

				}

			}
			else if (v->f(j,i+number)==INFINITE_ENERGY) {
				cur = v->f(i+1,j-1)+v->f(j,i+number)+erg1(i,j,i+1,j-1,ct,data);
				if (cur<=crit) {
					num++;
					heapi[num]=i;
					heapj[num]=j;
					energy[num] = crit;
					j = j+cntrl9+1;
					if (j>number) {
						i++;
						j=i+1;
					}

				}
				else {
					j++;
					if (j>number) {
						i++;
						j=i+1;
					}

				}

			}
			else {
				j++;
				if (j>number) {
					i++;
					j=i+1;
				}

			}

		}
		else {
			j++;
			if (j>number) {
				i++;
				j=i+1;
			}
		}
	}

	//sort the base pair list:
	///////////////////////////////////

	//make a heap:

	int q,up,ir;
	energy[0]=0;//energy[0] is never used, but gets tested in the case that up=0;
	for (q=2;q<=num;q++) {
		cur = q;
		up = cur/2;
		while ((energy[cur]<energy[up])&&up>=1) {
			swap(&heapi[cur],&heapi[up]);
			swap(&heapj[cur],&heapj[up]);
			swap(&energy[cur],&energy[up]);
			cur = cur/2;
			up = up/2;
		}
	}

	//sort the heap:

	for (ir=num-1;ir>=1;ir--) {
		swap(&heapi[ir+1],&heapi[1]);
		swap(&heapj[ir+1],&heapj[1]);
		swap(&energy[ir+1],&energy[1]);

		up =1 ;
		c = 2;
		while (c<=ir) {
			if (c!=ir) {
				if (energy[c+1]<energy[c]) c++;
			}
			if (energy[c]<energy[up]) {
				swap(&heapi[c],&heapi[up]);
				swap(&heapj[c],&heapj[up]);
				swap(&energy[c],&energy[up]);
				up=c;
				c=2*c;
			}
			else c = ir+1;
		}
	}

	cntr = num;

	//keep track if a structure was allocated, but has been rejected with failedprevious.
	bool failedprevious = false;
	while (flag) {
		//This is routine to select the region of the structure to be
		//	folded, it allows for sub-optimal structure predictions
		//err=0;
		//Select the next valid unmarked basepair
		while (mark[heapi[cntr]][heapj[cntr]]) {
			if (cntr==1) {
				flag=false;
				goto sub900;
			}
			cntr--;
		}
		iret = heapi[cntr];
		jret = heapj[cntr];
		//No prior structure was allocated, so add a structure now:
		if (!failedprevious) ct->AddStructure();
		//The prior failed structure needs to be cleaned, i.e. all base pairs removed:
		else ct->CleanStructure(ct->GetNumberofStructures());
		//Reset failedprevious as a new traceback is started:
		failedprevious=false;

		//Traceback to find best structure on included fragment (ie:iret to jret)
		if (flag) {
			for (count=1;count<=2;count++) {
				if (count==1) {
					ii=iret;
					ji=jret;
			
				}
				if (count==2) {
					ii=jret;
					ji=iret+(number);
				}

				if (mod[iret]||mod[jret]) {

					if (v->f(ii,ji)==INFINITE_ENERGY) {
				
						ii+=1;
						ji-=1;
					}
				}
#ifndef DYNALIGN_II
				if(trace(ct,data,ii,ji,v,w,wmb,w2,wmb2,lfce,fce,w3,w5,mod)!=0) tracebackerror=1;
#else
				if(trace(ct,data,ii,ji,v,w,wmb,w2,wmb2,lfce,fce,w3,w5,mod,NULL)!=0) tracebackerror=1;
#endif
				if (count==2) {

					ct->SetEnergy(ct->GetNumberofStructures(),energy[cntr]);
					//count the number of new base pairs not within window of existing
					//base pairs
					numbp = 0;
					for (k=1;k<=number;k++) {
						if (k<(ct->GetPair(k,ct->GetNumberofStructures()))) {
							if (!(mark[k][ct->GetPair(k,ct->GetNumberofStructures())])) numbp++;
						}
					}
					for (k=1;k<=(number);k++) {
						if (k<ct->GetPair(k,ct->GetNumberofStructures())) {
							//Mark "traced back" base pairs and also base pairs
							//	which are within a window of cntrl9
							mark[k][ct->GetPair(k,ct->GetNumberofStructures())]=true;
							if (cntrl9>0) {
								for (k1=max(1,k-cntrl9);k1<=min(number,k+cntrl9);k1++) {
									for (k2=max(k1,ct->GetPair(k,ct->GetNumberofStructures())-cntrl9);k2<=min(number,ct->GetPair(k,ct->GetNumberofStructures())+cntrl9);k2++) {
                        	
										mark[k1][k2]=true;
									}
								}
							}
						}
					}
					if (numbp<=cntrl9&&ct->GetNumberofStructures()>1) {
						failedprevious=true;

						//goto sub900;
					}
					else {

						//place the structure name (from ctlabel[1]) into each structure
						ct->SetCtLabel(ct->GetSequenceLabel(), ct->GetNumberofStructures());
				
#if defined(debugmode)
					
						strcpy(filename,"energydump");
						itoa(ct->GetNumberofStructures(),temp,10);
						strcat(filename,temp);
						strcat(filename,".out");
						energydump (ct, v, data, ct->GetNumberofStructures(),filename,iret,jret);
#endif

					}
					if (ct->GetNumberofStructures()==cntrl6&&!failedprevious) {
						flag=false;
					}
				}
			}
		sub900:
			continue;
		}
	}
	if (failedprevious) ct->RemoveLastStructure();

	de_allocate (mark,number+1);
	delete[] energy;
	delete[] heapi;
	delete[] heapj;
	return tracebackerror;
}


//force is used to enforce the folding constraints specified in a structure
//The following definitions are bitwise applied to the fce array:
//SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
//PAIR applies to any fce(i,j) where i is paired to j
//NOPAIR applies to any fce(i,j) where either i or j is paired to
//	another nucleotide or i and j are forbidden to pair
//DOUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
//INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
//	used for intermolecular folding
//INTER applies to any fce(i,i) s.t. i is the center of the virtual linker
//The above terms are defined in define.h


//mod[i] is a bool array that is set to true if i is a nuc accessible to chemical
//	modification
//lfce[i] is a bool array that is set to true if i is double stranded


//force also double the sequence s.t. ct->numseq[i+N] = ct->numseq[i] where
//	N is the number of nucleotides


void force(structure *ct,forceclass *fce, bool *lfce){
	int i,j;
	register int number;

	number = ct->GetSequenceLength();

	for (i=0;i<ct->GetNumberofSingles();i++) {
		
		if(ct->GetSingle(i)<=ct->GetSequenceLength())	forcesingle(ct->GetSingle(i),ct,fce);
	}

	for (i=0;i<ct->GetNumberofPairs();i++) {
		if(ct->GetPair5(i)<=ct->GetSequenceLength()&&ct->GetPair3(i)<=ct->GetSequenceLength()) {
			forcepair(ct->GetPair5(i),ct->GetPair3(i),ct,fce);
			forcedbl(ct->GetPair5(i),ct,fce,lfce);
			forcedbl(ct->GetPair3(i),ct,fce,lfce);
		}
	}

	for (i=0;i<ct->GetNumberofDoubles();i++) {
		if (ct->GetDouble(i)<=ct->GetSequenceLength()) forcedbl(ct->GetDouble(i),ct,fce,lfce);
	}


	//u's in gu pairs must be double stranded
	for (i=0;i<ct->GetNumberofGU();i++) {
		if (ct->GetGUpair(i)<=ct->GetSequenceLength()) forcedbl(ct->GetGUpair(i),ct,fce,lfce);
	}



	if (ct->intermolecular) {//this indicates an intermolecular folding
		//for (i=0;i<=3;i++) {
		// 	forcesingle(ct->inter[i])//don't allow the intermolecular indicators to pair
		//}
		for (i=0;i<3;i++) {
			forceinter(ct->inter[i],ct,fce);
		}

		fce->f(ct->inter[1],ct->inter[1]) = fce->f(ct->inter[1],ct->inter[1])|INTER;
	}

	for (i=0;i<ct->GetNumberofForbiddenPairs();i++) {
		if(ct->GetForbiddenPair5(i)<=ct->GetSequenceLength()&ct->GetForbiddenPair3(i)<=ct->GetSequenceLength()) fce->f(ct->GetForbiddenPair5(i),ct->GetForbiddenPair3(i)) = fce->f(ct->GetForbiddenPair5(i),ct->GetForbiddenPair3(i))|NOPAIR;
		if(ct->GetForbiddenPair5(i)<=ct->GetSequenceLength()&&ct->GetForbiddenPair3(i)<=ct->GetSequenceLength()) fce->f(ct->GetForbiddenPair3(i),ct->GetForbiddenPair5(i)+ct->GetSequenceLength())=fce->f(ct->GetForbiddenPair3(i),ct->GetForbiddenPair5(i)+ct->GetSequenceLength())|NOPAIR;
	}

	//Double up the sequence
	for (i=1;i<=number;i++) {
		ct->numseq[(number)+i] = ct->numseq[i];
	}

	//The next section handles the case where base pairs are not
	//not allowed to form between nucs more distant
	//than ct->GetPairingDistanceLimit()
	if (ct->DistanceLimited()) {

		if (!ct->templated) ct->allocatetem();

		for (j=minloop+2;j<=ct->GetSequenceLength();j++) {
			for (i=1;i<j;i++) {
				if (j-i>=ct->GetPairingDistanceLimit()) ct->tem[j][i]=false;
			}
		}
	}
}
#else
#endif

#ifndef DYNALIGN_II
	//This is the dynamic algorithm of Zuker:
	//cntrl6 = #tracebacks
	//cntrl8 = percent sort
	//cntrl9 = window

	//TProgressDialog is an interface for returning the progress of the calculation
	//Savfile is for creating a file with arrays and parameters for refolding with different 
	//suboptimal tracebacks
	//quickenergy indicates whether to find the lowest free energy for the sequence without a structure
#ifndef INSTRUMENTED
	int dynamic(structure* ct,datatable* data,int cntrl6, int cntrl8,int cntrl9,
			TProgressDialog* update, bool quickenergy, char* save, int maxinter, bool quickstructure)


#else //INSTRUMENTED IS DEFINED
		void dynamic(structure* ct,datatable* data,int cntrl6, int cntrl8,int cntrl9,
				arrayclass *v, arrayclass *vmb/*tracks MB loops*/, arrayclass *vext/*tracks exterior loops*/,
				TProgressDialog* update, bool quickenergy, char* save, int maxinter, bool quickstructure)
#endif //END of INSTRUMENTED iS DEFINED
		{		
			int number;		
			int vmin;
			int i,j,k,l;
			int m,n,o,p;

			bool *lfce,*mod;//[maxbases+1][maxbases+1];

			//integersize **work,**work2;
			integersize *w5,*w3;//**wca,**curE,**prevE,**tempE;

			int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
				{0,1,0,1,0,0},{0,0,0,0,0,0}};

			int tracebackerror=0;


			//array *v,*w;
			//inc is an array that saves time by showing which bases can pair before
			//	erg is called

			//number is the number of bases being folded
			//v[i][j] is the best energy for subsequence i to j when i and j are paired
			//	for i<j<n, it is the interior framgment between nucleotides i and j
			//	for i<n<j, it is the exterior structure fragmnet that is 5' to j-n and 3' to i
			//w[i][j] is the best energy for subsequence i to j

//<<<<<<< algorithm.cpp
	number = (ct->GetSequenceLength());//place the number of bases in a registered integer
//=======
//			number = (ct->numofbases);//place the number of bases in a registered integer
//>>>>>>> 1.40

#ifdef timer //time the algorithm execution
#include <time.h>
//<<<<<<< algorithm.cpp
	ofstream timeout;
	int seconds;
	char timerstring[100];
	char timelength[10];
	strcpy(timerstring,"time_pred_");
	sprintf(timelength,"%i",ct->GetSequenceLength());
	strcat(timerstring,timelength);
	strcat(timerstring,".out");

	timeout.open(timerstring);
	timeout<<time(NULL)<<"\n";
	seconds = time(NULL);
//=======
//			ofstream timeout;
//			int seconds;
//			char timerstring[100];
//			char timelength[10];
//			strcpy(timerstring,"time_pred_");
//			sprintf(timelength,"%i",ct->numofbases);
//			strcat(timerstring,timelength);
//			strcat(timerstring,".out");

//			timeout.open(timerstring);
//			timeout<<time(NULL)<<"\n";
//			seconds = time(NULL);
//>>>>>>> 1.40
#endif

			//allocate space for the v and w arrays:
			arrayclass w(number);
#ifndef INSTRUMENTED // if INSTRUMENTED compiler flag is defined
			arrayclass v(number);
#endif
			arrayclass wmb(number);
			forceclass fce(number);

			//add a second array for intermolecular folding:
			arrayclass *w2,*wmb2;
			if (ct->intermolecular) {
				w2 = new arrayclass(number);
				wmb2 = new arrayclass(number);
				//	work2 = new integersize *[2*number+3];

				//	for (i=0;i<2*number+3;i++) {
				//		work2[i] = new integersize [3];
				//		work2[i][0]=INFINITE_ENERGY;
				//		work2[i][1]=INFINITE_ENERGY;
				//		work2[i][2]=INFINITE_ENERGY;

				//	}

			}
			else {

				wmb2=NULL;
				w2=NULL;
				//	work2=NULL;
			}

			lfce = new bool [2*number+1];
			mod = new bool [2*number+1];

			for (i=0;i<=2*number;i++) {
				lfce[i] = false;
				mod[i] = false;
			}

			//for (i=1;i<=ct->nmod;i++) {

			//	if (ct->mod[i]>1&&ct->mod[i]<ct->numofbases) {
			//		mod[ct->mod[i]]=true;
			//		mod[ct->mod[i]+ct->numofbases]=true;
			//	}
			//}

			//	work = new integersize *[2*number+3];

//<<<<<<< algorithm.cpp
	for (i=0;i<ct->GetNumberofModified();i++) {
//=======
			//	for (i=0;i<2*number+3;i++) {
			//   		work[i] = new integersize [3];
			//		work[i][0]=INFINITE_ENERGY;
			//		work[i][1]=INFINITE_ENERGY;
			//		work[i][2]=INFINITE_ENERGY;
//>>>>>>> 1.40

//<<<<<<< algorithm.cpp
		if (ct->GetModified(i)>1&&ct->GetModified(i)<ct->GetSequenceLength()) {
			mod[ct->GetModified(i)]=true;
			mod[ct->GetModified(i)+ct->GetSequenceLength()]=true;
		}
	}
//=======
			//	}
			w5 = new integersize [number+1];
			w3 = new integersize [number+2];
//>>>>>>> 1.40

			for (i=0;i<=number;i++) {
				w5[i] = 0;
				w3[i] = 0;
			}

			w3[number+1] = 0;
			force(ct,&fce,lfce);
			vmin=INFINITE_ENERGY;

#ifndef INSTRUMENTED//If pre-compiler flag INSTRUMENTED is not defined, compile the following code
#ifndef DYNALIGN_II
			//perform the fill steps:(i.e. fill arrays v and w.)
			fill(ct, v, w, wmb, fce, vmin,lfce, mod,w5, w3, quickenergy, data, w2, wmb2, update, maxinter,quickstructure);
#else
                        fill(ct, v, w, wmb, fce, vmin,lfce, mod,w5, w3, quickenergy, data, w2, wmb2, NULL, update, maxinter);
#endif

//<<<<<<< algorithm.cpp
	if (save!=0) {
		ofstream sav(save,ios::binary);
	
		//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors
		short vers=safiversion;
		write(&sav,&vers); //save a version of the save file 
		//start with structure information
		int sequencelength=ct->GetSequenceLength();
		write(&sav,&(sequencelength));
		write(&sav,&(ct->intermolecular));

		int pairs = ct->GetNumberofPairs();
		write(&sav,&(pairs));
		for (i=0;i<ct->GetNumberofPairs();i++) {
			pairs=ct->GetPair5(i);
			write(&sav,&(pairs));
			pairs=ct->GetPair3(i);
			write(&sav,&(pairs));
		}

		pairs = ct->GetNumberofForbiddenPairs();
		write(&sav,&(pairs));
		for (i=0;i<ct->GetNumberofForbiddenPairs();i++) {
			pairs = ct->GetForbiddenPair5(i);
			write(&sav,&(pairs));
			pairs = ct->GetForbiddenPair3(i);
			write(&sav,&(pairs));
		}
		for (i=0;i<=ct->GetSequenceLength();i++) {
		
			write(&sav,&(ct->hnumber[i]));
			sav.write(&(ct->nucs[i]),1);

		}

		for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->numseq[i]));
	
		int doubles = ct->GetNumberofDoubles();
		write(&sav,&(doubles));
		for (i=0;i<ct->GetNumberofDoubles();i++) {
			doubles = ct->GetDouble(i);
			write(&sav,&(doubles));
		}

	
		if (ct->intermolecular) {
			for (i=0;i<3;i++) write(&sav,&(ct->inter[i]));
		
		}

		int singles = ct->GetNumberofSingles();
		write(&sav,&(singles));
		for (i=0;i<ct->GetNumberofSingles();i++) {
			singles = ct->GetSingle(i);	
			write(&sav,&(singles));
		}

		int modified;
		modified = ct->GetNumberofModified();
		write(&sav,&(modified));
		for (i=0;i<ct->GetNumberofModified();i++) {
			modified = ct->GetModified(i);
			write(&sav,&(modified));
		}
	
		modified = ct->GetNumberofGU();
		write(&sav,&(modified));
		for (i=0;i<ct->GetNumberofGU();i++) {
			modified = ct->GetGUpair(i);
			write(&sav,&(modified));
		}
	
		string label=ct->GetSequenceLabel();
		write(&sav,&(label));

		write(&sav,&(ct->templated));
		if (ct->templated) {
			for (i=0;i<=ct->GetSequenceLength();i++) {
				for (j=0;j<=i;j++) write(&sav,&(ct->tem[i][j]));	

			}

		}

		//write the SHAPE data (for pseudo-free energy constraints)
		write(&sav,&(ct->shaped));
		if (ct->shaped) {
			for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->SHAPE[i]));

		}

	
		//now write the array class data for v, w, and wmb:
		for (i=0;i<=ct->GetSequenceLength();i++) {
			write(&sav,&(w3[i]));
			write(&sav,&(w5[i]));
			for (j=0;j<=ct->GetSequenceLength();j++) {
				write(&sav,&(v.dg[i][j+i]));
				write(&sav,&(w.dg[i][j+i]));
				write(&sav,&(wmb.dg[i][j+i]));
				writesinglechar(&sav,&(fce.dg[i][j]));


				if (ct->intermolecular) {
					write(&sav,&(w2->dg[i][j+i]));
					write(&sav,&(wmb2->dg[i][j+i]));
				}
			}	
		}

		write(&sav,&(w3[ct->GetSequenceLength()+1]));
		for (i=0;i<=2*ct->GetSequenceLength();i++) {
			write(&sav,&(lfce[i]));
			write(&sav,&(mod[i]));
		}

		write(&sav,&vmin);

		//now write the thermodynamic data:
		for (i=0;i<5;i++) write(&sav,&(data->poppen[i]));
		write(&sav,&(data->maxpen));
		for (i=0;i<11;i++) write(&sav,&(data->eparam[i]));
		for (i=0;i<31;i++) {
			write(&sav,&(data->inter[i]));
			write(&sav,&(data->bulge[i]));
			write(&sav,&(data->hairpin[i]));

		}
		for (i=0;i<6;i++) {
			for (j=0;j<6;j++) {
				for (k=0;k<6;k++) {
					for (l=0;l<3;l++) {
						write(&sav,&(data->dangle[i][j][k][l]));	
					}
					for (l=0;l<6;l++) {
						write(&sav,&(data->stack[i][j][k][l]));
						write(&sav,&(data->tstkh[i][j][k][l]));
						write(&sav,&(data->tstki[i][j][k][l]));
						write(&sav,&(data->coax[i][j][k][l]));
						write(&sav,&(data->tstackcoax[i][j][k][l]));
						write(&sav,&(data->coaxstack[i][j][k][l]));
						write(&sav,&(data->tstack[i][j][k][l]));
						write(&sav,&(data->tstkm[i][j][k][l]));
						write(&sav,&(data->tstki23[i][j][k][l]));
						write(&sav,&(data->tstki1n[i][j][k][l]));
						for (m=0;m<6;m++) {
							for (n=0;n<6;n++) {
								write(&sav,&(data->iloop11[i][j][k][l][m][n]));
								for (o=0;o<6;o++) {
									if (inc[i][j]&&inc[n][o]) write(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
									for (p=0;p<6;p++) {
										if (inc[i][k]&&inc[j][l])
											write(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
									}
								}
							}
						}
					}
				}
			}
		}
		write(&sav,&(data->numoftloops));
		for (i=0;i<=data->numoftloops;i++) {
			for (j=0;j<2;j++) write(&sav,&(data->tloop[i][j]));

		}
		write(&sav,&(data->numoftriloops));
		for (i=0;i<=data->numoftriloops;i++) {
			for (j=0;j<2;j++) write(&sav,&(data->triloop[i][j]));

		}
		write(&sav,&(data->numofhexaloops));
		for (i=0;i<=data->numofhexaloops;i++) {
			for (j=0;j<2;j++) write(&sav,&(data->hexaloop[i][j]));

		}
		write(&sav,&(data->auend));
		write(&sav,&(data->gubonus));
		write(&sav,&(data->cint));
		write(&sav,&(data->cslope));
		write(&sav,&(data->c3));
		write(&sav,&(data->efn2a));
		write(&sav,&(data->efn2b));
		write(&sav,&(data->efn2c));
		write(&sav,&(data->init));
		write(&sav,&(data->mlasym));
		write(&sav,&(data->strain));
		write(&sav,&(data->prelog));
		write(&sav,&(data->singlecbulge));
	
		sav.close();
	}

	if (quickenergy) {
		//Don't do traceback, just return energy

		ct->SetEnergy(1,w5[number]);
		tracebackerror=0;
	}
	else if (quickstructure) {
		//Calculate only the lowest free energy structure
		tracebackerror=trace(ct,data,1,ct->GetSequenceLength(),&v,&w,&wmb,w2,wmb2,lfce,&fce,w3,w5,mod,true);

	}

	else tracebackerror=traceback(ct, data, &v, &w, &wmb, w2, wmb2,w3, w5, &fce, lfce, vmin, cntrl6, cntrl8, cntrl9,mod);

//<<<<<<< algorithm.cpp
//	if (quickenergy) {
		//Don't do traceback, just return energy
		
//		ct->SetEnergy(1,w5[number]);
//	}
	
//	else tracebackerror=traceback(ct, data, &v, &w, &wmb, w2, wmb2,w3, w5, &fce, lfce, vmin, cntrl6, cntrl8, cntrl9,mod);
	
//======
//>>>>>>> 1.40
#else//INSTRUMENTED IS DEFINED

//<<<<<<< algorithm.cpp
	fill(ct, *v, *vmb, *vext, w, wmb, fce, vmin,lfce, mod,w5, w3, quickenergy, data, w2, wmb2, update, maxinter);
	/*
	//write a dot plot file that tracks MB loops and exterior loops
	ofstream sav(save);
	double pairs;
	double mbpairs;
	double extpairs;
		
	for (i=1;i<=ct->GetSequenceLength();i++) {
		for (j=(i+1);j<=ct->GetSequenceLength();j++) {
				
			if ((v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))<0) {
					
				pairs= ((double)(v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0);
				mbpairs=min( ((double)(vmb->f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0), ((double)(v.f(i,j)+vmb->f(j,i+ct->GetSequenceLength()))/10.0));
				extpairs= ((double)(v.f(i,j)+vext->f(j,i+ct->GetSequenceLength()))/10.0);
				sav << i << "\t" << j << "\t" << pairs<<"\t"<< mbpairs<< "\t"<<extpairs<<"\n"; 
					
					
//=======
//			fill(ct, *v, *vmb, *vext, w, wmb, fce, vmin,lfce, mod,w5, w3, quickenergy, data, w2, wmb2, update, maxinter);
			/*
			//write a dot plot file that tracks MB loops and exterior loops
//			ofstream sav(save);
//			double pairs;
//			double mbpairs;
//			double extpairs;

//			for (i=1;i<=ct->GetSequenceLength();i++) {
//			for (j=(i+1);j<=ct->GetSequenceLength();j++) {

//			if ((v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))<0) {

//			pairs= ((double)(v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0);
//			mbpairs=min( ((double)(vmb->f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0), ((double)(v.f(i,j)+vmb->f(j,i+ct->GetSequenceLength()))/10.0));
//			extpairs= ((double)(v.f(i,j)+vext->f(j,i+ct->GetSequenceLength()))/10.0);
//			sav << i << "\t" << j << "\t" << pairs<<"\t"<< mbpairs<< "\t"<<extpairs<<"\n"; 


//>>>>>>> 1.40
			}
			}
			}

			sav.close();
			 */


#endif//END INSTRUMENTED 

			delete[] lfce;
			delete[] mod;

			delete[] w5;
			delete[] w3;


			if (ct->intermolecular) {
				delete w2;
				delete wmb2;
			}

#ifdef timer
			timeout << time(NULL)<<"\n";
			timeout << time(NULL) - seconds;
			timeout.close();
#endif

#ifndef INSTRUMENTED // if INSTRUMENTED is not defined
			return tracebackerror;
#endif // if INSTRUMENTED is defined


		}

#else
#endif

	//The fill routine is encapsulated in function fill.
	//This was separated from dynamic on 3/12/06 by DHM.  This provides greater flexibility
	//for use of the arrays for other tasks than secondary structure prediction, e.g. dot plots.
#if defined DYNALIGN_II
void fill(structure *ct, arrayclass &v, arrayclass &w, arrayclass &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
          integersize *w5, integersize *w3, bool quickenergy,
          datatable *data, arrayclass *w2, arrayclass *wmb2, arrayclass *we,TProgressDialog* update, int maxinter, bool quickstructure)

#elif !defined INSTRUMENTED//If pre-compiler flag INSTRUMENTED is not defined, compile the following code
	void fill(structure *ct, arrayclass &v, arrayclass &w, arrayclass &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
			integersize *w5, integersize *w3, bool quickenergy,
			datatable *data, arrayclass *w2, arrayclass *wmb2, TProgressDialog* update, int maxinter,bool quickstructure)

#else //IF DEFINED INSTRUMENTED
		void fill(structure *ct, arrayclass &v, arrayclass &vmb, arrayclass &vext, arrayclass &w, arrayclass &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
				integersize *w5, integersize *w3, bool quickenergy,
				datatable *data, arrayclass *w2, arrayclass *wmb2, TProgressDialog* update, int maxinter,bool quickstructure)

#endif //end !INTRUMENTED
		{
			//This is the fill routine:
			register int number, h, maximum;
			int d, maxj;
			integersize **wca,**curE,**prevE,**tempE;


			number = ct->GetSequenceLength();

			if (!ct->intermolecular) {
				//This code is needed for O(N^3) prediction of internal loops
				wca = new integersize *[number+1];
				curE= new integersize *[number+1];
				prevE= new integersize *[number+1];

				for (int locali=0;locali<=number;locali++) {
					wca[locali] = new integersize [number+1];
					curE[locali]= new integersize [number+1];
					prevE[locali]= new integersize [number+1];
					for (int localj=0;localj<=number;localj++){
						wca[locali][localj] = INFINITE_ENERGY;
						curE[locali][localj]=INFINITE_ENERGY;
						prevE[locali][localj]=INFINITE_ENERGY;
					}
				}
			}
			else {//intermolecular folding
				//initialize wca only
				wca = new integersize *[number+1];
				for (int locali=0;locali<=number;locali++) {
					wca[locali] = new integersize [number+1];

					for (int localj=0;localj<=number;localj++)	{
						wca[locali][localj] = INFINITE_ENERGY;
					}
				}
			}

			if (quickenergy||quickstructure) maximum = number;
			else maximum = (2*(number)-1);


			vmin=INFINITE_ENERGY;
			//h or (h-number+1) is the distance between i and j
			for (h=0;h<=maximum-1;h++){

				//d = j-i;
				d=(h<=(number-1))?h:(h-number+1);
				if (((h%10)==0)&&update) update->update((100*h)/(maximum+1));
				if (h==number&&!ct->intermolecular) {
					for(int locali=0;locali<=number;locali++) {
						for(int localj=0;localj<=number;localj++) {
							curE[locali][localj]=INFINITE_ENERGY;
							prevE[locali][localj]=INFINITE_ENERGY;
						}
					}
				}
				

				//These variables for start and end (as opposed the fancy syntax that was here before) are needed because of openmp.
				int startme,endme;//start and end for loop over i
				if (h<=(number-1)) {
					startme=1;
					endme=number-h;
				}
				else {
					startme=2*number-h;
					endme=number;

				}
#ifdef SMP
#pragma omp parallel for
#endif
				
				

				for (int i=startme;i<=endme;i++){
					int j=i+d;
					register int rarray;
					int dp,ll,jpf,jf,bl,maxasym;
					int before,after;
					int e[6];
					int k,p;
					int ip,jp,ii,jj,di;
					register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
						{0,1,0,1,0,0},{0,0,0,0,0,0}};
#ifndef disablecoax
					register int castack; //variable for coaxial stacking, not needed if coaxial stacking is disabled
#endif
					int maxi=((h<=(number-1))?(number-h):number);
					if (i==maxi)
						maxj=j;


					maxasym = maxinter;

					if (ct->templated) {
						if (i>ct->GetSequenceLength()) ii = i - ct->GetSequenceLength();
						else ii = i;
						if (j>ct->GetSequenceLength()) jj = j - ct->GetSequenceLength();
						else jj = j;
						if (jj<ii) {
							p = jj;
							jj = ii;
							ii = p;
						}
						if (!ct->tem[jj][ii]) goto sub2;
					}

					//Compute v[i][j], the minimum energy of the substructure from i to j,
					//inclusive, where i and j are base paired
					if (fce.f(i,j)&SINGLE) {
						//i or j is forced single-stranded
						v.f(i,j) = INFINITE_ENERGY + 50;
						goto sub2;
					}
					if (fce.f(i,j)&NOPAIR) {
						//i or j is forced into a pair elsewhere
						v.f(i,j)= INFINITE_ENERGY+50;

						goto sub2;
					}

					if (j<=(number)) {
						if ((j-i)<=minloop) goto sub3;

					}
					v.f(i,j) = INFINITE_ENERGY;

					if (inc[ct->numseq[i]][ct->numseq[j]]==0) goto sub2;

					//force u's into gu pairs
					for (ip=0;ip<ct->GetNumberofGU();ip++) {
						if (ct->GetGUpair(ip)==i) {
							if (ct->numseq[j]!=3) {
								v.f(i,j) = INFINITE_ENERGY;
								goto sub2;
							}
						}
						//else if ((ct->gu[ip]+number)==i) {
						// 	if (ct->numseq[j]!=3) {
						//   	v.f(i,j) = INFINITE_ENERGY;
						//      goto sub2;
						//   }
						//}
						else if (ct->GetGUpair(ip)==j) {
							if (ct->numseq[i]!=3) {
								v.f(i,j) = INFINITE_ENERGY;
								goto sub2;
							}
						}
						else if ((ct->GetGUpair(ip)+number)==j) {
							if (ct->numseq[i]!=3) {
								v.f(i,j) = INFINITE_ENERGY;
								goto sub2;
							}
						}

					}


					//now check to make sure that this isn't an isolated pair:
					//	(consider a pair separated by a bulge as not! stacked)

					//before = 0 if a stacked pair cannot form 5' to i
					before =0;
					if ((i>1&&j<(2*number)&&j!=number)) {
						if ((j>number&&((i-j+number)>minloop+2))||j<number) {
							before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
						}
					}

					//after = 0 if a stacked pair cannot form 3' to i
					if ((((j-i)>minloop+2)&&(j<=number)||(j>number+1))&&(i!=number)) {
						after = inc[ct->numseq[i+1]][ct->numseq[j-1]];

					}
					else after = 0;

					//if there are no stackable pairs to i.j then don't allow a pair i,j
					if ((before==0)&&(after==0)) {
						//v.f(i,j)= 0;
						goto sub2;
					}
					rarray = INFINITE_ENERGY;
					if (i==(number)||j==((number)+1)) goto sub1;


					//Perhaps i and j close a hairpin:
					rarray=min(rarray,erg3(i,j,ct,data,fce.f(i,j)));

					if ((j-i-1)>=(minloop+2)||j>(number)) {
						//Perhaps i,j stacks over i+1,j-1
						if (!mod[i]&&!mod[j])  //make sure this is not a site of chemical modification
							rarray=min(rarray,(erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1)));
						else {
							//allow G-U to be modified or a pair next to a G-U to be modified
							if ((ct->numseq[i]==3&&ct->numseq[j]==4)||(ct->numseq[i]==4&&ct->numseq[j]==3)) {
								rarray=min(rarray,(erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1)));	

							}
							else if ((ct->numseq[i+1]==3&&ct->numseq[j-1]==4)||(ct->numseq[i+1]==4&&ct->numseq[j-1]==3)) {

								rarray=min(rarray,(erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1)));

							}
							else if (i-1>0&&j+1<2*number) {
								if ((ct->numseq[i-1]==3&&ct->numseq[j+1]==4)||(ct->numseq[i-1]==4&&ct->numseq[j+1]==3)) {

									rarray=min(rarray,(erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1)));

								}
							}
						}
					}

					//Perhaps i,j closes an interior or bulge loop, search for the best
					//possibility

					//If this is intermolecular folding, or if running in parallel with openMP, revert to the old O(N^4) algorithm
					//For now, also use O(N^4) internal loops for applying SHAPE to single-stranded regions
#ifndef SMP
					if (ct->intermolecular||(ct->shaped&&(ct->SHAPEintercept_ss!=0||ct->SHAPEslope_ss!=0))||ct->ssoffset) {
#endif
						if (((j-i-1)>=(minloop+3))||(j>(number))) {
							for (di=(j-i-3);di>=1;di--) {
								for (ip=(i+1);ip<=(j-1-di);ip++) {
									jp = di+ip;
									if ((j-i-2-di)>(data->eparam[7])) goto sub1;
									if (abs(ip-i+j-jp)<=(maxinter)) {
										if (ip>(number)) {
											//if (jp<=number) {

											//	v.f(i,j)=min(v.f(i,j),(erg2(i,j,ip,jp,ct,data,fce[i][ip-number],
											//		fce[jp][j])+
											//		v.f(ip-(number),jp)));

											//}
											//else {
											rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip-number),
															fce.f(jp-number,j-number))+
														v.f(ip-(number),jp-(number))));

											if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]) {
												//ip or jp is modified

												rarray=min(rarray,erg2(i,j,ip,jp,ct,data,fce.f(i,ip-number),
															fce.f(jp-number,j-number))+
														v.f(ip-(number)+1,jp-(number)-1)+
														erg1(ip-number,jp-number,ip+1-number,jp-1-number,ct,data));

											}
											//}
										}
										else {
											if (jp<=number) {




												rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),
																fce.f(jp,j))+
															v.f(ip,jp)));

												if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]) {
													//i or j is modified
													rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),
																	fce.f(jp,j))+
																v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));

												}




											}
											else {


												rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp-number,j-number))+
															v.f(ip,jp)));


												if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]) {
													//i or j is modified
													rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp-number,j-number))+
																v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));

												}

											}




										}
									}



								}
							}
						}
#ifndef SMP
					}
					/* Perhaps i,j closes an interior or bulge loop, search for the best possibility
					   fill the interior loops' energy rarray first
					   calculate the small loop (size<=5) first
					   the larger loop is prefilled with curE[dp][i] (after sub2)

					   d= j-i, dp= jp-ip (interior loop)
					   i<ip<number<jp<j or i<ip<jp<j<number
					 */

					//Below is the O(N^3) algorithm
					else {
						if ((d-1)>=(minloop+3)||j>number) 
							for (dp=d-3;dp>=((j>number)?1:minloop+1);dp--) {
								ll=d-dp-2;//ll is the loop length in terms of number of unpaired nucs

								//calculate every ip,jp when ll <=5: this includes loops with special rules:
								//0x1,0x2,0x3,0x4,0x5,1x1,1x2,1x3,1x4,2x2,2x3
								if(ll>=1&&ll<=5)
								{

									for (ip=i+1;ip<=j-1-dp;ip++)
									{
										jp=ip+dp;
										if (inc[ct->numseq[ip]][ct->numseq[jp]])
											if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
											{
												//using jpf and jf and  instead of jp and j when j,jp>number
												jpf=( (jp<=number)?jp:jp-number);
												jf=((j<=number)?j:j-number);
												rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),
																fce.f(jpf,jf))+v.f(ip,jp)));
												//i or j is modified
												if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce.f(ip,jp)&SINGLE)) 
													rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jpf,jf))+
																v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));

											}

									}
								}
								//when size >=6 and <=30;
								//	else if (ll>=6&&ll<=(data->eparam[7]))
								else if (ll>=6&&ll<=maxinter)
								{

									//interior energy prefilled (after sub2:) + exterior engergy (the function for erg2in and erg2ex is found in rna_library_inter.cpp
									//The interior energy, curE, is the lowest free energy for a interior fragment with a pair closing dp nucs
									//the interior fragment and asymetry contriubutions are included
									rarray=min(rarray,(curE[dp][i]+erg2ex(i,j,ll,ct,data)));

									//considering loop 1x(n-1) (bl=1)and 0 x n(bulge)(bl=0) as stacking bonus on 1x(n-1) and bulge is not allowed
									for (bl=0;bl<=1;bl++)
									{
										ip=i+1+bl;
										jp=ip+dp;
										jpf=( (jp<=number)?jp:jp-number);
										jf=((j<=number)?j:j-number);
										if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
											if (inc[ct->numseq[ip]][ct->numseq[jp]])
												if (abs(ip-i+jp-j)<=maxasym) 
												{
													rarray = min(rarray, (erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jpf,jf))+v.f(ip,jp)));
													//i or j is modified
													if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce.f(ip,jp)&SINGLE)) 
														rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jpf,jf))+
																	v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));

												}

										jp=j-1-bl;
										ip=jp-dp;
										jpf=( (jp<=number)?jp:jp-number);
										jf=((j<=number)?j:j-number);
										if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
											if (inc[ct->numseq[ip]][ct->numseq[jp]])
												if (abs(ip-i+jp-j)<=maxasym)
												{
													rarray = min(rarray, (erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jpf,jf))+v.f(ip,jp)));
													if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce.f(ip,jp)&SINGLE)) 
														rarray=min(rarray,(erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jpf,jf))+
																	v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));

												} 
									}

								}
							}	  

					}
#endif


					//Perhaps i,j closes a multibranch or exterior loop, search for the best possibility



sub1:

#ifndef INSTRUMENTED//

#else//If INSTRUMENTED IS DEFINED
					v.f(i,j) = rarray;
					rarray = INFINITE_ENERGY;
#endif//END INSTRUMENTED


					if (((j-i-1)>=(2*minloop+4))||(j>(number))) {
						for (ii=1;ii<=4;ii++) e[ii]=INFINITE_ENERGY;




						//consider the exterior loop closed by i,j
						if (j>number) {
							rarray = min(rarray,w3[i+1] + w5[j-number-1] + penalty(i,j,ct,data));


							if (i!=number) rarray = min(rarray,erg4(i,j,i+1,1,ct,data,lfce[i+1]) + penalty(i,j,ct,data)+w3[i+2] + w5[j-number-1] );
							if (j!=(number+1)) rarray = min(rarray,erg4(i,j,j-1,2,ct,data,lfce[j-1]) +penalty(i,j,ct,data)+ w3[i+1] + w5[j-number-2] );
							if ((i!=number)&&(j!=(number+1))) {
								rarray = min(rarray,data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] 
										+checknp(lfce[i+1],lfce[j-1]) + w3[i+2] + w5[j-number-2]
										+penalty(i,j,ct,data)+ ct->SHAPEss_give_value(i+1)+ ct->SHAPEss_give_value(j-1));

							}


							//consider the coaxial stacking of a helix from i to ip onto helix ip+1 or ip+2 to j:

							//first consider a helix stacking from the 5' sequence fragment:
#ifndef disablecoax
							for (ip=j-number-minloop-1;ip>0;ip--) {
								//first consider flush stacking
								rarray = min(rarray,
										w3[i+1]+w5[ip-1]+penalty(i,j,ct,data)+penalty(j-number-1,ip,ct,data)+
										ergcoaxflushbases(ip,j-number-1,j-number,i,ct,data)+v.f(ip,j-number-1));

								if ((mod[ip]||mod[j-number-1])&&j-number-2>0&&!(fce.f(ip,j-number-1)&SINGLE)) {
									if (inc[ct->numseq[ip+1]][ct->numseq[j-number-2]]) {
										rarray = min(rarray,
												w3[i+1]+w5[ip-1]+penalty(i,j,ct,data)+penalty(j-number-1,ip,ct,data)+
												ergcoaxflushbases(ip,j-number-1,j-number,i,ct,data)+v.f(ip+1,j-number-2)+
												erg1(ip,j-number-1,ip+1,j-number-2,ct,data));
									}

								}


								if (j-number-2>0) {
									//now consider an intervening nuc
									if(i<number) {
										rarray = min(rarray,
												w3[i+2]+w5[ip-1]+penalty(i,j,ct,data)+penalty(ip,j-number-2,ct,data)+
												ergcoaxinterbases2(ip,j-number-2,j-number,i,ct,data)+v.f(ip,j-number-2)+checknp(lfce[i+1],lfce[j-number-1])
											    );


										if ((mod[ip]||mod[j-number-2])&&inc[ct->numseq[ip+1]][ct->numseq[j-number-3]]!=0&&!(fce.f(ip,j-number-2)&SINGLE)) {
											rarray = min(rarray,
													w3[i+2]+w5[ip-1]+penalty(i,j,ct,data)+penalty(ip,j-number-2,ct,data)+
													ergcoaxinterbases2(ip,j-number-2,j-number,i,ct,data)+v.f(ip+1,j-number-3)+
													erg1(ip,j-number-2,ip+1,j-number-3,ct,data)+checknp(lfce[i+1],lfce[j-number-1])
												    );


										}
									}


									//consider the other possibility for an intervening nuc
									rarray = min(rarray,
											w3[i+1]+w5[ip-1]+penalty(i,j,ct,data)+penalty(ip+1,j-number-2,ct,data)+
											ergcoaxinterbases1(ip+1,j-number-2,j-number,i,ct,data)+v.f(ip+1,j-number-2)+checknp(lfce[ip],lfce[j-number-1])
										    );


									if ((mod[ip+1]||mod[j-number-2])&&inc[ct->numseq[ip+2]][ct->numseq[j-number-3]]&&!(fce.f(ip+1,j-number-2)&SINGLE)) {
										rarray = min(rarray,
												w3[i+1]+w5[ip-1]+penalty(i,j,ct,data)+penalty(ip+1,j-number-2,ct,data)+
												ergcoaxinterbases1(ip+1,j-number-2,j-number,i,ct,data)+v.f(ip+2,j-number-3)
												+erg1(ip+1,j-number-2,ip+2,j-number-3,ct,data)+checknp(lfce[ip],lfce[j-number-1])
											    );
									}


								}


							}

							//now consider a helix stacking from the 3' sequence fragment:
							for (ip=i+minloop+1;ip<=number;ip++) {
								//first consider flush stacking
								rarray = min(rarray,
										w3[ip+1]+w5[j-number-1]+penalty(i,j,ct,data)+penalty(ip,i+1,ct,data)+
										ergcoaxflushbases(j-number,i,i+1,ip,ct,data)+v.f(i+1,ip));


								if (mod[i+1]||mod[ip]&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&!(fce.f(i+1,ip)&SINGLE)) {

									rarray = min(rarray,
											w3[ip+1]+w5[j-number-1]+penalty(i,j,ct,data)+penalty(ip,i+1,ct,data)+
											ergcoaxflushbases(j-number,i,i+1,ip,ct,data)+v.f(i+2,ip-1)
											+erg1(i+1,ip,i+2,ip-1,ct,data));

								}

								//now consider an intervening nuc
								if (j-number>1) {
									rarray = min(rarray,
											w3[ip+1]+w5[j-number-2]+penalty(i,j,ct,data)+penalty(ip,i+2,ct,data)+
											ergcoaxinterbases1(j-number,i,i+2,ip,ct,data)+v.f(i+2,ip)+checknp(lfce[i+1],lfce[j-number-1])
										    );


									if ((mod[i+2]||mod[ip])&&inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&!(fce.f(i+2,ip)&SINGLE)) {

										rarray = min(rarray,
												w3[ip+1]+w5[j-number-2]+penalty(i,j,ct,data)+penalty(ip,i+2,ct,data)+
												ergcoaxinterbases1(j-number,i,i+2,ip,ct,data)+v.f(i+3,ip-1)
												+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[j-number-1])
											    );

									}
								}

								//consider the other possibility for an intervening nuc
								rarray = min(rarray,
										w3[ip+1]+w5[j-number-1]+penalty(i,j,ct,data)+penalty(ip-1,i+2,ct,data)+
										ergcoaxinterbases2(j-number,i,i+2,ip-1,ct,data)+v.f(i+2,ip-1)+checknp(lfce[i+1],lfce[ip])
									    );

								if ((mod[i+2]||mod[ip-1])&&inc[ct->numseq[i+3]][ct->numseq[ip-2]]&&!(fce.f(i+2,ip-1)&SINGLE)) {

									rarray = min(rarray,
											w3[ip+1]+w5[j-number-1]+penalty(i,j,ct,data)+penalty(ip-1,i+2,ct,data)+
											ergcoaxinterbases2(j-number,i,i+2,ip-1,ct,data)+v.f(i+3,ip-2)
											+erg1(i+2,ip-1,i+3,ip-2,ct,data)+checknp(lfce[i+1],lfce[ip])
										    );
								}
							}
#endif //ifndef disablecoax
						}

#ifndef INSTRUMENTED

#else//IF INSTRUMENTED IS DEFINED
						v.f(i,j) = min(rarray,v.f(i,j));//QUESTION: DO WE NEED THIS? WAS IN algorithm.napss.cpp
						vext.f(i,j) = rarray;

						rarray = INFINITE_ENERGY;//QUESTION: DO WE NEED THIS? WAS IN algorithm.napss.cpp
#endif//END INSTRUMENTED

						//consider the multiloop closed by i,j
						if ((j-i)>(2*minloop+4)&&i!=number) {
							//no dangling ends on i-j pair:
							if (j-1!=number) {
								rarray = min(rarray,wmb.f(i+1,j-1)+data->eparam[5]+data->eparam[10]
										+ penalty(i,j,ct,data));


								//i+1 dangles on i-j pair:

								if (i+1!=number) rarray = min(rarray,erg4(i,j,i+1,1,ct,data,lfce[i+1]) + penalty(i,j,ct,data) +
										wmb.f(i+2,j-1) + data->eparam[5] + data->eparam[6] + data->eparam[10]);
							}
							if (j-2!=number) {
								//j-1 dangles
								if (j!=(number+1)) rarray = min(rarray,erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
										wmb.f(i+1,j-2) + data->eparam[5] + data->eparam[6] + data->eparam[10]);
								//both i+1 and j-1 dangle
								if ((i+1!=number)&&(j!=(number+1))) {
									rarray = min(rarray,data->tstkm[ct->numseq[i]][ct->numseq[j]]
											[ct->numseq[i+1]][ct->numseq[j-1]] +
											checknp(lfce[i+1],lfce[j-1])+
											wmb.f(i+2,j-2) + data->eparam[5] + 2*data->eparam[6] + data->eparam[10]
											+penalty(i,j,ct,data)+ ct->SHAPEss_give_value(i+1)+ ct->SHAPEss_give_value(j-1));
								}
							}

							//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
#ifndef disablecoax
							for (ip=i+1;(ip<j);ip++) {
								//first consider flush stacking

								//conditions guarantee that the coaxial stacking isn't considering an exterior loop 
								//if ((i!=number)/*&&(i+1!=number)*//*&&((j>number)||(ip!=number)&&(ip+1!=number))&&(j-1!=number)*/) 
								if (i!=number&&ip!=number&&j-1!=number) {
									
									if (inc[ct->numseq[i+1]][ct->numseq[ip]]) {
										//only proceed if i+1 can pair with ip

										rarray = min(rarray,penalty(i,j,ct,data)+v.f(i+1,ip)+
												penalty(i+1,ip,ct,data)+data->eparam[5]
												+2*data->eparam[10]+w.f(ip+1,j-1)+ergcoaxflushbases(j,i,i+1,ip,ct,data));


										if((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&!(fce.f(i+1,ip)&SINGLE)) {

											rarray = min(rarray,penalty(i,j,ct,data)+v.f(i+2,ip-1)+
													penalty(i+1,ip,ct,data)+data->eparam[5]
													+2*data->eparam[10]+w.f(ip+1,j-1)+ergcoaxflushbases(j,i,i+1,ip,ct,data)
													+erg1(i+1,ip,i+2,ip-1,ct,data));
										}
									}
									//if ((ip<j-1)&&(i+2!=number)) 
									if (inc[ct->numseq[i+2]][ct->numseq[ip]]&&ip+2<j-1&&i+1!=number&&ip+1!=number) {
										//now consider an intervening nuc
										if ((ip+2<j-1)/*&&(j>number||ip+2!=number)*/)
											rarray = min(rarray,penalty(i,j,ct,data)+v.f(i+2,ip)+
													penalty(i+2,ip,ct,data)+data->eparam[5]
													+2*data->eparam[6]+2*data->eparam[10]+w.f(ip+2,j-1)
													+ergcoaxinterbases2(j,i,i+2,ip,ct,data)+checknp(lfce[i+1],lfce[ip+1])
												    );

										if((mod[i+2]||mod[ip])&&inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&!(fce.f(i+2,ip)&SINGLE)) {

											rarray = min(rarray,penalty(i,j,ct,data)+v.f(i+3,ip-1)+
													penalty(i+2,ip,ct,data)+data->eparam[5]
													+2*data->eparam[6]+2*data->eparam[10]+w.f(ip+2,j-1)
													+ergcoaxinterbases2(j,i,i+2,ip,ct,data)
													+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[ip+1])
												    );

										}

										if (ip+1<j-2&&j-2!=number)
											rarray = min(rarray,penalty(i,j,ct,data)+v.f(i+2,ip)+
													penalty(i+2,ip,ct,data)+data->eparam[5]
													+2*data->eparam[6]+2*data->eparam[10]+w.f(ip+1,j-2)
													+ergcoaxinterbases1(j,i,i+2,ip,ct,data)+checknp(lfce[i+1],lfce[j-1])
												    );

										if((mod[i+2]||mod[ip])&&inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&!(fce.f(i+2,ip)&SINGLE)) {

											rarray = min(rarray,penalty(i,j,ct,data)+v.f(i+3,ip-1)+
													penalty(i+2,ip,ct,data)+data->eparam[5]
													+2*data->eparam[6]+2*data->eparam[10]+w.f(ip+1,j-2)
													+ergcoaxinterbases1(j,i,i+2,ip,ct,data)
													+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[j-1])
												    );

										}
									}//end (ip+2<j-1&&i+1!=number&&ip+1!=number)
								}
							}

							//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
							for (ip=j-1;ip>i;ip--) {

								//conditions guarantee that the coaxial stacking isn't considering an exterior loop
								//if ((i!=number)&&(i+1!=number)&&((j>number)||(ip!=number)&&(ip-1!=number))&&(j-1!=number)) {}
								if (j-1!=number&&ip-1!=number&&i!=number) {
									//first consider flush stacking

									if (inc[ct->numseq[j-1]][ct->numseq[ip]]) {
										//only proceed if j-1 can pair to ip

										rarray = min(rarray,penalty(i,j,ct,data)+v.f(ip,j-1)+
												penalty(j-1,ip,ct,data)+data->eparam[5]
												+2*data->eparam[10]+w.f(i+1,ip-1)+ergcoaxflushbases(ip,j-1,j,i,ct,data));


										if((mod[ip]||mod[j-1])&&inc[ct->numseq[ip+1]][ct->numseq[j-2]]&&!(fce.f(ip,j-1)&SINGLE)) {
											rarray = min(rarray,penalty(i,j,ct,data)+v.f(ip+1,j-2)+
													penalty(j-1,ip,ct,data)+data->eparam[5]
													+2*data->eparam[10]+w.f(i+1,ip-1)+ergcoaxflushbases(ip,j-1,j,i,ct,data)
													+erg1(ip,j-1,ip+1,j-2,ct,data));

										}
									}

									if  (inc[ct->numseq[j-2]][ct->numseq[ip]]&&j-2!=number) {
										//now consider an intervening nuc
										//if ((ip>i+1)&&(j>number||ip-2!=number))
										if (ip-2>i+1&&ip-2!=number) {
											rarray = min(rarray,penalty(i,j,ct,data)+v.f(ip,j-2)+
													penalty(j-2,ip,ct,data)+data->eparam[5]
													+2*data->eparam[6]+2*data->eparam[10]+w.f(i+1,ip-2)
													+ergcoaxinterbases1(ip,j-2,j,i,ct,data)+checknp(lfce[j-1],lfce[ip-1])
												    );



											if((mod[ip]||mod[j-2])&&inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&!(fce.f(ip,j-2)&SINGLE)) {
												rarray = min(rarray,penalty(i,j,ct,data)+v.f(ip+1,j-3)+
														penalty(j-2,ip,ct,data)+data->eparam[5]
														+2*data->eparam[6]+2*data->eparam[10]+w.f(i+1,ip-2)
														+ergcoaxinterbases1(ip,j-2,j,i,ct,data)
														+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[ip-1])
													    );

											}
										}

										if ((ip-1>i+2)&&i+1!=number) {
											rarray = min(rarray,penalty(i,j,ct,data)+v.f(ip,j-2)+
													penalty(j-2,ip,ct,data)+data->eparam[5]
													+2*data->eparam[6]+2*data->eparam[10]+w.f(i+2,ip-1)
													+ergcoaxinterbases2(ip,j-2,j,i,ct,data)+checknp(lfce[j-1],lfce[i+1])
												    );

											if((mod[ip]||mod[j-2])&&inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&!(fce.f(ip,j-2)&SINGLE)) {
												rarray = min(rarray,penalty(i,j,ct,data)+v.f(ip+1,j-3)+
														penalty(j-2,ip,ct,data)+data->eparam[5]
														+2*data->eparam[6]+2*data->eparam[10]+w.f(i+2,ip-1)
														+ergcoaxinterbases2(ip,j-2,j,i,ct,data)
														+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[i+1])
													    );

											}
										}
									}
								}
							}
#endif //ifndef disablecoax
							if (ct->intermolecular) {

								//intermolecular, so consider wmb2,
								//don't add the multiloop penalties because this is a exterior loop

								rarray = min(rarray,wmb2->f(i+1,j-1) + penalty(i,j,ct,data)+INFINITE_ENERGY);


								//i+1 dangles on i-j pair:
								if (i!=number) rarray = min(rarray,erg4(i,j,i+1,1,ct,data,lfce[i+1]) + penalty(i,j,ct,data) +
										wmb2->f(i+2,j-1)+INFINITE_ENERGY);
								//j-1 dangles
								if (j!=(number+1)) rarray = min(rarray,erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
										wmb2->f(i+1,j-2)+INFINITE_ENERGY);
								//both i+1 and j-1 dangle
								if ((i!=number)&&(j!=(number+1))) {
									rarray = min(rarray,
											data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
											checknp(lfce[i+1],lfce[j-1]) +
											wmb2->f(i+2,j-2) + penalty(i,j,ct,data)+INFINITE_ENERGY);

								}
							}
						}
					}

#ifndef INSTRUMENTED//IF INSTRUMENTED IS NOT DEFINED
					v.f(i,j) = rarray;
#else//IF INSTRUMENTED IS DEFINED
					vmb.f(i,j) = rarray;
					v.f(i,j) = min(rarray,v.f(i,j));//QUESTION: THIS IS DIFFERENT FROM THE MAIN algorithm.cpp (see #ifndef above)
#endif//END INSTRUMENTED
sub2:

					/*prefill curE[i] and prev[i] for the first two diagonals as ll =4 and 5	
					  As d =10, only fill curE (ll=4, d=10)
					  As d =11, only fill prevE (ll=4||5, d=11)
					  As d>11, fill curE(ll=4||5, d>11)
					  exchange curE and prevE after d >11  as curE[h][i][dp]=curE=curE[h-2][i+1][dp]
					  (curE[i][j][ll]=curE[i+1][j-1][ll-2])

note: dp = jp-ip, where ip and jp are the interior pairs in interior loop
d = j - i;
					 */


					if (((d-1)>=(minloop+3)||j>number)&&!ct->intermolecular) 
						for (dp=d-3;dp>=((j>number)?1:minloop+1);dp--)
						{
							ll=d-dp-2;
							//calculate every ip>ip+1,jp<jp-1 when ll ==5 ||4
							if (ll==4||ll==5) 
							{
								for (ip=i+2;ip<=j-2-dp;ip++)
								{
									jp=ip+dp;
									if (inc[ct->numseq[ip]][ct->numseq[jp]])
										if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
										{
											if(d==( (j>number)?7:10 )||d>( (j>number)?8:11 ))  
												//fill the first diagonal of d and first two of ll for every larger d 
											{
												curE[dp][i]=min(curE[dp][i],(erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+v.f(ip,jp)));

												//i or j is modified
												if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce.f(ip,jp)&SINGLE)) 
													curE[dp][i]=min(curE[dp][i],(erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));

											}


											else if ( d==((j>number)?8:11) )  //fill the second diagonal of d
											{
												prevE[dp][i]=min(prevE[dp][i],(erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+v.f(ip,jp)));

												if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce.f(ip,jp)&SINGLE)) 

													prevE[dp][i]=min(prevE[dp][i],(erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+ v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));
											}
										}
								}
							}

							//when size >=6 and <=30;
							//   else if (ll>=6&&ll<=(data->eparam[7]))
							else if (ll>=6&&ll<=maxinter){
								//calculate minimum curE[dp][i] of 1 x (n-1) for next step's 2 x (n-2)
								ip=i+2;
								jp=ip+dp;
								if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
									if (abs(ip-i+jp-j)<=maxasym) 
										if (inc[ct->numseq[ip]][ct->numseq[jp]]){
											curE[dp][i] = min(curE[dp][i], (erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+v.f(ip,jp))) ;	

											if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce.f(ip,jp)&SINGLE)) 
												curE[dp][i]=min(curE[dp][i],(erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+
															v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));


										}

								jp=j-2;
								ip=jp-dp;
								if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
									if (abs(ip-i+jp-j)<=maxasym) 
										if (inc[ct->numseq[ip]][ct->numseq[jp]]){
											curE[dp][i] = min(curE[dp][i], (erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+v.f(ip,jp))) ;	

											if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce.f(ip,jp)&SINGLE)) 
												curE[dp][i]=min(curE[dp][i],(erg2in(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+
															v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)));


										}
							}
						}	  

					//also block propagation of interior loops that contain nucleotides that need to be double-stranded:
					if ((lfce[i]||lfce[j])&&!ct->intermolecular) for (dp=1;dp<=d;dp++) curE[dp][i] = INFINITE_ENERGY;//QUESTION: THIS WASN'T IN THE algirithm.napss.cpp


					//Compute w[i][j]: best energy between i and j where i,j does not have
					//	to be a base pair
					//(an exterior loop when it contains n and 1 (ie:n+1)   )


					w.f(i,j)=INFINITE_ENERGY;

					if (fce.f(i,j)&PAIR)  {//force a pair between i and j
						w.f(i,j) = v.f(i,j)+data->eparam[10]+penalty(i,j,ct,data);

						goto sub3;
					}

					for (ii=1;ii<=5;ii++) e[ii] = INFINITE_ENERGY;

					if (i!=number) {
						//calculate the energy of i stacked onto the pair of i+1,j

						e[1] = v.f(i+1,j) + data->eparam[10] + data->eparam[6] +
							erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data);

						if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[j-1]]&&!(fce.f(i+1,j)&SINGLE)) {

							e[1] = min(e[1],v.f(i+2,j-1) + data->eparam[10] + data->eparam[6] +
									erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
									+erg1(i+1,j,i+2,j-1,ct,data));


						}
						if (!lfce[i]) {
							if (!(fce.f(i,i)&INTER))
								//add a nuc to an existing loop:
								e[4] = w.f(i+1,j) + data->eparam[6]+ ct->SHAPEss_give_value(i);
							//this is for when i represents the center of an intermolecular linker:
							else e[4] = w.f(i+1,j) + data->eparam[6] + INFINITE_ENERGY;
						}
					}
					if (j!=((number)+1)) {
						//calculate the energy of j stacked onto the pair of i,j-1
						if (j!=1) {
							e[2] = v.f(i,j-1) + data->eparam[10] + data->eparam[6] +
								erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data);

							if ((mod[i]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-2]]&&!(fce.f(i,j-1)&SINGLE)) {

								e[2] = min(e[2],v.f(i+1,j-2) + data->eparam[10] + data->eparam[6] +
										erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
										+erg1(i,j-1,i+1,j-2,ct,data));
							}

							if (!lfce[j]) {
								if (!(fce.f(j,j)&INTER)) {
									//add a nuc to an existing loop:
									e[5] = w.f(i,j-1) + data->eparam[6]+ ct->SHAPEss_give_value(j);
								}
								else e[5] = w.f(i,j-1) + data->eparam[6] + INFINITE_ENERGY;

							}
						}
					}
					if ((i!=(number))&&(j!=((number)+1))) {
						//calculate i and j stacked onto the pair of i+1,j-1
						if (j!=1&&!lfce[i]&&!lfce[j]) {
							e[3] = v.f(i+1,j-1) + data->eparam[10] + 2*(data->eparam[6]) +
								data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
								+penalty(j-1,i+1,ct,data)+ ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(j);



							if ((mod[i+1]||mod[j-1])&&(j-2>0)&&!(fce.f(i+1,j-1)&SINGLE)) {
								if(inc[ct->numseq[i+2]][ct->numseq[j-2]]) {

									e[3] = min(e[3],v.f(i+2,j-2) + data->eparam[10] + 2*(data->eparam[6]) +
											data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
											+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data)+ ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(j));

								}
							}
						}
					}

					//fragment with i paired to j
					e[1] = min(((data->eparam[10])+v.f(i,j)+penalty(j,i,ct,data)),e[1]);

					if ((mod[i]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[j-1]]&&!(fce.f(i,j)&SINGLE)) {

						e[1] = min(((data->eparam[10])+v.f(i+1,j-1)+penalty(j,i,ct,data)+erg1(i,j,i+1,j-1,ct,data)),e[1]);			

					}
					w.f(i,j) = min(e[1],e[2]);
					w.f(i,j) = min(w.f(i,j),e[3]);
					w.f(i,j) = min(w.f(i,j),e[4]);
					w.f(i,j) = min(w.f(i,j),e[5]);

					if (ct->intermolecular) {

						//wmb2[i][j%3] = INFINITE_ENERGY;
						//keep track of w2:
						for (ii=1;ii<=5;ii++) e[ii] = 2*INFINITE_ENERGY;

						if (i!=number) {
							//calculate the energy of i stacked onto the pair of i+1,j

							e[1] = v.f(i+1,j) +
								erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data);

							if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[j-1]]&&!(fce.f(i+1,j)&SINGLE)) {

								e[1] = min(e[1],v.f(i+2,j-1) +
										erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
										+erg1(i+1,j,i+2,j-1,ct,data));
							}

							//if (!lfce[i]) {
							//if (!(fce.f(i,i)&DUBLE))
							e[4] = w2->f(i+1,j);
							//this is for when i represents the center of an intermolecular linker:
							//else e[4] = w2->f(i+1,j) - INFINITE_ENERGY + data->init;
							//}
						}
						if (j!=((number)+1)) {
							//calculate the energy of j stacked onto the pair of i,j-1
							if (j!=1) {
								e[2] = v.f(i,j-1)   +
									erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data);

								if ((mod[i]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-2]]&&!(fce.f(i,j-1)&SINGLE)) {

									e[2] = min(e[2],v.f(i+1,j-2) +
											erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
											+erg1(i,j-1,i+1,j-2,ct,data));

								}

								//if (!lfce[j]) {
								//	if (!(fce[j][j]&DUBLE)) {
								e[5] = w2->f(i,j-1);
								//  }
								//else e[5] = w2->f(i,j-1) - INFINITE_ENERGY + data->init;

								//}
							}
						}
						if ((i!=(number))&&(j!=((number)+1))) {
							//calculate i and j stacked onto the pair of i+1,j-1
							if (j!=1) {
								e[3] = v.f(i+1,j-1)   +
									data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
									[ct->numseq[j]][ct->numseq[i]]
									+checknp(lfce[i],lfce[j])
									+penalty(j-1,i+1,ct,data);



								if ((mod[i+1]||mod[j-1])&&inc[ct->numseq[i+2]][ct->numseq[j-2]]&&!(fce.f(i+1,j-1)&SINGLE)) {

									e[3] = min(e[3],v.f(i+2,j-2) +
											data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
											+checknp(lfce[i],lfce[j])
											+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data));

								}
							}
						}

						e[1] = min(e[1],(v.f(i,j)+penalty(j,i,ct,data)));

						if (mod[i]||mod[j]&&inc[ct->numseq[i+1]][ct->numseq[j-1]]&&!(fce.f(i,j)&SINGLE)) {

							e[1] = min((v.f(i+1,j-1)+penalty(j,i,ct,data)+erg1(i,j,i+1,j-1,ct,data)),e[1]);			

						}

						w2->f(i,j) = min(e[1],e[2]);
						w2->f(i,j) = min(w2->f(i,j),e[3]);
						w2->f(i,j) = min(w2->f(i,j),e[4]);
						w2->f(i,j) = min(w2->f(i,j),e[5]);

					}

					/*if (((j-i-1)>(2*minloop+2))||(j>(number))) {
					// search for an open bifuraction:
					for (k=i;k<=j-1;k++) {
					if (k==(number)) w.f(i,j)=min(w.f(i,j),
					w3[i]+w5[j-(number)]);
					else w.f(i,j) = min(w.f(i,j),
					w.f(i,k)+work[k+1][j%3]);
					}
					}  */
					////fill wmb:

					if (((j-i-1)>(2*minloop+2))||j>number) {

						//search for an open bifurcation:
						
						//in this code, the w arrays are accessed directly, and not through the f function
						int end=min(number,j);
						for (k=i;k<end;++k) {
							wmb.f(i,j) = min(wmb.f(i,j),w.dg[i][k]+w.dg[k+1][j]);
						}		
						for (k=number;k<j;++k) {
							wmb.f(i,j) = min(wmb.f(i,j),w.dg[i][k]+w.dg[k+1-number][j-number]);
						}


						

						//for the sake of coaxial stacking, also consider the addition of nucs
						//to a previously calculated wmb
						if (i!=number)
							if (!lfce[i]) wmb.f(i,j) = min(wmb.f(i,j) ,wmb.f(i+1,j) +data->eparam[6]+ ct->SHAPEss_give_value(i));
						if (j!=number+1)
							if (!lfce[j]) wmb.f(i,j)  = min(wmb.f(i,j) ,wmb.f(i,j-1) +data->eparam[6]+ ct->SHAPEss_give_value(j));

						e[1]=INFINITE_ENERGY;
						e[2]=INFINITE_ENERGY;
						//also consider the coaxial stacking of two helixes
#ifndef disablecoax
						for (ip=i+minloop+1;ip<j-minloop-1;ip++) {
							//first consider flush stacking


							if (ip!=number) {
								if (inc[ct->numseq[i]][ct->numseq[ip]]&&inc[ct->numseq[j]][ct->numseq[ip+1]]) {
									//ony proceed if i and ip and ip+1 and j can pair
									e[1]=min(e[1],v.f(i,ip)+v.f(ip+1,j)+penalty(i,ip,ct,data)
											+penalty(ip+1,j,ct,data)+ergcoaxflushbases(i,ip,ip+1,j,ct,data));

									if (mod[i]||mod[ip]||mod[ip+1]||mod[j]) {

										if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
												&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&!(fce.f(i,ip)&SINGLE)&&!(fce.f(ip+1,j)&SINGLE)) {

											e[1]=min(e[1],v.f(i+1,ip-1)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
													+penalty(ip+1,j,ct,data)+ergcoaxflushbases(i,ip,ip+1,j,ct,data)
													+erg1(i,ip,i+1,ip-1,ct,data)+erg1(ip+1,j,ip+2,j-1,ct,data));
										}

										if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&!(fce.f(i,ip)&SINGLE)) {

											e[1]=min(e[1],v.f(i+1,ip-1)+v.f(ip+1,j)+penalty(i,ip,ct,data)
													+penalty(ip+1,j,ct,data)+ergcoaxflushbases(i,ip,ip+1,j,ct,data)
													+erg1(i,ip,i+1,ip-1,ct,data));
										}

										if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&!(fce.f(ip+1,j)&SINGLE)) {


											e[1]=min(e[1],v.f(i,ip)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
													+penalty(ip+1,j,ct,data)+ergcoaxflushbases(i,ip,ip+1,j,ct,data)
													+erg1(ip+1,j,ip+2,j-1,ct,data));

										}
									}
								}//end (inc[ct->numseq[i]][ct->numseq[ip]]&&inc[ct->numseq[j]][ct->numseq[ip+1]]) 

								if (ip+1!=number) {
									if (inc[ct->numseq[i]][ct->numseq[ip]]&&inc[ct->numseq[j-1]][ct->numseq[ip+2]]) {
										//only proceed if i and ip and j-1 and ip+2 can pair
										if (!lfce[ip+1]&&!lfce[j]) {
											//now consider an intervening mismatch
										
											e[2]=min(e[2],v.f(i,ip)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
													+penalty(ip+2,j-1,ct,data)+ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
													);

											if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
												if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
														&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&!(fce.f(i,ip)&SINGLE)&&!(fce.f(ip+2,j-1)&SINGLE)) {

													if(!lfce[ip+1]&&!lfce[j]) e[2]=min(e[2],v.f(i+1,ip-1)+v.f(ip+3,j-2)+penalty(i,ip,ct,data)
															+penalty(ip+2,j-1,ct,data)+ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
															+erg1(i,ip,i+1,ip-1,ct,data)+erg1(ip+2,j-1,ip+3,j-2,ct,data)
															);


												}

												if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&!(fce.f(i,ip)&SINGLE)) {

													if(!lfce[ip+1]&&!lfce[j]) e[2]=min(e[2],v.f(i+1,ip-1)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
															+penalty(ip+2,j-1,ct,data)+ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
															+erg1(i,ip,i+1,ip-1,ct,data));


												}

												if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&!(fce.f(ip+2,j-1)&SINGLE)) {


													if(!lfce[ip+1]&&!lfce[j]) e[2]=min(e[2],v.f(i,ip)+v.f(ip+3,j-2)+penalty(i,ip,ct,data)
															+penalty(ip+2,j-1,ct,data)+ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
															+erg1(ip+2,j-1,ip+3,j-2,ct,data));

												}
											}
										}//end (!lfce[ip+1]&&!lfce[j])
									}//end (inc[ct->numseq[i]][ct->numseq[ip]]&&inc[ct->numseq[j-1]][ct->numseq[ip+2]])


									if (inc[ct->numseq[i+1]][ct->numseq[ip]]&&inc[ct->numseq[j]][ct->numseq[ip+2]]) {
										//only proceed if i and ip j and ip+2 can pair
										if(!lfce[i]&&!lfce[ip+1]&&i!=number) {
											e[2]=min(e[2],v.f(i+1,ip)+v.f(ip+2,j)+penalty(i+1,ip,ct,data)
													+penalty(ip+2,j,ct,data)+ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
												);

											if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
												if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
														&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&!(fce.f(i+1,ip)&SINGLE)&&!(fce.f(ip+2,j)&SINGLE)) {

													e[2]=min(e[2],v.f(i+2,ip-1)+v.f(ip+3,j-1)+penalty(i+1,ip,ct,data)
															+penalty(ip+2,j,ct,data)+ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
															+erg1(i+1,ip,i+2,ip-1,ct,data)+erg1(ip+2,j,ip+3,j-1,ct,data)
														);	



												}
												if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&!(fce.f(i+1,ip)&SINGLE)) {

													e[2]=min(e[2],v.f(i+2,ip-1)+v.f(ip+2,j)+penalty(i+1,ip,ct,data)
															+penalty(ip+2,j,ct,data)+ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
															+erg1(i+1,ip,i+2,ip-1,ct,data)
														);


												}

												if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&!(fce.f(ip+2,j)&SINGLE)) {


													e[2]=min(e[2],v.f(i+1,ip)+v.f(ip+3,j-1)+penalty(i+1,ip,ct,data)
															+penalty(ip+2,j,ct,data)+ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
															+erg1(ip+2,j,ip+3,j-1,ct,data)
														);

												}
											}
										}
									}//end (inc[ct->numseq[i+1]][ct->numseq[ip]]&&inc[ct->numseq[j]][ct->numseq[ip+2]]) 
								}
							}



						}
#endif //ifndef disablecoax

						wmb.f(i,j) =min(wmb.f(i,j) ,e[1]+2*data->eparam[10]);
						wmb.f(i,j) =min(wmb.f(i,j) ,e[2]+2*data->eparam[10]+2*data->eparam[6]);

						if (j<=number)
							wca[i][j] = min(e[1],e[2]);


						w.f(i,j) = min(w.f(i,j),wmb.f(i,j) );

						if (ct->intermolecular) {
							//intermolecular folding:

							//search for an open bifurcation:
							for (k=i;k<=j;k++) {

								if (k!=number) wmb2->f(i,j) = min(wmb2->f(i,j),w2->f(i,k)+w2->f(k+1,j));
							}

							if (i!=number) {
								if (!(fce.f(i,i)&INTER)) wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) );
								else  wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) + data->init - INFINITE_ENERGY);
							}

							if (j!=number+1) {
								if (!(fce.f(j,j)&INTER)) wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1));
								else wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1) +data->init-INFINITE_ENERGY);
							}


							w2->f(i,j) = min(w2->f(i,j),wmb2->f(i,j) );

						}


					}






sub3:

					//Calculate vmin, the best energy for the entire sequence
					if (j>(number)) {
						vmin = min(vmin,v.f(i,j)+v.f(j-(number),i));
					}

					//Compute w5[i], the energy of the best folding from 1->i, and
					//w3[i], the energy of the best folding from i-->GetSequenceLength()
					if (i==1&&j<=number)
					{  
						if (j<=minloop+1) {
							if (lfce[j]) w5[j]= INFINITE_ENERGY;
							else  w5[j] = w5[j-1]+ ct->SHAPEss_give_value(j);
						}

						else {
							if (lfce[j]) w5[j] = INFINITE_ENERGY;

							else w5[j] = w5[j-1] + ct->SHAPEss_give_value(j);

							for (k=1;k<=5;k++) e[k] = INFINITE_ENERGY;//e[k]=0;
#ifndef disablecoax
							castack = INFINITE_ENERGY;
#endif //ifndef disablecoax
							for (k=0;k<=(j-4);k++) {

								e[1] = min(e[1],(w5[k]+v.f(k+1,j)+penalty(j,k+1,ct,data)));

								if (mod[k+1]||mod[j]&&inc[ct->numseq[k+2]][ct->numseq[j-1]]&&!(fce.f(k+1,j)&SINGLE)) {

									e[1] = min(e[1],(w5[k]+v.f(k+2,j-1)+penalty(j,k+1,ct,data)
												+erg1(k+1,j,k+2,j-1,ct,data)));
								}



								e[2] = min(e[2],(w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v.f(k+2,j)+penalty(j,k+2,ct,data)));

								if(mod[k+2]||mod[j]&&inc[ct->numseq[k+3]][ct->numseq[j-1]]&&!(fce.f(k+2,j)&SINGLE)) {
									e[2] = min(e[2],(w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v.f(k+3,j-1)
												+penalty(j,k+2,ct,data)+erg1(k+2,j,k+3,j-1,ct,data)));

								}


								e[3] = min(e[3],(w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v.f(k+1,j-1)+penalty(j-1,k+1,ct,data)));

								if (mod[k+1]||mod[j-1]&&inc[ct->numseq[k+2]][ct->numseq[j-2]]&&!(fce.f(k+1,j-1)&SINGLE)) {

									e[3] = min(e[3],(w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v.f(k+2,j-2)
												+penalty(j-1,k+1,ct,data)+erg1(k+1,j-1,k+2,j-2,ct,data)));
								}



								e[4] = min(e[4],(w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]]
											[ct->numseq[j]][ct->numseq[k+1]] 
											+checknp(lfce[j],lfce[k+1]) + v.f(k+2,j-1)+
											penalty(j-1,k+2,ct,data))+ ct->SHAPEss_give_value(j)+ ct->SHAPEss_give_value(k+1));

								if (mod[k+2]||mod[j-1]&&inc[ct->numseq[k+3]][ct->numseq[j-2]]&&!(fce.f(k+2,j-1)&SINGLE)) {

									e[4] = min(e[4],(w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]]
												[ct->numseq[j]][ct->numseq[k+1]] 
												+checknp(lfce[j],lfce[k+1]) + v.f(k+3,j-2)+
												penalty(j-1,k+2,ct,data)+erg1(k+2,j-1,k+3,j-2,ct,data))+ ct->SHAPEss_give_value(j)+ ct->SHAPEss_give_value(k+1));

								}







#ifndef disablecoax
								castack = min(castack,w5[k]+wca[k+1][j]);
#endif //ifndef disablecoax


							}
							w5[j] = min(w5[j],e[1]);
							w5[j] = min(w5[j],e[2]);
							w5[j] = min(w5[j],e[3]);
							w5[j] = min(w5[j],e[4]);
#ifndef disablecoax
							w5[j] = min(w5[j],castack);
#endif //ifndef disablecoax

						}		
					}	

					if (j==number)
					{	
						w3[0] = 0;
						w3[number+1] = 0;
						if (i>=(number-minloop)) {    //number+1 ... number-minloop
							if (lfce[i]) w3[i] = INFINITE_ENERGY;
							else w3[i]=w3[i+1]+ ct->SHAPEss_give_value(i);
						}
						//w3[i]=0;
						if (i>=1&&i<=((number)-minloop-1)) {

							if (lfce[i]) w3[i] = INFINITE_ENERGY;

							else w3[i] = w3[i+1]+ ct->SHAPEss_give_value(i);

							for (k=1;k<=5;k++) e[k] = INFINITE_ENERGY;
#ifndef disablecoax
							castack = INFINITE_ENERGY;
#endif //ifndef disablecoax
							for (k=((number)+1);k>=(i+4);k--) {
								e[1] = min(e[1],(v.f(i,k-1)+w3[k]+penalty(k-1,i,ct,data)));

								if((mod[i]||mod[k-1])&&inc[ct->numseq[i+1]][ct->numseq[k-2]]&&!(fce.f(i,k-1)&SINGLE)) {

									e[1] = min(e[1],(v.f(i+1,k-2)+w3[k]+penalty(k-1,i,ct,data)+erg1(i,k-1,i+1,k-2,ct,data)));

								}

								e[2] = min(e[2],(v.f(i+1,k-1)+erg4(k-1,i+1,i,2,ct,data,lfce[i])+penalty(k-1,i+1,ct,data) + w3[k])
									  );

								if((mod[i+1]||mod[k-1])&&inc[ct->numseq[i+2]][ct->numseq[k-2]]&&!(fce.f(i+1,k-1)&SINGLE)) {

									e[2] = min(e[2],(v.f(i+2,k-2)+erg4(k-1,i+1,i,2,ct,data,lfce[i])+
												penalty(k-1,i+1,ct,data) + w3[k]+erg1(i+1,k-1,i+2,k-2,ct,data)));

								}

								e[3] = min(e[3],(v.f(i,k-2)+erg4(k-2,i,k-1,1,ct,data,lfce[k-1]) + penalty(k-2,i,ct,data) + w3[k])
									  );

								if((mod[i]||mod[k-2])&&inc[ct->numseq[i+1]][ct->numseq[k-3]]&&!(fce.f(i,k-2)&SINGLE)) {
									e[3] = min(e[3],(v.f(i+1,k-3)+erg4(k-2,i,k-1,1,ct,data,lfce[k-1]) + 
												penalty(k-2,i,ct,data) + w3[k]+erg1(i,k-2,i+1,k-3,ct,data)));

								}

								if (!lfce[i]&&!lfce[k-1]) {
									e[4] = min(e[4],(v.f(i+1,k-2)+data->tstack[ct->numseq[k-2]][ct->numseq[i+1]]
												[ct->numseq[k-1]][ct->numseq[i]]
												+checknp(lfce[k-1],lfce[i])+w3[k]+
												penalty(k-2,i+1,ct,data))+ ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(k-1));



									if((mod[i+1]||mod[k-2])&&inc[ct->numseq[i+2]][ct->numseq[k-3]]&&!(fce.f(i+1,k-2)&SINGLE)) {
										e[4] = min(e[4],(v.f(i+2,k-3)+data->tstack[ct->numseq[k-2]][ct->numseq[i+1]]
													[ct->numseq[k-1]][ct->numseq[i]]
													+checknp(lfce[k-1],lfce[i])+w3[k]+
													penalty(k-2,i+1,ct,data)+erg1(i+1,k-2,i+2,k-3,ct,data))
												+ ct->SHAPEss_give_value(i)+ ct->SHAPEss_give_value(k-1));


									}
								}

								//also consider coaxial stacking:
#ifndef disablecoax
								for (ip=k+minloop+1;ip<=number+1;ip++) {


									//first consider flush stacking:
									castack=min(castack,v.f(i,k-1)+v.f(k,ip-1)+w3[ip]+
											penalty(i,k-1,ct,data)+penalty(k,ip-1,ct,data)+
											ergcoaxflushbases(i,k-1,k,ip-1,ct,data));

									if(mod[i]||mod[k-1]||mod[k]||mod[ip-1]) {

										if ((mod[i]||mod[k-1])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[i+1]][ct->numseq[k-2]]
												&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&!(fce.f(i,k-1)&SINGLE)&&!(fce.f(k,ip-1)&SINGLE)) {
											castack=min(castack,v.f(i+1,k-2)+v.f(k+1,ip-2)+w3[ip]+
													penalty(i,k-1,ct,data)+penalty(k,ip-1,ct,data)+
													ergcoaxflushbases(i,k-1,k,ip-1,ct,data)
													+erg1(i,k-1,i+1,k-2,ct,data)+erg1(k,ip-1,k+1,ip-2,ct,data));
										}
										if((mod[i]||mod[k-1])&&inc[ct->numseq[i+1]][ct->numseq[k-2]]&&!(fce.f(i,k-1)&SINGLE)) {

											castack=min(castack,v.f(i+1,k-2)+v.f(k,ip-1)+w3[ip]+
													penalty(i,k-1,ct,data)+penalty(k,ip-1,ct,data)+
													ergcoaxflushbases(i,k-1,k,ip-1,ct,data)
													+erg1(i,k-1,i+1,k-2,ct,data));

										}

										if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&!(fce.f(k,ip-1)&SINGLE)) {

											castack=min(castack,v.f(i,k-1)+v.f(k+1,ip-2)+w3[ip]+
													penalty(i,k-1,ct,data)+penalty(k,ip-1,ct,data)+
													ergcoaxflushbases(i,k-1,k,ip-1,ct,data)
													+erg1(k,ip-1,k+1,ip-2,ct,data));
										}

									}


									//now consider an intervening mismatch:
									if (i>0&&!lfce[i]&&!lfce[k-1]) {
										castack=min(castack,v.f(i+1,k-2)+v.f(k,ip-1)+w3[ip]+
												penalty(i+1,k-2,ct,data)+penalty(k,ip-1,ct,data)+
												ergcoaxinterbases1(i+1,k-2,k,ip-1,ct,data));

										if(mod[i+1]||mod[k-2]||mod[k]||mod[ip-1]){

											if((mod[i+1]||mod[k-2])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[i+2]][ct->numseq[k-3]]
													&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&!(fce.f(i+1,k-2)&SINGLE)&&!(fce.f(k,ip-1)&SINGLE)){
												castack=min(castack,v.f(i+2,k-3)+v.f(k+1,ip-2)+w3[ip]+
														penalty(i+1,k-2,ct,data)+penalty(k,ip-1,ct,data)+
														ergcoaxinterbases1(i+1,k-2,k,ip-1,ct,data)
														+erg1(i+1,k-2,i+2,k-3,ct,data)+erg1(k,ip-1,k+1,ip-2,ct,data)
													   );

											}

											if((mod[i+1]||mod[k-2])&&inc[ct->numseq[i+2]][ct->numseq[k-3]]&&!(fce.f(i+1,k-2)&SINGLE)) {
												castack=min(castack,v.f(i+2,k-3)+v.f(k,ip-1)+w3[ip]+
														penalty(i+1,k-2,ct,data)+penalty(k,ip-1,ct,data)+
														ergcoaxinterbases1(i+1,k-2,k,ip-1,ct,data)
														+erg1(i+1,k-2,i+2,k-3,ct,data)
													   );

											}
											if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&!(fce.f(k,ip-1)&SINGLE)) {
												castack=min(castack,v.f(i+1,k-2)+v.f(k+1,ip-2)+w3[ip]+
														penalty(i+1,k-2,ct,data)+penalty(k,ip-1,ct,data)+
														ergcoaxinterbases1(i+1,k-2,k,ip-1,ct,data)
														+erg1(k,ip-1,k+1,ip-2,ct,data)
													   );	


											}


										}

									}
									if (!lfce[k-1]&&!lfce[ip-1]) {

										castack = min(castack,v.f(i,k-2)+v.f(k,ip-2)+w3[ip]+
												penalty(i,k-2,ct,data)+penalty(k,ip-2,ct,data)+
												ergcoaxinterbases2(i,k-2,k,ip-2,ct,data));

										if (mod[i]||mod[k-2]||mod[k]||mod[ip-2]) {

											if ((mod[i]||mod[k-2])&&(mod[k]||mod[ip-2])&&inc[ct->numseq[i+1]][ct->numseq[k-3]]
													&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&!(fce.f(i,k-2)&SINGLE)&&!(fce.f(k,ip-2)&SINGLE)) {

												castack = min(castack,v.f(i+1,k-3)+v.f(k+1,ip-3)+w3[ip]+
														penalty(i,k-2,ct,data)+penalty(k,ip-2,ct,data)+
														ergcoaxinterbases2(i,k-2,k,ip-2,ct,data)
														+erg1(i,k-2,i+1,k-3,ct,data)+erg1(k,ip-2,k+1,ip-3,ct,data)
													     );
											}

											if ((mod[i]||mod[k-2])&&inc[ct->numseq[i+1]][ct->numseq[k-3]]&&!(fce.f(i,k-2)&SINGLE)) {

												castack = min(castack,v.f(i+1,k-3)+v.f(k,ip-2)+w3[ip]+
														penalty(i,k-2,ct,data)+penalty(k,ip-2,ct,data)+
														ergcoaxinterbases2(i,k-2,k,ip-2,ct,data)
														+erg1(i,k-2,i+1,k-3,ct,data)
													     );
											}

											if ((mod[k]||mod[ip-2])&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&!(fce.f(k,ip-2)&SINGLE)) {

												castack = min(castack,v.f(i,k-2)+v.f(k+1,ip-3)+w3[ip]+
														penalty(i,k-2,ct,data)+penalty(k,ip-2,ct,data)+
														ergcoaxinterbases2(i,k-2,k,ip-2,ct,data)
														+erg1(k,ip-2,k+1,ip-3,ct,data)
													     );
											}

										}
									}

								}
#endif //ifndef disablecoax

							}
							w3[i] = min(w3[i],e[1]);
							w3[i] = min(w3[i],e[2]);
							w3[i] = min(w3[i],e[3]);
							w3[i] = min(w3[i],e[4]);
#ifndef disablecoax
							w3[i] = min(w3[i],castack);
#endif //ifndef disablecoax
						}



					}



#ifdef DYNALIGN_II

			if(j<=number){
				for(size_t it=1;it<=5;it++)e[it]=0;
						
			castack=0;

			//	we.f(i,j)=INFINITE_ENERGY;
			if(j-i+1<=minloop+1){
			  if(lfce[j]||lfce[i])we->f(i,j)=INFINITE_ENERGY;
			  else we->f(i,j)=min(we->f(i,j-1),we->f(i+1,j));
			  
			}
			 
			else{
			  if(lfce[j]||lfce[i])we->f(i,j)=INFINITE_ENERGY;
			  else we->f(i,j)=min(we->f(i,j-1),we->f(i+1,j));

			  for(k=1;k<=5;k++)e[k]=INFINITE_ENERGY;
			  
			  castack=INFINITE_ENERGY;

			  e[1]=min(e[1],v.f(i,j)+penalty(j,i,ct,data));
			  
			  if (mod[i]||mod[j]&&inc[ct->numseq[i+1]][ct->numseq[j-1]]&&!(fce.f(i,j)&SINGLE)) {
			    e[1] = min(e[1],(v.f(i+1,j-1)+penalty(j,i,ct,data)
					     +erg1(i,j,i+1,j-1,ct,data)));
			  }
			  
			  e[2]=min(e[2],(v.f(i+1,j)+penalty(j,i+1,ct,data)+erg4(j,i+1,i,2,ct,data,lfce[i])));//yinghan:erg4, how to set up the fourth parameter? sep.27
			  
			  if (mod[i+1]||mod[j]&&inc[ct->numseq[i+2]][ct->numseq[j-1]]&&!(fce.f(i+1,j)&SINGLE)) {
			    e[2] = min(e[2],(v.f(i+2,j-1)+penalty(j,i+1,ct,data)+erg4(j,i+1,i,2,ct,data,lfce[i])
					     +erg1(i+1,j,i+2,j-1,ct,data)));
			  }//yinghan:shoudln't they check mod first? sep.27

			  e[3]=min(e[3],(v.f(i,j-1)+penalty(j-1,i,ct,data)+erg4(j-1,i,j,1,ct,data,lfce[j])));

			  if (mod[i]||mod[j-1]&&inc[ct->numseq[i+1]][ct->numseq[j-2]]&&!(fce.f(i,j-1)&SINGLE)){
			    e[3] = min(e[3],(v.f(i+1,j-2)+penalty(j-1,i,ct,data)+erg4(j-1,i,j,1,ct,data,lfce[j])
					     +erg1(i,j-1,i+1,j-2,ct,data)));
			  }

			  e[4]=min(e[4],(data->tstack[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
					 +checknp(lfce[j],lfce[i])+v.f(i+1,j-1)+penalty(j-1,i+1,ct,data)));

			  if(mod[i+1]||mod[j-1]&&inc[ct->numseq[i+2]][ct->numseq[j-2]]&&!(fce.f(i+1,j-1)&SINGLE)){
			    e[4]=min(e[4],(data->tstack[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
					   +checknp(lfce[j],lfce[i])+v.f(i+2,j-2)+erg1(i+1,j-1,i+2,j-2,ct,data)
					   +penalty(j-1,i+1,ct,data)));
			  }
			    

			  for (k=i;k<=j;k++)e[5]=min(e[5],we->f(i,k)+we->f(k+1,j));

			  
			  castack=min(castack,wca[i][j]);
			}

			we->f(i,j)=min(we->f(i,j),e[1]);
			we->f(i,j)=min(we->f(i,j),e[2]);
			we->f(i,j)=min(we->f(i,j),e[3]);
			we->f(i,j)=min(we->f(i,j),e[4]);
			we->f(i,j)=min(we->f(i,j),e[5]);
			we->f(i,j)=min(we->f(i,j),castack);
			
			}
#else
#endif
			


				}
#ifndef SMP
				if (!ct->intermolecular) {
					if (d>(maxj>number?8:11))
					{
						tempE=curE;
						curE=prevE;
						prevE=tempE;
					}	  
					if(d>(maxj>number?7:10))
						for (int localdp=1;localdp<=d-1;localdp++) {
							for (int locali=((h<=(number-2))?1:(2*number-h-1));locali<=((h<=(number-2))?(number-h-1):number);locali++)
							{
								if (locali<number)
									curE[localdp][locali]=curE[localdp][locali+1];
							}
						}
				}
#endif //ifndef SMP
			}



			//////////////////////////
			//output V, W, WMB, and W2V:
#if defined (debugmode)
			ofstream foo;
			foo.open("arrays.out");
			foo << "i" << "\t"<<"j"<<"\t"<<"v.f(i,j)"<<"\t"<<"w.f(i,j)"<<"\t"<<"wmb.f(i,j)"<<"\t"<<"v.f(j,i+number)"<<"\t"<<"w.f(j,i+number)"<<"\t"<<"wmb.f(j,i+number)"<<"\n";
			for (j=0;j<=number;j++) {
				for (i=0;i<=j;i++) {

					foo << i << "\t"<<j<<"\t"<<v.f(i,j)<<"\t"<<w.f(i,j)<<"\t"<<wmb.f(i,j)<<"\t"<<v.f(j,i+number)<<"\t"<<w.f(j,i+number)<<"\t"<<wmb.f(j,i+number)<<"\n";

				}	
			}

			foo <<"\n\n\n";
			foo << "i" << "\t" << "w5[i]" << "\t" << "w3[i]" << "\n";
			for (i=0;i<=number;i++) foo << i << "\t" << w5[i] << "\t" << w3[i] << "\n";

			foo.close();

#endif


			//clean up memory use:
			for (int locali=0;locali<=number;locali++)
				delete[] wca[locali];
			delete[] wca;

			if (!ct->intermolecular) {
				for (int locali=0;locali<=number;locali++) {
					delete[] curE[locali];
					delete[] prevE[locali];
				}
				delete[] curE;
				delete[] prevE;

			}


		}//end fill

#ifndef INSTRUMENTED//IF INSTRUMENTED IS NOT DEFINED


	//this function is used to set up the fce array for a base, x, that should be single-stranded
	void forcesingle(int x,structure* ct,forceclass *v) {
		int i;

		for (i=x;i<x+(ct->GetSequenceLength());i++) {
			v->f(x,i)=v->f(x,i)|SINGLE;
		}
		for (i=1;i<=x;i++) {
			v->f(i,x)=v->f(i,x)|SINGLE;
		}
		for (i=x+1;i<=ct->GetSequenceLength();i++) {
			v->f(i,x+ct->GetSequenceLength()) = v->f(i,x+ct->GetSequenceLength())|SINGLE;
		}
	}


	void forcepair(int x,int y,structure *ct,forceclass *v) {
		int i,j;
		v->f(x,y) = v->f(x,y)|PAIR;
		v->f(y,x+ct->GetSequenceLength())=v->f(y,x+ct->GetSequenceLength())|PAIR;
		for (i=y+1;i<=x-1+ct->GetSequenceLength();i++) {
			v->f(x,i) = v->f(x,i)|NOPAIR;
		}
		for (i=x;i<=y-1;i++) {
			v->f(x,i) = v->f(x,i)|NOPAIR;
		}
		for (i=1;i<=x-1;i++) {
			v->f(i,y) = v->f(i,y)|NOPAIR;
		}
		for (i=x+1;i<=y;i++) {
			v->f(i,y) = v->f(i,y)|NOPAIR;
		}
		for (i=1;i<=x-1;i++) {
			v->f(i,x) = v->f(i,x)|NOPAIR;
		}
		for (i=y+1;i<=ct->GetSequenceLength();i++) {
			v->f(i,y+ct->GetSequenceLength())=v->f(i,y+ct->GetSequenceLength())|NOPAIR;
		}
		for (i=y;i<=x-1+(ct->GetSequenceLength());i++) {
			v->f(y,i) = v->f(y,i)|NOPAIR;
		}
		for (i=(ct->GetSequenceLength())+x+1;i<=(ct->GetSequenceLength())+y-1;i++) {
			v->f(y,i) = v->f(y,i)|NOPAIR;
		}
		for (i=x+1;i<=y-1;i++) {
			v->f(i,x+ct->GetSequenceLength()) = v->f(i,x+ct->GetSequenceLength())|NOPAIR;
		}
		for (i=y+1;i<=ct->GetSequenceLength();i++) {
			v->f(i,x+ct->GetSequenceLength()) = v->f(i,x+ct->GetSequenceLength())|NOPAIR;
		}
		for (i=1;i<=x-1;i++) {
			for (j = x+1;j<=y-1;j++){
				v->f(i,j) = v->f(i,j)|NOPAIR;
			}
		}
		for (i=x+1;i<=y-1;i++) {
			for (j=y+1;j<=(ct->GetSequenceLength())+x-1;j++) {
				v->f(i,j) = v->f(i,j)|NOPAIR;
			}
		}
		for (i=y+1;i<=ct->GetSequenceLength();i++) {
			for (j=(ct->GetSequenceLength())+x+1;j<=(ct->GetSequenceLength())+y-1;j++) {
				v->f(i,j) = v->f(i,j)|NOPAIR;
			}
		}
	}



	void forcedbl(int dbl,structure* ct,forceclass *w,bool *v) {
		int i,j;


		v[dbl] = true;
		v[dbl+ct->GetSequenceLength()] = true;


		for(i=dbl+1;i<=ct->GetSequenceLength();i++) {
			for (j=1;j<dbl;j++) {
				w->f(j,i) = w->f(j,i)|DUBLE;
			}
		}
		for(j=(dbl+(ct->GetSequenceLength())-1);j>ct->GetSequenceLength();j--) {
			for (i=dbl+1;i<=ct->GetSequenceLength();i++) {
				w->f(i,j) = w->f(i,j)|DUBLE;
			}
		}


	}

	void forceinter(int dbl,structure* ct,forceclass *w) {
		int i,j;


		for(i=dbl+1;i<=ct->GetSequenceLength();i++) {
			for (j=1;j<dbl;j++) {
				w->f(j,i) = w->f(j,i)|INTER;
			}
		}
		for(j=(dbl+(ct->GetSequenceLength())-1);j>ct->GetSequenceLength();j--) {
			for (i=dbl+1;i<=ct->GetSequenceLength();i++) {
				w->f(i,j) = w->f(i,j)|INTER;
			}
		}
		for(i=dbl+1+ct->GetSequenceLength();i<=2*ct->GetSequenceLength();i++) {
			for (j=ct->GetSequenceLength();j<dbl+ct->GetSequenceLength();j++) {
				w->f(j,i) = w->f(j,i)|INTER;
			}
		}





	}

	void forceinterefn(int dbl,structure* ct,forceclass *w) {
		int i,j;


		for(i=dbl+1;i<=ct->GetSequenceLength();i++) {
			for (j=1;j<dbl;j++) {
				w->f(j,i) = w->f(j,i)|INTER;
			}
		}
	}

void filter(structure* ct, int percent, int max, int window) {

	//structure temp;
	short int i,j,k1,k2,crit;
	bool** mark;
	//bool keep;

	mark = new bool *[ct->GetSequenceLength()+1];
	for (i=0;i<=ct->GetSequenceLength();i++) {
		mark[i] = new bool [ct->GetSequenceLength() + 1];
	}
	for (i=1;i<=ct->GetSequenceLength();i++) {
		for (j=i;j<=ct->GetSequenceLength();j++) {
			mark[i][j] = false;
		}
	}

 


	crit = (short) (ct->GetEnergy(1) + abs((int)((float)ct->GetEnergy(1) *((((float)percent)/100.0)))));
	//temp.numofstructures = 0;

	//check each structure:
	for (i=1;i<=ct->GetNumberofStructures();++i) {
		if (ct->GetEnergy(i) > crit) {
			//none of the remaining structures should be kept bcs the free
			
			for (j=ct->GetNumberofStructures();j>=i;j--) ct->RemoveLastStructure();

			//clean up memory:
			de_allocate(mark,ct->GetSequenceLength()+1);
			return; 
			//energy is > than % sort
		}
		else if (i>max) {
			//none of the remaining structures should be kept bcs the max #
			//of structures has been reached

			for (j=ct->GetNumberofStructures();j>=i;j--) ct->RemoveLastStructure();

			de_allocate(mark,ct->GetSequenceLength()+1);
			return;
		}

		//now map the baspairs within size < window and verify whether this structure
		//should be kept or discarded
//<<<<<<< algorithm.cpp
		//keep = false;
		int newpairs=0;

		for (j=1;j<=ct->GetSequenceLength();j++) {
			if (ct->GetPair(j,i)>j) {
				if (!mark[j][ct->GetPair(j,i)]) {
/*||||||| ../../RNAstructure-before-fork/src/algorithm.cpp
		keep = false;
		for (j=1;j<=ct->numofbases;j++) {
			if (ct->basepr[i][j]>j) {
				if (!mark[j][ct->basepr[i][j]]) {
=======
		keep = false;
		for (j=1;j<=ct->GetSequenceLength();j++) {
			if (ct->GetPair(j,i)>j) {
				if (!mark[j][ct->GetPair(j,i)]) {
>>>>>>> ../../RNAstructure-object/src/algorithm.cpp*/
					//this base has not been marked so keep the structure:
					++newpairs;
				}
			}
		}
		for (j=1;j<=ct->GetSequenceLength();j++) {
			if (ct->GetPair(j,i)>j) {
				//now mark the basepairs:
				for (k1=j-window;k1<=j+window;k1++) {
					for (k2=ct->GetPair(j,i)-window;k2<=ct->GetPair(j,i)+window;k2++) {
						if ((k1>0)&&(k2>0)&&(k1<=ct->GetSequenceLength())&&(k2<=ct->GetSequenceLength())) {
							mark[k1][k2] = true;
						}
					}
				}
			}
		}

/*//<<<<<<< algorithm.cpp
		if (newpairs>window) { //structure needs to be kept, copy it over to temp
			(ct->numofstructures)++;
			ct->energy[ct->numofstructures] = ct->energy[i];
			for (j=1;j<=ct->numofbases;j++) {
				ct->basepr[ct->numofstructures][j] = ct->basepr[i][j];
			}
||||||| ../../RNAstructure-before-fork/src/algorithm.cpp
		if (keep) { //structure needs to be kept, copy it over to temp
			(ct->numofstructures)++;
			ct->energy[ct->numofstructures] = ct->energy[i];
			for (j=1;j<=ct->numofbases;j++) {
				ct->basepr[ct->numofstructures][j] = ct->basepr[i][j];
			}
=======*/
		if (newpairs<=window) { //structure needs to be kept, copy it over to temp
			//Get rid of this structure:
			ct->RemoveStructure(i);
			--i;//decrement the counter so that the new ith structure is not missed
//>>>>>>> ../../RNAstructure-object/src/algorithm.cpp*/
		}
	}



	//clean up memory use:
	de_allocate(mark,ct->GetSequenceLength()+1);
}

	void cctout( structure *ct, char *filename) {
	int i, j;
	ofstream out(filename);
    
    
	out << "-100\n";
	out << ct->GetSequenceLength()<<"\n";
	out << ct->GetNumberofStructures()<<" ";

	out << ct->GetCtLabel(1).c_str();

	for (i=1;i<=ct->GetSequenceLength();i++) {
		out << ct->numseq[i]<<"\n";
	}
	for (i=1;i<=ct->GetNumberofStructures();i++) {
		out << ct->GetEnergy(1)<<"\n";
		for (j=1;j<=ct->GetSequenceLength();j++) {
			out << ct->GetPair(j,i)<<"\n";
		}
	}



}


void calcpnum(dotarray *dots, int *pnum, int increment, short int numofbases,
              TProgressDialog *PD) {
	short int i,j;


	for (i=1;i<=numofbases;i++) {
		pnum[i] = 0;
		//count the dots in the ith column
		for (j=i+1;j<=numofbases;j++) {
			if (dots->dot(i,j)<=increment) pnum[i]++;
		}
		//count the dots in the ith row
		for (j=1;j<i;j++) {
			if (dots->dot(j,i)<=increment) pnum[i]++;

		}

	}



}


//This function will give an energy breakdown for a structure, int n,
//stored in ct, it takes the arrays below as input

void energydump (structure *ct, arrayclass *v, datatable *data, int n,char *filename,int ii, int ji) {
	short int stack[500],stackpos,i,j,k,temp,count,helix,stacke[500],stackpose,start;
	ofstream out;
	char number[15];
	bool exterior;

   
	
	out.open(filename);
	stackpos = 0;

   
	sprintf(number,"%f",(float (ct->GetEnergy(n)))/conversionfactor);
	out << "Structure:  "<<n<<"\n";
	out <<"\n# "<<n<<"  Total Energy = "<<number<<"\n\n";
	out << "pair of origin: "<<ii<<" - "<<ji<<"\n\n";
   
 
	stackpos=0;

	//Analyze the exterior loop
	i = 0;
	temp=0;
	while (i<ct->GetSequenceLength()) {
		i++;
		if (ct->GetPair(i,n)>0) {
       	
			if (i<=ii&&ct->GetPair(i,n)>=ji) {
				temp = temp + v->f(ct->GetPair(i,n),i+ct->GetSequenceLength());
			  
			}
			else temp = temp-v->f(i,ct->GetPair(i,n));
		  
			//stackpos++;
			//stack[stackpos] = i;
			
			i = ct->GetPair(i,n);
		}
	}

	sprintf(number,"%f",(float (temp))/conversionfactor);
	out << "Exterior loop energy = "<<number<<"\n";



   

	//Place the forced pair on the stack
	stackpose =1;
	stacke[1] = ii;
	//Trace out the exterior fragment
	while (stackpose>0) {
		i = stacke[stackpose];
		stackpose--;

		helix = 0;
		//follow the helix:
		while (ct->GetPair(i-1,n)==ct->GetPair(i,n)+1) {
			//helix continues:
			temp = erg1(ct->GetPair(i,n),i,ct->GetPair(i,n)+1,i-1,ct,data);
			//note: v->f(i,ct->basepr[n][i])-v->f(i+1,ct->basepr[n][i+1]);
			//does not work here for cases when a loop is open in min free energy structure
			//but closed in the sub optimal structure being traced
			sprintf(number,"%f",(float (temp))/conversionfactor);
			helix = helix + temp;
			out << "\tStack energy = "<<number<<"  for "<<ct->GetPair(i,n)<<"-"
			    <<i<<" onto "<<ct->GetPair(i,n)+1<<"-"<<i-1<<"\n";
			i--;
		}
		sprintf(number,"%f",(float (helix))/conversionfactor);
		out << "Helix energy = "<<number<<"\n";

      


		//now we've come to a loop, what type?
		j = ct->GetPair(i,n);

	
		start = i;
		exterior=false;
		//check for exterior loop
      
	  
		


		temp = v->f(j,i+ct->GetSequenceLength());
		
		while (i!=j+1&&i!=j) {
			i--;
			if (i==0) {
				exterior=true;
				i = ct->GetSequenceLength()+1;
			 
			}
			else if (ct->GetPair(i,n)>0) {
				//we've found another helix
			
				k = ct->GetPair(i,n);
				
            
				if (i<=start&&k<i) {
					stackpos++;
					stack[stackpos] = ct->GetPair(i,n);
					temp = temp - v->f(k,i);
				}
				else if (i>start) {
					stackpos++;
					stack[stackpos] = ct->GetPair(i,n);
					temp = temp - v->f(k,i);
				} 
				else {
					stackpose++;
					stacke[stackpose] = i;
					temp = temp - v->f(k,i+ct->GetSequenceLength());

				}
				i = k;

			}
		 

		}


		if (!exterior) {

		
		  

     
			//multi loop
			sprintf(number,"%f",(float (temp))/conversionfactor);
			out << "Multibranch loop energy = "<<number<<"  for closure by "<<
				ct->GetPair(j,n)<<"-"<<j<<"\n";

		}





	}
   

	//Place the forced pair on the stack
	stackpos++;
	stack[stackpos] = ii;


	//Trace out the interior fragment
	while (stackpos>0) {
		i = stack[stackpos];
		stackpos--;

		helix = 0;
		//follow the helix:
		while (ct->GetPair(i+1,n)==ct->GetPair(i,n)-1) {
			//helix continues:
			temp = erg1(i,ct->GetPair(i,n),i+1,ct->GetPair(i+1,n),ct,data);
			//note: v->f(i,ct->basepr[n][i])-v->f(i+1,ct->basepr[n][i+1]);
			//does not work here for cases when a loop is open in min free energy structure
			//but closed in the sub optimal structure being traced
			sprintf(number,"%f",(float (temp))/conversionfactor);
			helix = helix + temp;
			out << "\tStack energy = "<<number<<"  for "<<(i+1)<<"-"
			    <<ct->GetPair(i+1,n)<<" onto "<<i<<"-"<<ct->GetPair(i,n)<<"\n";
			i++;
		}
		sprintf(number,"%f",(float (helix))/conversionfactor);
		out << "Helix energy = "<<number<<"\n";

      


		//now we've come to a loop, what type?
		j = ct->GetPair(i,n);

	

		temp = v->f(i,j);
		count = 0;
		while (i<j-1) {
			i++;
			if (ct->GetPair(i,n)>0) {
				//we've found another helix
				count++;
				k = ct->GetPair(i,n);
				temp = temp - v->f(i,k);
				stackpos++;
				stack[stackpos] = i;
				i = k;

			}

		}
		if (count ==0) {
			//hairpin:
			sprintf(number,"%f",(float (temp))/conversionfactor);
			out << "Hairpin energy = "<<number<<"  for closure by "<<
				ct->GetPair(j,n)<<"-"<<j<<"\n";

		}
		else if (count==1) {
			//bulge or internal loop:
			sprintf(number,"%f",(float (temp))/conversionfactor);
			out << "Bulge/Internal loop energy = "<<number<<"  for closure by "<<
				ct->GetPair(j,n)<<"-"<<j<<"\n";

		}
		else {
			//multi loop
			sprintf(number,"%f",(float (temp))/conversionfactor);
			out << "Multibranch loop energy = "<<number<<"  for closure by "<<
				ct->GetPair(j,n)<<"-"<<j<<"\n";

		}





	}

	
   
	out.close();

}


//This function will give an energy breakdown for a structure, int n,
//stored in ct, it takes the arrays below as input
//note:this is not fully debugged
void energydump (structure *ct, datatable *data,arrayclass *v, int n,char *filename) {
	int stack[500],stackpos,i,j,k,temp,count,helix,auaddition;
	ofstream out;
	char number[6],auend[6];
	int inc[6][6] = {{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,0,0,0},{0,0,0,0,1,0}
	                 ,{0,1,0,1,0,0},{0,0,0,0,0,0}};
	bool sbulge;

	sbulge = false;
	out.open(filename);
	stackpos = 0;

  
	sprintf(auend,"%f",(float (data->auend))/conversionfactor);

   
	sprintf(number,"%f",(float (ct->GetEnergy(n)))/conversionfactor);
	out << "Structure:  "<<n<<"\n";
	out <<"\n# "<<n<<"  Total Energy = "<<number<<"\n\n";

	temp = ct->GetEnergy(n);

	//Analyze the exterior loop
	i = 0;
	while (i<ct->GetSequenceLength()) {
		i++;
		if (ct->GetPair(i,n)>0) {
			stackpos++;
			stack[stackpos] = i;
			if (inc[ct->numseq[i]][ct->numseq[ct->GetPair(i,n)]]) {
				temp = temp - data->auend;
			}
			temp = temp - v->f(i,ct->GetPair(i,n));
			i =ct->GetPair(i,n);
		}
	}
	sprintf(number,"%f",(float (temp))/conversionfactor);
	out << "Exterior loop energy = "<<number<<"\n";

	while (stackpos>0) {
		i = stack[stackpos];
		stackpos--;
		helix = 0;

		if (inc[ct->numseq[i]][ct->numseq[ct->GetPair(i,n)]]&&!sbulge) {
			out << "Non-GC end = "<<auend<<"\n";
			helix = helix + data->auend;
		}
		sbulge = false;

		//follow the helix:
		while (ct->GetPair(i+1,n)==ct->GetPair(i,n)-1) {
			//helix continues:
			temp = erg1(i,ct->GetPair(i,n),i+1,ct->GetPair(i+1,n),ct,data);
			//note: v->f(i,ct->basepr[n][i])-v->f(i+1,ct->basepr[n][i+1]);
			//does not work here for cases when a loop is open in min free energy structure
			//but closed in the sub optimal structure being traced
			sprintf(number,"%f",(float (temp))/conversionfactor);
			helix = helix + temp;
			out << "Stack energy = "<<number<<"  for "<<(i+1)<<"-"
			    <<ct->GetPair(i+1,n)<<" onto "<<i<<"-"<<ct->GetPair(i,n)<<"\n";
			i++;
		}




		//now we've come to a loop, what type?
		auaddition = 0;
		j = ct->GetPair(i,n);
		temp = v->f(i,j);
		count = 0;
		while (i<j-1) {
			i++;
			if (ct->GetPair(i,n)>0) {
				//we've found another helix
				if (inc[ct->numseq[i]][ct->numseq[ct->GetPair(i,n)]])
					auaddition = auaddition - data->auend;
				count++;
				k = ct->GetPair(i,n);
				temp = temp - v->f(i,k);
				stackpos++;
				stack[stackpos] = i;
				i = k+1;

			}

		}
		//check for single bulge:
		if (count==1) {
			if (ct->GetPair(j-1,n)>0) {
				if (ct->GetPair(ct->GetPair(j,n)+1,n)==0) {
					if (ct->GetPair(ct->GetPair(j,n)+2,n)>0) {
						sbulge = true;
						auaddition = 0;
					}
				}

			}
			else if (ct->GetPair(ct->GetPair(j,n)+1,n)>0) {
				if (ct->GetPair(j-1,n)==0) {
					if (ct->GetPair(j-2,n)>0) {
						sbulge = true;
						auaddition = 0;
					}
				}
			}

		}

		if (inc[ct->numseq[j]][ct->numseq[ct->GetPair(j,n)]]&&!sbulge) {
			out << "Non-GC end = "<<auend<<"\n";
			helix = helix + data->auend;
		}


     
		sprintf(number,"%f",(float (helix))/conversionfactor);
		out << "\tHelix energy = "<<number<<"\n";

		if (count ==0) {
			//hairpin:
			temp = temp + auaddition;
			sprintf(number,"%f",(float (temp))/conversionfactor);
			out << "Hairpin energy = "<<number<<"  for closure by "<<
				ct->GetPair(j,n)<<"-"<<j<<"\n";

		}
		else if (count==1) {
			//bulge or internal loop:
			temp = temp + auaddition;
			sprintf(number,"%f",(float (temp))/conversionfactor);
			out << "Bulge/Internal loop energy = "<<number<<"  for closure by "<<
				ct->GetPair(j,n)<<"-"<<j<<"\n";

		}
		else {
			//multi loop
			temp = temp + auaddition;
			sprintf(number,"%f",(float (temp))/conversionfactor);
			out << "Multibranch loop energy = "<<number<<"  for closure by "<<
				ct->GetPair(j,n)<<"-"<<j<<"\n";

		}





	}
	out.close();





}

//cntrl6 = #tracebacks
//cntrl8 = percent sort
//cntrl9 = window
void readsav(const char *filename, structure *ct, arrayclass *w2, arrayclass *wmb2, 
             integersize *w5, integersize *w3, bool *lfce, bool *mod, datatable *data,
             arrayclass *v, arrayclass *w, arrayclass *wmb, forceclass *fce, int *vmin) {

	int i,j,k,l,m,n,o,p;

	short vers;//save file version

	register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	                        {0,1,0,1,0,0},{0,0,0,0,0,0}};

	ifstream sav(filename,ios::binary);
	
	//read the save file information so that the sequence can be re-folded,
	//include thermodynamic data to prevent traceback errors


	//read the file version first
	read(&sav,&vers);

	//start with structure information
	int sequencelength;
	read(&sav,&(sequencelength));
	read(&sav,&(ct->intermolecular));

	
	int pairnumber;
	read(&sav,&(pairnumber));

	int pair5,pair3;
	for (i=0;i<pairnumber;i++) {
		read(&sav,&(pair5));
		
		read(&sav,&(pair3));
		ct->AddPair(pair5,pair3);
	}
	read(&sav,&(pairnumber));
	for (i=0;i<pairnumber;i++) {
		read(&sav,&(pair5));
		read(&sav,&(pair3));

		ct->AddForbiddenPair(pair5,pair3);
	}
	for (i=0;i<=ct->GetSequenceLength();i++) {
		
		read(&sav,&(ct->hnumber[i]));
		sav.read(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->numseq[i]));
	
	//Read the number of nucleotides forced double stranded
	read(&sav,&(pairnumber));

	//Read and set double stranded nucleotides
	for (i=0;i<pairnumber;i++) {
		int doubles;
		read(&sav,&(doubles));

		ct->AddDouble(doubles);

	}

	
	if (ct->intermolecular) {
		w2 = new arrayclass(ct->GetSequenceLength());
		wmb2 = new arrayclass(ct->GetSequenceLength());
		
		for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));
		
	}

	//Read the number of nucleotides that are not allowed to pair
	read(&sav,&(pairnumber));

	//Read the nucleotides that cannot pair and store the information.
	for (i=0;i<pairnumber;i++) {
		int singles;
		read(&sav,&(singles));

		ct->AddSingle(singles);

	}


	//Read the number of chemically modified nucleotides
	read(&sav,&(pairnumber));

	//Read and store the chemically modified nucleotides
	for (i=0;i<pairnumber;i++) {

		int modified;
		read(&sav,&(modified));

		ct->AddModified(modified);

	}
	
	//Read the number of Us in GU pairs
	read(&sav,&(pairnumber));

	//Read the Us in GU pairs, and store information
	for (i=0;i<pairnumber;i++) {
		
		int GU;
		read(&sav,&(GU));  

		ct->AddGUPair(GU);

	}
	
	string label;
	read(&sav,&(label));
	ct->SetSequenceLabel(label);

	read(&sav,&(ct->templated));
	if (ct->templated) {

		ct->allocatetem();
		for (i=0;i<=ct->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) read(&sav,&(ct->tem[i][j]));	

		}

	}

	read(&sav,&ct->shaped);
	if (ct->shaped) {
		ct->SHAPE = new double [2*ct->GetSequenceLength()+1];
		for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->SHAPE[i]));

	}
	

	//now read the array class data for v, w, and wmb:
	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct->GetSequenceLength();i++) {
		read(&sav,&(w3[i]));
		read(&sav,&(w5[i]));
		for (j=0;j<=ct->GetSequenceLength();j++) {
			read(&sav,&(v->dg[i][j+i]));
			read(&sav,&(w->dg[i][j+i]));
			read(&sav,&(wmb->dg[i][j+i]));
			readsinglechar(&sav,&(fce->dg[i][j]));
			if (ct->intermolecular) {
				read(&sav,&(w2->dg[i][j+i]));
				read(&sav,&(wmb2->dg[i][j+i]));

			}
		}
	}

	read(&sav,&(w3[ct->GetSequenceLength()+1]));
	for (i=0;i<=2*ct->GetSequenceLength();i++) {
		read(&sav,&(lfce[i]));
		read(&sav,&(mod[i]));
	}

	read(&sav, vmin);

	//now open the data files:
	for (i=0;i<5;i++) read(&sav,&(data->poppen[i]));
	read(&sav,&(data->maxpen));
	for (i=0;i<11;i++) read(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		read(&sav,&(data->inter[i]));
		read(&sav,&(data->bulge[i]));
		read(&sav,&(data->hairpin[i]));
	}

	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					read(&sav,&(data->dangle[i][j][k][l]));	
				}
				for (l=0;l<6;l++) {
					read(&sav,&(data->stack[i][j][k][l]));
					read(&sav,&(data->tstkh[i][j][k][l]));
					read(&sav,&(data->tstki[i][j][k][l]));
					read(&sav,&(data->coax[i][j][k][l]));
					read(&sav,&(data->tstackcoax[i][j][k][l]));
					read(&sav,&(data->coaxstack[i][j][k][l]));
					read(&sav,&(data->tstack[i][j][k][l]));
					read(&sav,&(data->tstkm[i][j][k][l]));
					read(&sav,&(data->tstki23[i][j][k][l]));
					read(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							read(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<6;o++) {
								if (inc[i][j]&&inc[n][o]) read(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								else data->iloop21[i][j][k][l][m][n][o]=INFINITE_ENERGY;
								for (p=0;p<6;p++) {
									if (inc[i][k]&&inc[j][l])
										read(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
									else data->iloop22[i][j][k][l][m][n][o][p]=INFINITE_ENERGY;
								}
							}
							

						}
					}
				}
			}
		}
	}
	read(&sav,&(data->numoftloops));
	for (i=0;i<=data->numoftloops;i++) {
		for (j=0;j<2;j++) read(&sav,&(data->tloop[i][j]));

	}
	read(&sav,&(data->numoftriloops));
	for (i=0;i<=data->numoftriloops;i++) {
		for (j=0;j<2;j++) read(&sav,&(data->triloop[i][j]));

	}
	read(&sav,&(data->numofhexaloops));
	for (i=0;i<=data->numofhexaloops;i++) {
		for (j=0;j<2;j++) read(&sav,&(data->hexaloop[i][j]));

	}
	read(&sav,&(data->auend));
	read(&sav,&(data->gubonus));
	read(&sav,&(data->cint));
	read(&sav,&(data->cslope));
	read(&sav,&(data->c3));
	read(&sav,&(data->efn2a));
	read(&sav,&(data->efn2b));
	read(&sav,&(data->efn2c));
	read(&sav,&(data->init));
	read(&sav,&(data->mlasym));
	read(&sav,&(data->strain));
	read(&sav,&(data->prelog));
	read(&sav,&(data->singlecbulge));



	sav.close();


}



//opensav will open a save file created by the fill algorithm (function dynamic)
//	it then runs the traceback routine with the parameters provided
void opensav(char* filename, structure* ct, int cntrl6, int cntrl8,int cntrl9) {
	int i;

	//register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	//{0,1,0,1,0,0},{0,0,0,0,0,0}};
		
	arrayclass *w2,*wmb2;
	integersize *w5,*w3;
	int vmin;
	bool *lfce,*mod;
	datatable *data;

	data = new datatable;
	short vers;



	//peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
	ifstream sav(filename,ios::binary);
	
	read(&sav,&vers);//could add functionality here to make sure that the correct version is being read
	int sequencelength;
	read(&sav,&(sequencelength));
	read(&sav,&(ct->intermolecular));
	sav.close();


	//allocate everything
	ct->allocate(sequencelength);
	
	arrayclass w(ct->GetSequenceLength());
	arrayclass v(ct->GetSequenceLength());
	arrayclass wmb(ct->GetSequenceLength());	
	forceclass fce(ct->GetSequenceLength());


	lfce = new bool [2*ct->GetSequenceLength()+1];
	mod = new bool [2*ct->GetSequenceLength()+1];
	
	w5 = new integersize [ct->GetSequenceLength()+1];
	w3 = new integersize [ct->GetSequenceLength()+2];

	if (ct->intermolecular) {
		w2 = new arrayclass(ct->GetSequenceLength());
		wmb2 = new arrayclass(ct->GetSequenceLength());
		
		for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));
		
	}
	else {
		w2 = NULL;
		wmb2 = NULL;
	}
	
	readsav(filename, ct, w2, wmb2, w5, w3, lfce, mod, data,
	        &v, &w, &wmb, &fce, &vmin);

	traceback(ct, data, &v, &w, &wmb, w2, wmb2,w3, w5, &fce, lfce, vmin, cntrl6, cntrl8, cntrl9,mod);


	delete[] lfce;
	delete[] mod;
	



	delete[] w5;
	delete[] w3;


	if (ct->intermolecular) {
		delete w2;
		delete wmb2;
	}

	delete data;

}

	//returns true if pair i and j is not a GU pair and the most dajacent pairs are not GU
	//	used with chemical modification data to make sure that some contributions are not over-counted
	bool notgu(int i, int j, structure *ct) {

		if (ct->numseq[i]==3&&ct->numseq[j]==4) return false;
		else if (ct->numseq[i]==4&&ct->numseq[j]==3) return false;
		else if (ct->numseq[i+1]==3&&ct->numseq[j-1]==4) return false;
		else if (ct->numseq[i+1]==4&&ct->numseq[j-1]==3) return false;
		else if (ct->numseq[i-1]==3&&ct->numseq[j+1]==4) return false;
		else if (ct->numseq[i-1]==4&&ct->numseq[j+1]==3) return false;
		else return true;

	}

#endif //END INSTRUMENTED
