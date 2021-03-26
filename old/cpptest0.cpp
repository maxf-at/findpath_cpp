//g++ cpptest0.cpp -o cpptest0 -lRNA -lm

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <algorithm> 
#include <iomanip>
#include <cstring>


using namespace std;

extern "C" {
//#include "../fold.h"

#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/fold.h"
//#include "ViennaRNA/cofold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/structures.h"
//#include "ViennaRNA/landscape/findpath.h"


}



int main(int argc, char** argv)
{
	
  /* The RNA sequence */
  std::string sequence = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA";
  
  char* seq = const_cast<char*>(sequence.c_str());
  
  /* allocate memory for MFE structure (length + 1) */
  char  *structure = (char *)vrna_alloc(sequence.length()+1);
  //std::string structure;
  
  /* predict Minmum Free Energy and corresponding secondary structure */
  float mfe = vrna_fold(seq, structure);
  /* print sequence, structure and MFE */
  //printf("%s\n%s [ %6.2f ]\n", seq, structure, mfe);
  
  short * pt = vrna_ptable(structure);
  
  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
  
    // # fc.eval_move_pt(pt,-1,-33)
    // # fc.eval_move_pt(pt,1,33)
    // fc.eval_move(sa,-1,-33)
    // fc.eval_move(sb,1,33)
    // # fc.eval_structure(sa)
    // # fc.eval_structure(sb)
    // # fc.eval_structure_pt(pa)
    // # fc.eval_structure_pt(pb)
  float test_A;
  float test_B;
  float test_C;
  float test_D;
  clock_t begin, end;

  // begin = clock();
  // for (int i=0; i<1000000; i++)
  //   test_A = vrna_eval_structure(fc, structure);  // 0.83
  // end = clock();  

  // begin = clock();
  // for (int i=0; i<1000000; i++)
  //   test_B = vrna_eval_move(fc, structure, -3, -20); // 0.18
  // end = clock();  

  // begin = clock();
  // for (int i=0; i<1000000; i++)
  //   test_C = vrna_eval_structure_pt(fc, pt); // 0.745
  // end = clock();  

  begin = clock();
  for (int i=0; i<1000000; i++)
    test_D = vrna_eval_move_pt(fc, pt, -3, -20); // 0.09
  end = clock();  

  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  

  // float test_D = vrna_eval_structure_pt(fc, pt) / 100.0;  
  // float test_C = vrna_eval_structure(fc, "...((((........)))).(((((((...))))))).....");
  
  cout << seq << "\n" << structure << mfe << " / " << test_D << " / "<< test_A  <<"\n";
  cout << test_B << " / "<< test_C  << " / " << time_spent << "\n";

  // 3,19
  
  /* cleanup memory */
  //free(structure);
  return 0;
}
