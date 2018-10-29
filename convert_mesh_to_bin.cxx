#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <omp.h>
using namespace std;
int main (int argc, char** argv) {
   std::ifstream input_file(argv[1]);
   int npat_, nnod_;
   fstream outFile("out.bin", ios::out | ios::binary);
   input_file >> nnod_ >> npat_;
   outFile.write((char*)&nnod_, sizeof(int));
   outFile.write((char*)&npat_, sizeof(int));    
   size_t npat = (size_t)npat_;
   size_t nnod = (size_t)nnod_;
#if HUGE_NORM
   double sum[3] = {0.0, 0.0, 0.0};
   double sunod[3];
   std::cout<< "allocating sunod array of size " << nnod * 3 << std::endl; 
   double* sunod_arr = new double [nnod * 3];
   std::cout<< "done allocating sunod array" << std::endl; 
   
   std::cout<< "reading sunod array of size " << nnod * 3 << std::endl; 
   for (size_t i = 0; i < nnod; ++i) {            
     input_file >> sunod[0] >> sunod[1] >> sunod[2];                  
     sunod_arr[i*3+0]= sunod[0];
     sunod_arr[i*3+1]= sunod[1];
     sunod_arr[i*3+2]= sunod[2];
     sum[0]+=sunod[0];
     sum[1]+=sunod[1];
     sum[2]+=sunod[2];
   }
   std::cout<< "done reading sunod array" << std::endl;  
   std::cout<< "normalizing sunod array" << std::endl;  
   sum[0]/=nnod; 
   sum[1]/=nnod; 
   sum[2]/=nnod; 
#pragma omp parallel for   
   for (size_t i = 0; i < nnod; ++i) {           
     sunod_arr[i*3+0]-= sum[0];
     sunod_arr[i*3+1]-= sum[1];
     sunod_arr[i*3+2]-= sum[2];
   }
   std::cout<< "done normalizing sunod array" << std::endl; 
   std::cout<< "writing sunod array to file" << std::endl; 
   outFile.write((char*)sunod_arr, sizeof(double)*nnod*3);
   std::cout<< "done writing sunod array to file" << std::endl; 
   delete sunod_arr;
   int nsupan[6];
   std::cout<< "allocating nsupan array of size " << npat * 6 << std::endl; 
   int* nsupan_arr = new int[npat * 6];
   std::cout<< "done allocating nsupan array of size " << npat * 6 << std::endl; 
  
   std::cout<< "reading nsupan array" << std::endl; 
   for (size_t i = 0; i < npat; ++i) {
     input_file >> nsupan_arr[i*6+0] >> nsupan_arr[i*6+1]  >> nsupan_arr[i*6+2] >> nsupan_arr[i*6+3] >> nsupan_arr[i*6+4] >> nsupan_arr[i*6+5];
   }
   std::cout<< "done reading nsupan array" << std::endl; 
   
   std::cout<< "writing nsupan array to file" << std::endl; 
   outFile.write ((char*)nsupan_arr,sizeof(int)*npat*6);
   std::cout<< "done writing nsupan array to file" << std::endl; 
#elif HUGE_NORM_PART
   double sum[3] = {0.0, 0.0, 0.0};
   double sunod[3];
   size_t sunod_arr_size = nnod * 3;
   std::cout<< "allocating sunod array of size " << sunod_arr_size<< std::endl; 
   double* sunod_arr = new double [sunod_arr_size];
   std::cout<< "done allocating sunod array" << std::endl; 
   
   std::cout<< "reading sunod array of size " << sunod_arr_size << std::endl; 
   for (size_t i = 0; i < nnod; ++i) {            
     input_file >> sunod[0] >> sunod[1] >> sunod[2];                  
     sunod_arr[i*3+0]= sunod[0];
     sunod_arr[i*3+1]= sunod[1];
     sunod_arr[i*3+2]= sunod[2];
     sum[0]+=sunod[0];
     sum[1]+=sunod[1];
     sum[2]+=sunod[2];
   }
   std::cout<< "done reading sunod array" << std::endl;  
   std::cout<< "normalizing sunod array" << std::endl;  
   sum[0]/=nnod; 
   sum[1]/=nnod; 
   sum[2]/=nnod; 
#pragma omp parallel for   
   for (size_t i = 0; i < nnod; ++i) {           
     sunod_arr[i*3+0]-= sum[0];
     sunod_arr[i*3+1]-= sum[1];
     sunod_arr[i*3+2]-= sum[2];
   }
   std::cout<< "done normalizing sunod array" << std::endl; 
#if CHUNKED_WRITE
   size_t nsupan_arr_size = 6;
   std::cout<< "allocating nsupan array of size " << nsupan_arr_size << std::endl; 
   int* nsupan_arr = new int[nsupan_arr_size];
   std::cout<< "done allocating nsupan array of size " << nsupan_arr_size << std::endl; 
   std::cout<< "reading and writing mesh" << std::endl; 
   for (size_t i = 0; i < npat; ++i) {
     input_file >> nsupan_arr[0] >> nsupan_arr[5]  >> nsupan_arr[2] >> nsupan_arr[4] >> nsupan_arr[1] >> nsupan_arr[3];
     std::cout << "outputting triangle: " << i << std::endl;
     for(size_t j = 0; j < 6; ++j) {
       nsupan_arr[j]--;
       outFile.write((char*)&(sunod_arr[nsupan_arr[j]*3 + 0]),sizeof(double));
       outFile.write((char*)&(sunod_arr[nsupan_arr[j]*3 + 1]),sizeof(double));
       outFile.write((char*)&(sunod_arr[nsupan_arr[j]*3 + 2]),sizeof(double));
       outFile.flush();
       std::cout<< "output triangular point: " << j << "  " ;
     }
     std::cout<< std::endl;
   }
   std::cout<< "done writing mesh array" << std::endl; 
   input_file.close();   
   outFile.close();
   delete sunod_arr;
   delete nsupan_arr;
#else 
   size_t nsupan_arr_size = npat * 6;
   std::cout<< "allocating nsupan array of size " << nsupan_arr_size << std::endl; 
   int* nsupan_arr = new int[nsupan_arr_size];
   std::cout<< "done allocating nsupan array of size " << nsupan_arr_size << std::endl; 
   std::cout<< "reading nsupan array" << std::endl; 
   for (size_t i = 0; i < npat; ++i) {
     input_file >> nsupan_arr[i*6+0] >> nsupan_arr[i*6+5]  >> nsupan_arr[i*6+2] >> nsupan_arr[i*6+4] >> nsupan_arr[i*6+1] >> nsupan_arr[i*6+3];
     for(size_t j = 0; j < 6; ++j) nsupan_arr[i*6+j]--;
   }
   input_file.close();
   std::cout<< "done reading nsupan array" << std::endl; 
   std::cout<< "filling out new file content " << std::endl;
   double* file_content = new double [npat*6*3];
   for (size_t i = 0; i < npat; ++i) {
     for(size_t j = 0; j < 6; ++j) {
       file_content[i*6*3 + j*3 + 0] = sunod_arr[nsupan_arr[i*6+j]*3 + 0];
       file_content[i*6*3 + j*3 + 1] = sunod_arr[nsupan_arr[i*6+j]*3 + 1];
       file_content[i*6*3 + j*3 + 2] = sunod_arr[nsupan_arr[i*6+j]*3 + 2];
     } 
   }
   std::cout<< "done filling out new file content " << std::endl;
   delete sunod_arr;
   delete nsupan_arr;
   std::cout<< "writing out nsupan array to file" << std::endl; 
   outFile.write ((char*)file_content,sizeof(double)*npat*6*3);
   std::cout<< "done writing nsupan array to file" << std::endl; 
   outFile.close();
   delete file_content;
#endif
#else
   double sunod[3];
   for (size_t i = 0; i < nnod; ++i) {            
     input_file >> sunod[0] >> sunod[1] >> sunod[2];                  
     outFile.write ((char*)sunod, 3*sizeof(double));
   }
   int nsupan[6];
   for (size_t i = 0; i < npat; ++i) {
     input_file >> nsupan[0] >> nsupan[1]  >> nsupan[2] >> nsupan[3] >> nsupan[4] >> nsupan[5];
     outFile.write ((char*)nsupan, 6*sizeof(int));
   }
#endif
  outFile.close();
  return 0;
}
