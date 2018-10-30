#ifndef SELF_METADATA
#define SELF_METADATA
#include "utils.h"
#include "global_data.h"
namespace bemfmm {
  class self_metadata { 
  public:
#if USE_PART || USE_MMAP
    d_vector const& x;
    d_vector const& y;
    d_vector const& z;
#else
    int** const& nsupan;
    double** const& sunod;
#endif
    integral_data const& int_data;
    d_vector_2d const& vcsio;
    d_complex_t const& cjvk; 
    double const* ipolymatrix;
    d_vector const& wwo;
    d_vector const& xxo;
    d_vector const& wws;
    d_vector const& xxs;
    int const& nlqp;    
    
#if USE_PART || USE_MMAP
   self_metadata(d_vector const& _x, d_vector const& _y, d_vector const& _z, integral_data const& _int_data, 
      d_vector_2d const& _vcsio, d_complex_t const& _cjvk, double const* _ipolymatrix, int const& _nlqp) :
      x(_x), y(_y), z(_z), int_data(_int_data), vcsio(_vcsio), cjvk(_cjvk), 
      ipolymatrix(_ipolymatrix), wwo(_int_data.wwo), xxo(_int_data.xxo), wws(_int_data.wws), xxs(_int_data.xxs), nlqp(_nlqp) { }
#else    
    self_metadata(int** const& _nsupan, double** const& _sunod, integral_data const& _int_data, 
      d_vector_2d const& _vcsio, d_complex_t const& _cjvk, double const* _ipolymatrix, int const& _nlqp) :
      nsupan(_nsupan), sunod(_sunod), int_data(_int_data), vcsio(_vcsio), cjvk(_cjvk), 
      ipolymatrix(_ipolymatrix), wwo(_int_data.wwo), xxo(_int_data.xxo), wws(_int_data.wws), xxs(_int_data.xxs), nlqp(_nlqp) { }
#endif



  };
}
#endif
