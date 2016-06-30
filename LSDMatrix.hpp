///LSDMatrix.hpp
///
/// This Matrix class is a simple, fast [Citation needed], implementation
/// of a dynamically sized matrix object with rows and columns [i][j], where the
/// dimensions can be set at run time. It is an experiment for comparison
/// with the TNT::Arrayxx class we often use in the LSD package.
///
/// It is put together from recipe 11.14 in the "C++ Cookbook" and Stroustrup's
/// "The C++ Programming Langauge", and some other stuff I read on internet forums,
/// so I take no real credit for it!
///
/// Please note it is entirely experimental at this stage..


#ifndef LSDMATRIX_HPP
#define LSDMATRIX_HPP


// Require recipe 11.12!

#include <valarray>
#include <numeric>
#include <algorithm>

template<class Value_T>
class LSDMatrix
{
public:
  typedef Value_T value_type;
  typedef LSDMatrix self;
  typedef value_type* iterator;
  typedef const value_type* const_iterator;
  typedef Value_T* row_type;
  typedef stride_iter<value_type*> col_type; // need to implement stride iter!
  typedef const value_type* const_row_type;
  typedef stride_iter<const value_type*> const_col_type;

  // CONSTRUCTORS
  LSDMatrix() : nrows(0), ncols(0, m() { }

  LSDMatrix(int r, int c)
  : nrows(r), ncols(c), m(r*c) { }

  LSDMatrix(const self& x)
  : m(x.m), nrows(x.rows), ncols(x.cols) { }

  template<typename T>
  explicit LSDMatrix(const valarray<T> & x)
  : m(x.size() + 1), nrows(x.size()), ncols(1)
  {
    for(int i=0; i<x.size(); i++)
    {
      m[i] = x[i];
    }
  }

  // Allow construction from matrices of other types
  template<typename T>
  explicit LSDMatrix(const LSDMatrix<T>& x)
  : m(x.size() + 1), nrows(x.rows), ncols(x.cols)
  {
    copy(x.begin(), x.end(), m.begin());
  }

  // Public Functions
  int rows() const { return nrows; }
  int cols() const { return ncols; }
  int size() const { return nrows * ncols; }

  // Element access
  row_type row_begin(int n)
  {
    return &m[n * ncols];
  }

  row_type row_end


};

#endif // LSDMATRIX_HPP

