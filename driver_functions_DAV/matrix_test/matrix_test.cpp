#include <memory>
#include <iostream>
#include <string>
#include <chrono>
#include "../../TNT/tnt.h"

#include <eigen3/Eigen/Core>


class LSDGrid2D
{
    size_t _rows;
    size_t _columns;
    std::unique_ptr<int[]> data;

public:

    LSDGrid2D(size_t rows, size_t columns)
        : _rows{rows}, _columns{columns}, data{std::make_unique<int[]>(rows * columns)} {
    }


    size_t rows() const {
        return _rows;
    }

    size_t columns() const {
        return _columns;
    }

    int * operator[](size_t row) {
        return row * _columns + data.get();
    }

};


template<class T>
class LSDMatrix2D
{
public:
  LSDMatrix2D(unsigned rows, unsigned ncols);

  // For size being zero, we should throw error
  class BadSize { };

  // Law of big 3
  ~LSDMatrix2D();
  LSDMatrix2D(const LSDMatrix2D<T>& m);
  LSDMatrix2D& operator= (const LSDMatrix2D<T>& m);

  // Array access methods to get element with (i,j) notation
  T& operator() (unsigned i, unsigned j);
  const T& operator() (unsigned i, unsigned j) const;

  // Throw bounds violation if i, j too big
  class BoundsViolation { };

private:
  T* data_;
  unsigned nrows_, ncols_;
};

// Access matrix elements with m(i,j) notation
template<class T>
inline T& LSDMatrix2D<T>::operator() (unsigned row, unsigned col)
{
  if (row >= nrows_ || col >= ncols_) throw BoundsViolation();
  return data_[row*ncols_ + col];
}

// Access matrix elements with m(i,j) notation (constatnt)
template<class T>
inline const T& LSDMatrix2D<T>::operator() (unsigned row, unsigned col) const
{
  if (row >= nrows_ || col >= ncols_) throw BoundsViolation();
  return data_[row*ncols_ + col];
}

// Declare matrix
template<class T>
inline LSDMatrix2D<T>::LSDMatrix2D(unsigned nrows, unsigned ncols)
  : data_ (new T[nrows * ncols]),
    nrows_ (nrows),
    ncols_ (ncols)
{
  if (nrows == 0 || ncols == 0)
    throw BadSize();
}

// Clean up after we're done with our matrix!
template<class T>
inline LSDMatrix2D<T>::~LSDMatrix2D()
{
  delete[] data_;
}




class LSDArray2D
{
  int* array;
  int m_width;
public:
  LSDArray2D( int w, int h )
    :m_width( w ), array( new int[ w * h ] ) {}

  ~LSDArray2D() { delete[] array; }

  int at( int x, int y )
  const { return array[ index( x, y ) ]; }

  void set(int x, int y, int v)
  const { array[index (x, y) ] = v ; }

protected:
  int index( int x, int y )
  const { return x + m_width * y; }

};




int main()
{
  size_t imax;
  size_t jmax;

  std::cout << "Enter number of rows: " << std::endl;
  std::cin >> imax;

  std::cout << "Enter number of cols: " << std::endl;
  std::cin >> jmax;

//  std::cout << "LSD Version: " << std::endl;
//  auto start = std::chrono::high_resolution_clock::now();
//  LSDGrid2D raster(imax, jmax);

//  for (int i=0; i<imax; ++i)
//  {
//    for (int j=0; j<jmax; ++j)
//    {
//      raster[i][j] = 4;
//      raster[i][j] *= 7;
//       //std::cout << raster[i][j] << ' ';
//    }
//    //std::cout << std::endl;
//  }
//  auto elapsed = std::chrono::high_resolution_clock::now() - start;
//  long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
//  std::cout << "Time: " << microseconds << std::endl;

  std::cout << "C++FAQ Version: " << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  LSDMatrix2D<double> raster1(imax, jmax);

  for (int i=0; i<imax; ++i)
  {
    for (int j=0; j<jmax; ++j)
    {
      raster1(i, j) = 4;
      raster1(i, j) /= 7;
       //std::cout << raster[i][j] << ' ';
    }
    //std::cout << std::endl;
  }
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  std::cout << "Time: " << microseconds << std::endl;

  std::cout << "Eigen Version: " << std::endl;
  auto start1 = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd rasterEigen(imax,jmax);

  for (int j=0; j<jmax; ++j)
  {
    for (int i=0; i<jmax; ++i)
    {
      rasterEigen(i,j) = 4;
      rasterEigen(i,j) /= 7;
       //std::cout << rasterTNT[i][j] << ' ';
    }
    //std::cout << std::endl;
  }
  auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  long long microseconds1 = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  std::cout << "Time: " << microseconds1 << std::endl;


  std::cout << "TNT Version: " << std::endl;
  auto start2 = std::chrono::high_resolution_clock::now();
  TNT::Array2D<double> rasterTNT(imax, jmax);

  for (int i=0; i<imax; ++i)
  {
    for (int j=0; j<jmax; ++j)
    {
      rasterTNT[i][j] = 4;
      rasterTNT[i][j] /= 7;
       //std::cout << rasterTNT[i][j] << ' ';
    }
    //std::cout << std::endl;
  }
  auto elapsed2 = std::chrono::high_resolution_clock::now() - start2;
  long long microseconds2 = std::chrono::duration_cast<std::chrono::microseconds>(elapsed2).count();
  std::cout << "Time: " << microseconds2 << std::endl;







//  std::cout << "Array2D Version: " << std::endl;
//  auto start3 = std::chrono::high_resolution_clock::now();
//  LSDArray2D raster2D(imax, jmax);

//  for (int i=0; i<imax; ++i)
//  {
//    for (int j=0; j<jmax; ++j)
//    {
//      raster2D.set(i,j, 0);
//      int val = raster2D.at(i,j) + 7;
//      raster2D.set(i,j, val);
//       //std::cout << rasterTNT[i][j] << ' ';
//    }
//    //std::cout << std::endl;
//  }
//  auto elapsed3 = std::chrono::high_resolution_clock::now() - start3;
//  long long microseconds3 = std::chrono::duration_cast<std::chrono::microseconds>(elapsed3).count();
//  std::cout << "Time: " << microseconds3 << std::endl;



}
