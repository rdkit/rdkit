// file stddefs.H
// D Cosgrove
// Zeneca Pharms
// 27th June 1995
//
// This file has all the standard definitions etc used by the DACLibrary
// and other things

#ifndef DAC_STDDEFS
#define DAC_STDDEFS

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include <cmath>

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

// ****************************************************************************
// cactvs also define YES and NO!
typedef enum { DAC_YES , DAC_NO , DAC_PASS } DAC_YES_NO;

namespace DACLIB {

template <class T> inline T square( T x ) { return x * x; }
template <class T> inline T cube( T x ) { return x * x * x; }

// ****************************************************************************
template <class T> inline T length( const T vector[3] ) {

  T ret_value;
  ret_value = square( vector[0] ) + square( vector[1] ) +
    square( vector[2] );
  return sqrt( ret_value );

}

// ****************************************************************************
template <class T> inline void normalise( T vec[3] ) {

  T len = length( vec );
  vec[0] /= len;
  vec[1] /= len;
  vec[2] /= len;

}

// ****************************************************************************
template <class T> inline void cross_product( const T vec1[3] , const T vec2[3] ,
                                              T cp[3] ) {

  cp[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  cp[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  cp[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

}

// ****************************************************************************
// compute the normalised cross-product pair of vectors
template <class T> inline void norm_cross_product( const T vec1[3] ,
                                                   const T vec2[3] ,
                                                   T cp[3] ) {

  T l;

  // take cross-product
  cross_product( vec1 , vec2 , cp );
  // normalise it
  l = length( cp );

  cp[0] /= l;
  cp[1] /= l;
  cp[2] /= l;

}

// ****************************************************************************
// compute the dot product
template <class T> inline T dot_product( const T vec1[3] , const T vec2[3] ) {

  return( vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2] );

}

// ****************************************************************************
// compute the vector joining the two given points
template <class T> inline void join_vector( const T vec1[3] , const T vec2[3] ,
                                            T vec12[3] ) {

  vec12[0] = vec2[0] - vec1[0];
  vec12[1] = vec2[1] - vec1[1];
  vec12[2] = vec2[2] - vec1[2];

}

// ****************************************************************************
// function to find the cosine of the angle between the first
// vector and the second.  The vector lengths are passed in as well.
template <class T> inline T cos_angle( const T vec1[3] , T len1 ,
                                       const T vec2[3] , T len2 ) {

  T cos_theta;

  cos_theta = dot_product( vec1 , vec2 );
  cos_theta /= ( len1 * len2 );

  return cos_theta;

}

// ****************************************************************************
// function to find the sine of the angle between the first vector and
// the second. Vector lengths passed in as well.
template <class T> inline T sin_angle( const T vec1[3] , T len1 ,
                                       const T vec2[3] , T len2 ) {

  T sin_theta;
  T cp[3] , len_cp;

  cross_product( vec1 , vec2 , cp );
  len_cp = length( cp );
  sin_theta = len_cp / ( len1 * len2 );

  return sin_theta;

}

// ****************************************************************************
// function to find the angle between the two vectors, lengths given,
// relative to the first vector - with ascending values in anti-clockwise
// direction from direction of 1st vector.
template <class T> inline T angle( const T vec1[3] , T len1 ,
                                   const T vec2[3] , T len2 ) {

  T sin_theta = sin_angle( vec1 , len1 , vec2 , len2 );
  T cos_theta = cos_angle( vec1 , len1 , vec2 , len2 );

  // correct for floating point errors - sometimes creeps ever so slightly
  // above 1.0 or below -1.0
  if( sin_theta > 1.0 ) {
    sin_theta = 1.0;
  } else if( sin_theta < -1.0 ) {
    sin_theta = -1.0;
  }

  T theta = fabs( asin( sin_theta ) );

  // find the appropriate quadrant - if both positive, nothing to do
  if( sin_theta > 0.0 && cos_theta < 0.0 )
    // 2nd quadrant - subtract from 180
    theta = M_PI - theta;
  else if( sin_theta < 0.0 && cos_theta < 0.0 )
    // 3rd quadrant, add 180
    theta += M_PI;
  else if( sin_theta < 0.0 && cos_theta > 0.0 )
    // 4th quadrant, subtract from 360
    theta = 2 * M_PI - theta;

  return theta;

}

// ****************************************************************************
// the squared distance between the 2 3D points
template <class T> inline T sq_distance( const T *vec1 , const T *vec2 ) {

  T dist;

  dist = square( vec1[0] - vec2[0] ) + square( vec1[1] - vec2[1] ) +
    square( vec1[2] - vec2[2] );

  return dist;

}

// ****************************************************************************
// the squared distance between 2 vectors of 3D points
template <class T> inline T sq_distance( const T *vec1 , const T *vec2 ,
					 int num_points ) {

  T dist = 0;
  const T *v1 = vec1 , *v2 = vec2;

  int      i;
  for( i = 3 * num_points ; i  ; i-- ) {
    dist += square( *v1 - *v2 );
    v1++;
    v2++;
  }

  return dist;

}

// ****************************************************************************
// the distance between the 2 3D points
template <class T> inline T distance( const T *vec1 , const T *vec2 ) {

  return sqrt( sq_distance( vec1 , vec2 ) );

}

// ****************************************************************************
// the distance between 2 vectors of 3D points
template <class T> inline T distance( const T *vec1 , const T *vec2 ,
                                      int num_points ) {

  return sqrt( sq_distance( vec1 , vec2 , num_points ) );

}

// ****************************************************************************
// rotate about a point
template <class T , class U> inline void rotate( U coords[3] ,
                                                 T rot_matrix[3][3] ,
                                                 T rot_centre[3] ) {

  U cds[3];

  coords[0] -= rot_centre[0];
  coords[1] -= rot_centre[1];
  coords[2] -= rot_centre[2];

  cds[0] = rot_matrix[0][0] * coords[0] +
    rot_matrix[0][1] * coords[1] + rot_matrix[0][2] * coords[2];
  cds[1] = rot_matrix[1][0] * coords[0] +
    rot_matrix[1][1] * coords[1] + rot_matrix[1][2] * coords[2];
  cds[2] = rot_matrix[2][0] * coords[0] +
    rot_matrix[2][1] * coords[1] + rot_matrix[2][2] * coords[2];

  coords[0] = cds[0] + rot_centre[0];
  coords[1] = cds[1] + rot_centre[1];
  coords[2] = cds[2] + rot_centre[2];

}

// ****************************************************************************
// rotate directly
template <typename T , typename U> inline void rotate( U coords[3] ,
                                                       T rot_matrix[3][3] ) {

  U cds[3];

  cds[0] = rot_matrix[0][0] * coords[0] +
    rot_matrix[0][1] * coords[1] + rot_matrix[0][2] * coords[2];
  cds[1] = rot_matrix[1][0] * coords[0] +
    rot_matrix[1][1] * coords[1] + rot_matrix[1][2] * coords[2];
  cds[2] = rot_matrix[2][0] * coords[0] +
    rot_matrix[2][1] * coords[1] + rot_matrix[2][2] * coords[2];

  coords[0] = cds[0]; coords[1] = cds[1]; coords[2] = cds[2];

}

// ****************************************************************************
// translate coords
template <class T , class U> inline void translate( U coords[3] ,
                                                    T x_trans , T y_trans ,
                                                    T z_trans ) {

  coords[0] += x_trans;
  coords[1] += y_trans;
  coords[2] += z_trans;

}

// **************************************************************************
// does the same as the Unix touch command - open a file for writing and
// closes it, this creating it if it doesn't already exist, and updating
// the last access time if it does. Returns true if successful, false otherwise
inline bool touch_file( const std::string filename ) {

  std::ofstream tf( filename.c_str() );
  if( tf && tf.good() )
    return true;
  else
    return false;

}

// **************************************************************************
// write a 3D vector to cout, for debugging use.
template <class T> inline void vec_print( const T *vec ,
                                          bool add_new_line = true ) {
  std::cout << vec[0] << " , " << vec[1] << " , " << vec[2];
  if( add_new_line )
    std::cout << std::endl;

}

// **************************************************************************
// swap two things
template <class T> inline void swap( T* thing1 , T* thing2 ) {
  T *temp_thing = thing1;
  thing1 = thing2;
  thing2 = temp_thing;
}

// **************************************************************************
// make sure the given array is big enough to hold the number sent in
template <class T> inline void make_buffer_big_enough( T *&buff ,
						       int &new_num ,
						       int &curr_size ) {

  if( new_num >= curr_size ) {
    delete [] buff;
    buff = new T[new_num];
    curr_size = new_num;
  }

}

// **************************************************************************
// calculate normal to 2D vector
template <class T> inline void calc_normal_to_2D_vec( T bond[2] ,
                                                      T normal[2] ) {

  T      length , bond_hat[2] , ell[2];

  // check if bond is parallel to an axis
  if( fabs( 0.0 - bond[0] ) < 1.0e-10 ) {
    normal[0] = 1.0;
    normal[1] = 0.0;
    return;
  } else if( fabs( 0.0 - bond[1] ) < 1.0e-10 ) {
    normal[0] = 0.0;
    normal[1] = 1.0;
    return;
  }

  length = bond[0] * bond[0] + bond[1] * bond[1];
  length = sqrt( length );

  // normalise the bond length
  bond_hat[0] = bond[0] / length;
  bond_hat[1] = bond[1] / length;

  // calculate ell, the length of the project of ( the projection of
  // the bond onto the x axis ) back onto the bond
  length = bond[0] * bond_hat[0];
  ell[0] = length * bond_hat[0];
  ell[1] = length * bond_hat[1];

  // calculate the normal to the bond, which is bond_hat scaled
  // by ell minus the projection of the bond on the x axis
  normal[0] = ell[0] - bond[0];
  normal[1] = ell[1];
  length = sqrt( normal[0] * normal[0] + normal[1] * normal[1] );
  normal[0] /= length;
  normal[1] /= length;

}

// **************************************************************************
// make a square matrix of the given size in a memory efficient way. Must be
// deleted in two steps : delete [] t[0]; delete [] t;
template <class T> inline void make_square_matrix( T **&t , int mat_size ) {

  if( mat_size <= 0 )
    t = 0;
  else {
    t = new T *[mat_size];
    t[0] = new T[mat_size * mat_size];
    for( int i = 1 ; i < mat_size ; ++i )
      t[i] = t[i-1] + mat_size;
  }

}

// **************************************************************************
template <class T> inline void destroy_square_matrix( T **&t ) {

  if( t ) {
    delete [] t[0];
    delete [] t;
  }

}

// ***********************************************************************
// make and destroy 2D matrix - more general case of above.
template <class T> inline T** make_2d_matrix( int num_x , int num_y ) {

  if( num_x <= 0 || num_y <= 0 )
    return 0;

  T **ret_mat = new T*[num_x];
  ret_mat[0] = new T[num_x * num_y];
  for( int i = 1 ; i < num_x ; ++i )
    ret_mat[i] = ret_mat[i-1] + num_y;

  return ret_mat;

}

// ***********************************************************************
template <class T> inline void destroy_2d_matrix( T **&in_mat ) {

  delete [] in_mat[0];
  delete [] in_mat;

}

// *******************************************************************
// make and destroy 3D matrix - required only 3 allocations of memory
// and all the values are in 1 continguous array of T elements, so
// that 2 matrices can be compared very rapidly.
template <class T> inline T*** make_3d_matrix( int num_x , int num_y ,
                                               int num_z ) {

  if( num_x <= 0 || num_y <= 0 || num_z <= 0 )
    return 0;

  T ***ret_mat = new T**[num_x];
  ret_mat[0] = new T*[num_x * num_y];
  ret_mat[0][0] = new T[num_x * num_y * num_z];

  for( int j = 1 ; j < num_y ; ++j )
    ret_mat[0][j] = ret_mat[0][j-1] + num_z;

  for( int i = 1 ; i < num_x ; ++i ) {
    ret_mat[i] = ret_mat[i-1] + num_y;
    ret_mat[i][0] = ret_mat[i-1][0] + num_y * num_z;
    for( int j = 1 ; j < num_y ; ++j )
      ret_mat[i][j] = ret_mat[i][j-1] + num_z;
  }

  return ret_mat;

}

//***************************************************************************
template <class T> inline void destroy_3d_matrix( T ***&in_mat ) {

  if( !in_mat )
    return;

  delete [] in_mat[0][0];
  delete [] in_mat[0];
  delete [] in_mat;

  in_mat = 0;

}

//***************************************************************************
// factorials
inline double fac( int n ) {

  double fac = 1.0;
  if( n < 30 ) {
    for( int i = 1 ; i <= n ; ++i )
	fac *= i;
  } else {
    // uses Gosper's approximation, more accurate than Stirling's
    double fltn = double( n );
    double facn = sqrt( ( ( 2 * fltn ) + 1.0 / 3.0 ) * M_PI );
    facn *= pow( fltn , fltn ) * exp( -fltn );
    fac = round( facn );
  }

  return fac;

}

// **********************************************************************
template <class T> inline T radians_to_degrees( T rads ) {
  return rads * 180.0 / M_PI;
}

// **********************************************************************
template <class T> inline T degrees_to_radians( T degs ) {
  return degs * M_PI / 180.0;
}

} // end of namespace

// for easy output, syntax of use:
// copy( rara.begin() , rara.end() , intOut )
static std::ostream_iterator<char> charOut( std::cout , " " );
static std::ostream_iterator<int> intOut( std::cout , " " );
static std::ostream_iterator<unsigned int> uintOut( std::cout , " " );
static std::ostream_iterator<float> floatOut( std::cout , " " );
static std::ostream_iterator<double> doubleOut( std::cout , " " );
static std::ostream_iterator<std::string> stringOut( std::cout , " " );

#endif
