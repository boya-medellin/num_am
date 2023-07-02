#include "lib.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// LINEAR SYSTEMS ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// returns the matrix with the i-th line multiplied by x and added to the j-th line
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Elimination (M mymatrix, int i, int j, double x){

  if ( mymatrix[i].size() == mymatrix[j].size() ){
    for(int it = 0; it < mymatrix[i].size() ; ++it){
      mymatrix[j][it] += mymatrix[i][it] * x;
    }
    return mymatrix;
  }
}


// returns the upper triangular matrix
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Upper_triangular(M mymatrix){
  /* i is collumn iterator , j is lines iterator*/
  for(int i = 0; i < mymatrix.size()  ; ++i){
    mymatrix = Pivot(mymatrix, i, i);
    for(std::size_t j = i + 1; j < mymatrix.size(); ++j){
      double l = - ( mymatrix[j][i] / mymatrix[i][i] );
      mymatrix = Elimination(mymatrix, i, j, l);
    }
  }
  return mymatrix;
}


// performs the Gauss-Jordan method
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Gauss_jordan(M mymatrix){
   for(std::size_t j = 0; j < mymatrix.size()  ; ++j){
    mymatrix = Pivot(mymatrix, j , j);
    for(int i = 0; i < mymatrix.size()  ; ++i){
      if ( i != j){
        double l = - ( mymatrix[i][j] / mymatrix[j][j] );
        mymatrix = Elimination(mymatrix, j, i, l);
      }
    }
  }
  return mymatrix;
}


template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Cramer_replace(M mymatrix, V b, int j){
  if ( j > mymatrix.size() ){
	  std::cout << "cramer_replace::error ' j > mymatrix.size() ' \n";
    return mymatrix;
  }
  for(int i = 0; i < mymatrix.size(); ++i){
    mymatrix[i][j] = b[i];
  }
  return mymatrix;
}


// solves the A*x = b system using the cramer method
template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Cramer(M mymatrix, V b){
  /* only works with upper_tri_det and NOT std_det */
  V x;
  for(std::size_t j = 0 ; j < b.size() ; ++j){
    M Bmatrix = Cramer_replace(mymatrix, b, j);
    x.push_back( Upper_tri_det(Bmatrix) / Upper_tri_det(mymatrix) );
  }
  return x;
}


// solves the A*x = b system
template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Solve(M mymatrix, V b){   
  if ( mymatrix.size() == b.size() ){
    for(int i = 0 ; i < mymatrix.size() ; ++i){   // push b into mymatrix
      mymatrix[i].push_back(b[i]);
    }
    mymatrix = Upper_triangular(mymatrix);        // upper triangular

    V x(mymatrix.size());        
    for(int i = mymatrix.size()-1; i >= 0; --i){
      double sum = 0 ;
      for(std::size_t j = i + 1 ; j<mymatrix.size() ; ++j){
        sum += mymatrix[i][j] * x[j];
      }
      x[i] = ( (mymatrix[i][mymatrix.size()] - sum)/(mymatrix[i][i]) );
    }
    return x;
  } else {
	  std::cout<< "A and b of different length\n";
  }
}


template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Solve_diagonally_dominant(M mymatrix, V b){
  if ( Is_diagonally_dominant(mymatrix) ){
	V x(mymatrix.size());
	for(int i = 0 ; i < x.size(); ++i){
	  x[i] = 0;
	}
	for(std::size_t k = 0 ; k < 25 ; ++k){
	  for(int i = 0 ; i < mymatrix.size() ; ++i){
		double sum1 = 0;
		double sum2 = 0;
		for(std::size_t j = 0 ; j <= i-1; ++j){
		  sum1 += mymatrix[i][j]*x[j] ;
		}
		for(std::size_t j = i + 1; j < mymatrix.size() ; ++j){
		  sum2 += mymatrix[i][j]*x[j] ;
		}
		x[i] = (1/mymatrix[i][i]) * (b[i] - sum1 - sum2 ) ;
	  }
	}
	return x;
  }else{
	  std::cout<< "Not diagonally dominant matrix, returning x vector\n";
	  return b;
  }
}


template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Jacobi(M mymatrix, V b){
  if ( Is_diagonally_dominant(mymatrix) ){
    V x(mymatrix.size());
    for(int i = 0 ; i < x.size(); ++i){
      x[i] = 0;
    }
    for(std::size_t k = 0 ; k < 25 ; ++k){
      for(int i = 0 ; i < mymatrix.size() ; ++i){
        double sum = 0;
        for(std::size_t j = 0; j < mymatrix.size() ; ++j){
          sum += mymatrix[i][j]*x[j] ;
        }
        x[i] += (1/mymatrix[i][i]) * (b[i] - sum) ;
      }
    }
    return x;
  } else {
	  std::cout<< "Not diagonally dominant matrix, returning x vector\n";
	  return b;
  }
}


template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Gauss_seidel(M mymatrix, V b){
  /* Οι απαιτούμενες πράξεις για τον υπολογισμό του x(k+1) στις επα-
  ναληπτικές μεθόδους που παρουσιάσαμε, είναι n2, όπου n η διάσταση του x.
  Η μέθοδος Gauss–Seidel μπορεί να εφαρμοστεί και να συγκλίνει οπωσδήποτε,
  εκτός από τα γραμμικά συστήματα με κυρίαρχη διαγώνιο, και σε συστήματα στα
  οποία ο πίνακας των συντελεστών είναι συμμετρικός θετικά ορισμένος
  */

  if ( Is_diagonally_dominant(mymatrix) ){
    V x(mymatrix.size());
    for(int i = 0 ; i < x.size(); ++i){
      x[i] = 0;
    }
    for(std::size_t k = 0 ; k < 25 ; ++k){
      for(int i = 0 ; i < mymatrix.size() ; ++i){
        double sum1 = 0;
        double sum2 = 0;

        for(std::size_t j = 0 ; j <= i-1; ++j){
          sum1 += mymatrix[i][j]*x[j] ;
        }
        for(std::size_t j = i ; j < mymatrix.size() ; ++j){
          sum2 += mymatrix[i][j]*x[j] ;
        }
        x[i] += (1/mymatrix[i][i]) * (b[i] - sum1 - sum2 ) ;
      }
    }
    return x;
  } else {
	  std::cout<< "Not diagonally dominant matrix, returning x vector\n";
	  return b;
  }
}


// Finds one eigenvalue for mymatrix using the secant method
template<typename T, typename V, typename M, typename F>
T linear<T,V,M,F>::Eigenvalue(M mymatrix, T x0, T x1){
  T constexpr toler{1e-8};
  T x, fnval;
  T fa{Characteristic(mymatrix, x0)}, fb{Characteristic(mymatrix, x1)};

  x = (x1*fa - x0*fb)/(fa - fb);
  fnval = Characteristic(mymatrix, x);

  while (std::abs(fnval) > toler){
	x0 = x1;
	fa = fb;
	x1 = x;
	fb = Characteristic(mymatrix, x);
	x = (x1*fa - x0*fb)/(fa - fb);
	fnval = Characteristic(mymatrix, x);
  }
  return x;
}


// reads matrix from matrix.txt
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Read(){
  M mymatrix; 
  std::ifstream file("matrix.txt", std::ios::in);
  if (file.good()){
	  std::string str;
    while( getline(file, str) ){
		  std::stringstream ss(str);
      T num;
      V eq;
      while (ss >> num){
        eq.push_back(num);
      }
      mymatrix.push_back(eq);
    }
  }
  return mymatrix;
}


// returns an identity matrix of specified size
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Identity_matrix(int size){
  M mymatrix;

  for(int i = 0; i < size; ++i){
    mymatrix.push_back(V() );
    for(std::size_t j = 0; j < size; ++j){
      if ( i == j){
        mymatrix[i].push_back(1);
      } else {
        mymatrix[i].push_back(0);
      }
    }
  }
  return mymatrix;
}


// returns the matrix reduced by the i-th line and j-th column.
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Reduced(M mymatrix, int i, int j){ 
 
  if (mymatrix.size() == 1){
    return mymatrix;
  }

  mymatrix.erase(mymatrix.begin() + i);
  for (typename M::iterator it = mymatrix.begin(); it != mymatrix.end(); ++it){
    it->erase(it->begin() + j);
  }
  return mymatrix;
}


// pushes matrix2 into matrix1
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Push_matrix(M mymatrix1, M mymatrix2){ 

  // TODO : implement size check
  for(int i = 0 ; i < mymatrix1.size() ; ++i){
    for(std::size_t j = 0 ; j < mymatrix1.size() ; ++j){
      mymatrix1[i].push_back( mymatrix2[i][j]);
    }
  }
  return mymatrix1;
}


// returns the inverse matrix
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Inverse(M mymatrix){  

  M identity = Identity_matrix(mymatrix.size());

  mymatrix = Push_matrix(mymatrix, identity);
  mymatrix = Upper_triangular(mymatrix);

  M leftmatrix;
  M rightmatrix;
  for(int i = 0 ; i < mymatrix.size() ; ++i){
    leftmatrix.push_back(V ());
    rightmatrix.push_back(V ());
    for(std::size_t j = 0 ; j < mymatrix.size(); ++j){
      leftmatrix[i].push_back( mymatrix[i][j]);
    }
    for(std::size_t j = mymatrix.size() ; j < 2*mymatrix.size() ; ++j){
      rightmatrix[i].push_back( mymatrix[i][j]);
    }
  }
  M inverse_matrix;

  for(int i = 0 ; i < mymatrix.size() ; ++i){
    V x = Solve( leftmatrix, Column(rightmatrix, i)) ;
    inverse_matrix.push_back(x);
  }

  inverse_matrix = Transpose(inverse_matrix);
  return inverse_matrix;
}


// returns the transposed matrix
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Transpose(M mymatrix){
  M newmatrix;
  for(std::size_t k = 0 ; k  < mymatrix.size() ; ++k){
    newmatrix.push_back(V ());
  }
  for(int i = 0 ; i < mymatrix.size() ; ++i){
    for(std::size_t j = 0 ; j < mymatrix.size() ; ++j){
      newmatrix[j].push_back(mymatrix[i][j]);
    }
  }
  return newmatrix;
}


// returns the matrix with the j-th line added to the i-th line
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Add_lines(M mymatrix, int i, int j){  
  if ( mymatrix[i].size() == mymatrix[j].size() ){
    for(int it = 0 ; it < mymatrix[i].size() ; ++it){
      mymatrix[i][it] += mymatrix[j][it];
    }
  }
  return mymatrix;
}


// returns the matrix with the i-th line switched with the j-the line
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Switch_lines(M mymatrix, int i, int j){ 
  mymatrix[i].swap(mymatrix[j]);
  return mymatrix;
}


// returns matrix with i-th line multiplied by x
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Multiply_lines(M mymatrix, int i, T x){
  for (typename V::iterator it = mymatrix[i].begin(); it != mymatrix[i].end(); ++it){
    *it *= x;
  }
  return mymatrix;
}


// returns the product of two matrices
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Multiply_matrices(M mymatrix1 , M mymatrix2 ){
  M newmatrix;
  if ( mymatrix1[0].size() == mymatrix2.size() ) {   // mxn , nxp
    for(int i = 0; i < mymatrix1.size(); ++i){
	  newmatrix.push_back(V ());
      for(std::size_t j = 0 ; j < mymatrix2.size(); ++j){
        double sum = 0 ;
        for(std::size_t k = 0 ; k < mymatrix2.size(); ++k){
          sum += mymatrix1[i][k] * mymatrix2[k][j];
        }
        newmatrix[i].push_back(sum);
      }
    }
  } else {
	std::cout << "Matrices of invalid sizes\n";
  }
  return newmatrix;
}


// returns the matrix pivoted in the j-th column and i-th line
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Pivot(M mymatrix, int j , int i){
  if ( i == mymatrix.size() - 1 ){
    return mymatrix;
  }
  if ( mymatrix[i][j] < mymatrix[i+1][j]){
    return Switch_lines(mymatrix, i, i+1);
  } else {
    return mymatrix;
  }
}


// returns the matrix pivoted in the j-th line and the i-th column.
template<typename T, typename V, typename M, typename F>
M linear<T,V,M,F>::Reverse_pivot(M mymatrix, int j, int i){

  if ( i == mymatrix.size() -1 ){
    return mymatrix;
  }
  if ( mymatrix[i][j] < mymatrix[i-1][j]){
    return Switch_lines(mymatrix, i -1 , i);
  } else {
    return mymatrix;
  }
}


template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Multiply_matrix_vec(M matrix, V vec){
  V prod;
  for(int i = 0 ; i < vec.size() ; ++i){
    V temp = matrix[i];
    prod.push_back( Multiply_vec(vec, temp)) ;
  }
  return prod;
}

// returns specified line
template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Line(M mymatrix, int n){
  V eq;
  eq = mymatrix[n];
  return eq;
}


// returns specified column
template<typename T, typename V, typename M, typename F>
V linear<T,V,M,F>::Column(M mymatrix, int n){
  V column_vec;
  for (typename M::iterator it=mymatrix.begin(); it != mymatrix.end(); ++it){
    int count = 0;
    for (typename V::iterator it2= it->begin(); it2 != it->end() ; ++it2){
      if (count == n ){
        column_vec.push_back(*it2);
        break;
      }else{
      ++count;
      }
    }   
  }
  return column_vec;
}


// prints matrix
template<typename T, typename V, typename M, typename F>
void linear<T,V,M,F>::Print_matrix(M mymatrix){
  for (typename M::iterator it=mymatrix.begin(); it != mymatrix.end(); ++it){
    for (typename V::iterator it2= it->begin(); it2 != it->end() ; ++it2){
	  std::cout<< *it2 << ' ' ;
    }
	  std::cout << '\n';
  }
}


// prints vector
template<typename T, typename V, typename M, typename F>
void linear<T,V,M,F>::Print_vec(V vec){
  for (typename V::iterator it = vec.begin(); it != vec.end(); ++it){
	  std::cout<< *it << " ";
  }
  std::cout << '\n';
}


// multiplies two vectors and returns product
template<typename T, typename V, typename M, typename F>
T linear<T,V,M,F>::Multiply_vec(V vec1, V vec2){
    if (vec1.size() == vec2.size() ){
        T result = 0;
        for(int i = 0 ; i < vec1.size(); ++i){
            result += vec1[i]*vec2[i];
        }
        return result;
    } else {
	  std::cerr << "Vectors of invalide sizes\n";
    }
}


// returns the determinant of a M using the standard method
template<typename T, typename V, typename M, typename F>
T linear<T,V,M,F>::Std_det(M mymatrix){
  static int count = 0;

  if (mymatrix.size() == 1){
    return mymatrix[0][0];
  }

  static double sum = 0 ;
  static int j = 0;
 
  if (mymatrix.size() == 2) {
    return (( mymatrix[0][0] * mymatrix[1][1] ) - (mymatrix[0][1] * mymatrix[1][0]) ) ;
  } else {
    for(int i = 0 ; i < mymatrix.size() ; ++i){
      if ( (i + j) % 2 == 0 ){    // TODO : use pow()
        sum +=  mymatrix[i][j] * Std_det(Reduced(mymatrix, i, j)); 
        ++count ;
      } else {
        sum -=  mymatrix[i][j] * Std_det(Reduced(mymatrix, i, j));
        ++count;
      }
    }
  }
  if ( count%2 == 0 ){
    return sum;
  }else{
    return -sum;
  }
}


// returns the determinant of a marix using the upper triangular M method
template<typename T, typename V, typename M, typename F>
T linear<T,V,M,F>::Upper_tri_det(M mymatrix){
  mymatrix = Upper_triangular(mymatrix);
  double sum = 1;
  for(int i = 0; i < mymatrix.size(); ++i){
    sum *= mymatrix[i][i];
  }
  return sum;
}


// Characteristic Equation - eigenvalue problem polynomial
template<typename T, typename V, typename M, typename F>
T linear<T,V,M,F>::Characteristic(M mymatrix, T x){
  M matrix_ = mymatrix;
  for(int i=0; i < mymatrix.size(); ++i){
	matrix_[i][i] -= x;
  }
  return Upper_tri_det(matrix_);
}


// returns true for square matrix, false for non-square
template<typename T, typename V, typename M, typename F>
bool linear<T,V,M,F>::Is_square(M mymatrix){
  if (mymatrix.size() == (mymatrix[0].size() ) ){
    return true;
  } else {
    return false;
  }
}


//template<typename T, typename V, typename M, typename F>
//  bool linear<T,V,M,F>::is_symetric(M mymatrix){
// 
//    bool cond1;
//    for(int i = 0 ; i < mymatrix.size() ; ++i){
//      for(std::size_t j = 0 ; j < mymatrix.size(); ++j){
//        for(std::size_t k = i ; k > 0 ; ++k){
//          if
//        }
//      }
//    }
//    }
//  }


// returns true for diagonally dominant M, false for non-diagonally dominant M
template<typename T, typename V, typename M, typename F>
bool linear<T,V,M,F>::Is_diagonally_dominant(M mymatrix){

  bool condition1 = false;             
  // all diagonal elements are larger or equal to the rest of the elements in the column
  bool condition2 = false;
  // at least one diagonal element is larger that the rest of the elements in the column

  // i : collumn iterator, j : row iterator
  for(int i = 0 ; i < mymatrix.size(); ++i)
    for(std::size_t j = 0 ; j < mymatrix.size(); ++j){
      if ( i != j ){
        if ( mymatrix[i][i] >= mymatrix[j][i] ){
          condition1 = true;
          if (mymatrix[i][i] > mymatrix[j][i] ){
            condition2 = true;
          }
        } else {
          return false;
        }
      }
    }
  return (condition1 && condition2);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// NON LINEAR SYSTEMS ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T, typename V, typename M, typename F>
bool nonlinear<T,V,M,F>::Valid(){
  double fa{Function(x0)}, fb{Function(x1)}; 
  if ( x0 > x1 ){
	std::swap(x0, x1);
  }
  if ((fa*fb >= 0) || ( x0 == x1 ) || (fa == fb)){
	std::cerr << "Invalid function or boundary values.\n";
	return false;
  }
  return true;
}


template<typename T, typename V, typename M, typename F>
T nonlinear<T,V,M,F>::Bisection(){
  std::cout << "Bisection Method: \n";
  double constexpr toler{1e-8};
  double x, fnval;
  int iter{0};

  double fa{Function(x0)}, fb{Function(x1)};

  if ( !Valid() ){ return 0; }

  x = (x0+x1)/2.;
  fnval = Function(x);

  while( std::abs(fnval) > toler ){
	++iter;
	if ( std::signbit(fa) == std::signbit(fnval)){
	  x0 = x;
	  fa = fnval;
	}
	else {
	  x1 = x;
	  fb = fnval;
	}
	x = (x0 + x1)/2.;
	fnval = Function(x);
  }
  if (this->verbose){
	std::cout.precision(precision);
	std::cout << "The zero-posistion is " << std::fixed << x <<
	  "\nThe value of the function is " << std::fixed << fnval << 
	  "\nIterations: " << iter << '\n';
  }
  return x;
}


template<typename T, typename V, typename M, typename F>
T  nonlinear<T,V,M,F>::False_position(){
  std::cout << "False Position Method: \n";
  double constexpr toler{1e-8};
  double x, fnval;
  int iter{0};

  double fa{Function(x0)}, fb{Function(x1)};

  if ( !Valid()){ return 0; }

  x = (x1*fa - x0*fb)/(fa - fb);
  fnval = Function(x);

  while( std::abs(fnval) > toler ){
	++iter;
	if ( std::signbit(fa) == std::signbit(fnval)){
	  x0 = x;
	  fa = fnval;
	}
	else {
	  x1 = x;
	  fb = fnval;
	}
	x = (x1*fa - x0*fb)/(fa - fb);
	fnval = Function(x);
  }

  if (verbose){
	std::cout.precision(precision);
	std::cout << "The zero-posistion is " << std::fixed << x <<
	  "\nThe value of the function is " << std::fixed << fnval << 
	  "\nIterations: " << iter << '\n';
  }
  return x;
}


template<typename T, typename V, typename M, typename F>
T nonlinear<T,V,M,F>::Illinois(){
  std::cout << "Illinois Method: \n";
  double constexpr toler{1e-8};
  double x, fnval;

  double fa{Function(x0)}, fb{Function(x1)};

  int iter{0};
  int counta{0}, countb{0};

  if ( !Valid()){ return 0; }

  x = (x1*fa - x0*fb)/(fa - fb);
  fnval = Function(x);

  while(std::abs(fnval) > toler ){
	++iter;
	if ( std::signbit(fa) == std::signbit(fnval)){
	  x0 = x;
	  fa = fnval;
	  ++counta;
	  countb = 0;
	}
	else {
	  x1 = x;
	  fb = fnval;
	  ++countb;
	  counta = 0;
	}
	if ( counta == 2 ){
	  x = (2*x1*fa - x0*fb)/(2*fa - fb);
	} else if ( countb == 2 ){
	  x = (x1*fa - 2*x0*fb)/(fa - 2*fb);
	} else {
	  x = (x1*fa - x0*fb)/(fa - fb);
	}
	fnval = Function(x);
  }

  if (verbose){
	std::cout.precision(precision);
	std::cout << "The zero-posistion is " << std::fixed << x <<
	  "\nThe value of the function is " << std::fixed << fnval << 
	  "\nIterations: " << iter << '\n';
  }
  return x;
}


template<typename T, typename V, typename M, typename F>
T nonlinear<T,V,M,F>::Secant(){
  std::cout << "Secant Method: \n";
  double constexpr toler{1e-8};
  double x, fnval;
  double fa{Function(x0)}, fb{Function(x1)};

  int iter{0};

  x = (x1*fa - x0*fb)/(fa - fb);
  fnval = Function(x);

  while (std::abs(fnval) > toler){
	++iter;
	x0 = x1;
	fa = fb;
	x1 = x;
	fb = Function(x);
	x = (x1*fa - x0*fb)/(fa - fb);
	fnval = Function(x);
  }

  if (verbose){
	std::cout.precision(precision);
	std::cout << "The zero-posistion is " << std::fixed << x <<
	  "\nThe value of the function is " << std::fixed << fnval << 
	  "\nIterations: " << iter << '\n';
  }
  return x;
}


template<typename T, typename V, typename M, typename F>
T nonlinear<T,V,M,F>::Muller(){
  std::cout << "Muller Method: \n";
  double constexpr toler{1e-8};
  double fnval;
  int iter{0};

  auto w0 = [this](T, T){ return (( function(x2) - function(x0))/(x2 - x0)); };
  auto w1 = [this](T, T, T){ return (( function(x2) - function(x1))/(x2 - x1)); };
  auto a  = [&, this](T, T, T){ return (( w1(x0,x1,x2) - w0(x0, x2) )/(x1-x0)); };
  auto b  = [&, this](T, T, T){ return (w0(x0, x2) + a(x0,x1,x2)*(x2-x0)); };
  auto c  = [&, this](T){ return function(x2); };
  auto d  = [&, this](T, T, T){ double pot1 =( b(x0, x1, x2) + sqrt(std::pow(b(x0,x1,x2),2) - 4.*a(x0,x1,x2)*c(x2)) );
							  double pot2 =( b(x0, x1, x2) - sqrt(std::pow(b(x0,x1,x2),2) - 4.*a(x0,x1,x2)*c(x2)) );
							  return ( (std::abs(pot1) > std::abs(pot2)) ? pot1 : pot2);
							  };

  T x3 = x2 - (2.*c(x2)/d(x0,x1,x2));

  //main loop
  while (std::abs(Function(x3)) >= toler ){
	++iter;
	x0 = x1;
	x1 = x2;
	x2 = x3;
	x3 = x2 - (2.*c(x2)/d(x0,x1,x2));
  }
  fnval = Function(x3);

  if (verbose){
	std::cout.precision(precision);
	std::cout << "The zero-posistion is " << std::fixed << x3 <<
	  "\nThe value of the function is " << std::fixed << fnval << 
	  "\nIterations: " << iter << '\n';
  }
  return x3;
}

template<typename T, typename V, typename M, typename F>
T nonlinear<T,V,M,F>::Fixed_point(){
  std::cout << "Fixed Point Method: \n";
  double constexpr toler{1e-8};
  int iter{0};
  T fnval, x{0.5};

  // auto f = [&](T){ return std::pow(x,2) - 6*x + 5; };
  // auto g = [&](T){ return (std::pow(x,2) + 5)/6; };

  while ( std::abs(Function(x)) >= toler) {
	++iter;
	x = Function2(x);
  }

  fnval = Function(x);
  if (verbose){
	std::cout.precision(precision);
	std::cout << "The zero-posistion is " << std::fixed << x <<
	  "\nThe value of the function is " << std::fixed << fnval << 
	  "\nIterations: " << iter << '\n';
  }
  return x;
}

template<typename T, typename V, typename M, typename F>
T nonlinear<T,V,M,F>::Newton_raphson(){
  std::cout << "Newton-Raphson Method: \n";
  double constexpr toler{1e-8};
  int iter{0};
  T fnval, x{1};

  // auto f = [&](T){ return 3*x*std::exp(x) - 1; };
  // auto fder = [&](T){ return 3*std::exp(x) + 3*x*std::exp(x); };

  while ( std::abs( Function(x) ) >= toler) {
	++iter;
    x -= ( Function(x) )/( Function2(x) ) ;
  }
  fnval = Function(x);

  if (verbose){
	std::cout.precision(precision);
	std::cout << "The zero-posistion is " << std::fixed << x <<
	  "\nThe value of the function is " << std::fixed << fnval << 
	  "\nIterations: " << iter << '\n';
  }
  return x;
}



///////////////////////////////////////////////////////////////////////////////////////
/////////////////// F U N C T I O N  A P P R O X I M A T I O N ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

template<typename T, typename V, typename M, typename F>
func<T,V,M,F>::func(F function, T start, T end, int intervals, bool logs, bool verbose):
  function{function}, start{start}, end{end}, intervals{intervals}, logs{logs}, verbose{verbose}
{
  int n = 10 ;
  T step = ((end - start)/n) ;

  // Chebyshev :
  double constexpr pi = 3.14159265359;
  auto x = [&](int i){ return ( (end - start) / 2.) * std::cos((i + 0.5)*pi/(intervals + 1.)) + (end + start)/2.; };
  T x_;
  for(int i = 0; i < intervals; ++i){
	x_ = x(i);
	xdata.push_back(x_);
	data.push_back(function(x_));
  }

  // reverse vectors
  std::reverse(xdata.begin(), xdata.end());
  std::reverse(data.begin(), data.end());

  // example from the book : f(x) = 1/x
  // xdata = {2.0, 2.5, 4.0};
  // data = {0.5, 0.4, 0.25};

  if (logs){
	std::ofstream file("data.dat", std::ios::out);
	if ( file.is_open() ){
	  for(int i = 0 ; i <= intervals ; ++i){
		file << xdata[i] << " " << data[i] << '\n';
	  }
	  file.close();
	}
  }
}


template<typename T, typename V, typename M, typename F>
func<T,V,M,F>::func(V xdata, V data, bool logs, bool verbose):
  xdata{xdata}, data{data}, logs{logs}, verbose{verbose}
{}


template<typename T, typename V, typename M, typename F>
M func<T,V,M,F>::Vandermonde(V x){
  M v;
  for(int i = 0; i < x.size(); ++i){
    v.push_back( V() );
    for(std::size_t j = 0 ; j < x.size() ; ++j){
      v[i].push_back( std::pow(x[i], j)) ;
    }
  }
  return v;
}


template<typename T, typename V, typename M, typename F>
V func<T,V,M,F>::Polynomial(){
  M v = Vandermonde(xdata); 	// Vandermonde matrix
  V coeff, approx; 		 		// Coefficient vector, Approximantion Vector
								// v * coeff = data
								// approx = v * coeff

  coeff = this->Solve(v, data) ;
  approx = this->Multiply_matrix_vec(v, coeff);

  if (logs){
	std::string filename = "poly_method" ;
	Print_to_file(approx, filename) ;
	Abs_error(approx, filename) ;
  }
  if (verbose){
	std::cout << "Polynomial Method : \n";
	std::cout << "Mean Error: " << Mean_error(approx) << '\n';
  }
  return coeff;
}


template<typename T, typename V, typename M, typename F>
T func<T,V,M,F>::Poly_eval(V coeff, T x){
  V x_;
  for(int i = 0; i < coeff.size(); ++i){
	x_.push_back( std::pow(x, i) );
  }
  return this->Multiply_vec(coeff, x_);
}


template<typename T, typename V, typename M, typename F>
V func<T,V,M,F>::Newton_polynomial(){
  M v;
  V coeff, approx;

  for(int i = 0 ; i < xdata.size(); ++i){
	v.push_back( V() ) ;
	int k = 0;
	for(std::size_t j = 0 ; j < xdata.size() ; ++j){
	  if (j == 0) {
		v[i].push_back(1);
	  } else {
		if (j <= i){
		  v[i].push_back( v[i].back() * (xdata[i]-xdata[k]) );
		  ++k;
		} else {
		  v[i].push_back(0);
		}
	  }
	}
  }

  coeff = this->Solve(v, data) ;
  approx = this->Multiply_matrix_vec(v, coeff);

  if (logs){
	std::string filename = "newton_method";
	Print_to_file(approx, filename);
	Abs_error(approx, filename);
  }
  if (verbose){
	std::cout << "Newton Method : \n";
	std::cout << "Mean Error: " << Mean_error(approx) << '\n';
  }
  return coeff;
}


template<typename T, typename V, typename M, typename F>
T func<T,V,M,F>::Newton_eval(V coeff, T x){
  V x_;
  x_.push_back(1);
  for(int i = 1; i < coeff.size(); ++i){
	x_.push_back( x_.back() * (x - xdata[i]) );
  }
  return this->Multiply_vec(coeff, x_);
}


template<typename T, typename V, typename M, typename F>
T func<T,V,M,F>::Lagrance_sum(int i, int start, int stop){
  int n = xdata.size();
  T sum = 0;
  if ( stop > n){
	return 1.;
  }
  for(std::size_t k = start; k < stop; ++k){
	if ( k != i ){
	  sum += Lagrance_sum(i, start+1, stop+1) * xdata[k];
	}
  }
  return sum; 
};

template<typename T, typename V, typename M, typename F>
V func<T,V,M,F>::Lagrance_polynomial(){
  V approx, coeff;
  M l_coeff, v;

  // Vieta formula for calculating the coefficients
  int n = xdata.size();
  T a_n;

  auto sign = [](int i){
	if ( i%2 == 0 ) { return 1; } else { return -1;} 
  };

  auto first_term = [n, this](int i){ 
	T prod = 1;
	for(std::size_t j = 0; j < n; ++j){
	  if ( i != j ){
		prod *= 1/(xdata[i] - xdata[j]);
	  }
	}
	return prod; };

  auto last_term = [n, this](int i){
	T prod = 1;
	for(std::size_t j = 0; j < n; ++j){
	  if ( i != j){
		prod *= xdata[j];
	  }
	}
	return prod;
  };


  for(int i = 0; i < n; ++i){
	l_coeff.push_back( V () );
	a_n = first_term(i);
	l_coeff[i].push_back( a_n );
	for(std::size_t j = 1; j < n; ++j){
	  if ( j == n-1 ){
		l_coeff[i].push_back( sign(j) * a_n * last_term(i));
	  } else {
		l_coeff[i].push_back( sign(j) * a_n * Lagrance_sum(i, 0, n-j+1) );
	  }
	}
  }
  coeff = this->Multiply_matrix_vec(this->Transpose(l_coeff), data);

  // approximation
  v = Vandermonde(xdata);
  approx = this->Multiply_matrix_vec(v, coeff);

  if (logs){
	std::string filename = "lagrance_method";
	Print_to_file(approx, filename);
	Abs_error(approx, filename);
  }
  if (verbose){
	std::cout << "Lagrance Method : \n";
	std::cout << "Mean Error: " << Mean_error(approx) << '\n';
  }
  return coeff;
}


template<typename T, typename V, typename M, typename F>
T func<T,V,M,F>::Lagrance_eval(V coeff, T x){
  V x_;
  for(int i = 0; i < coeff.size(); ++i){
	x_.push_back( std::pow(x, i) );
  }
  std::reverse(x_.begin(), x_.end());
  return this->Multiply_vec(x_, coeff);
}


template<typename T, typename V, typename M, typename F>
V func<T,V,M,F>::Rational_poly(int m, int n){
  M rational_matrix;
  V p_coeff, q_coeff, approx;

  // fill
  for(int i = 0; i < xdata.size(); ++i){
	rational_matrix.push_back(V ());
	for(std::size_t j = 0; j <= m; ++j){
	  rational_matrix[i].push_back( std::pow(xdata[i], j) );
	}
	for(std::size_t k = 1; k <= n; ++k){
	  rational_matrix[i].push_back( -data[i]*std::pow(xdata[i], k) );
	}
  }

  // get coefficients
  V vector = this->Solve(rational_matrix, data);

  // approximation
  p_coeff = std::vector<T>(vector.begin() , vector.begin() + m + 1);
  q_coeff = std::vector<T>(vector.begin() + m + 1, vector.end());
  q_coeff.insert(q_coeff.begin(), 1.);

  auto p = [p_coeff](T x){
	T sum = 0;
	for(int i = 0; i < p_coeff.size(); ++i){
	  sum += p_coeff[i] * std::pow(x, i);
	}
	return sum;
  };
  auto q = [q_coeff](T x){
	T sum = 0;
	for(int i = 0; i < q_coeff.size(); ++i){
	  sum += q_coeff[i] * std::pow(x, i);
	}
	return sum;
  };
  for(int i = 0; i < xdata.size(); ++ i){
	approx.push_back( p(xdata[i]) / q(xdata[i]) );
  }

  if (logs){
	std::string filename = "rational_polynomial_method";
	Print_to_file(approx, filename);
	Abs_error(approx, filename);
  }
  if (verbose){
	std::cout << "Rational Polynomial Method : \n";
	std::cout << "Mean Error: " << Mean_error(approx) << '\n';
  }

  return vector;
}


template<typename T, typename V, typename M, typename F>
V func<T,V,M,F>::Spline(){
  M spline, A;
  V B, coeff;
  T a, b, c, d;
  std::size_t n = xdata.size() - 1;

  // order data
  if ( xdata.size()%2 != 0 ){
	xdata.pop_back();
	data.pop_back();
  }

  // initialize A and B
  for (int i = 0; i < 3*n; ++i){
	A.push_back( V() );
	for (std::size_t j = 0; j < 3*n; ++j){
	  A[i].push_back(0);
	}
  }

  // fill A and b
  for (int i = 0; i < n; ++i){
	std::size_t const j = 3*i;
	A[i][j] 	= std::pow((xdata[i+1] - xdata[i]), 2);
	A[i][j+1] 	= xdata[i+1] - xdata[i];
	A[i][j+2] 	= 1.0;

	B.push_back( (data[i+1] - data[i])/(xdata[i+1] - xdata[i]) );
  }
  for (int i = 0; i < n-1; ++i){
	std::size_t const I = n+i;
	std::size_t const j = 3*i;
	A[I][j] 	= 3.0 * std::pow( (xdata[i+1] - xdata[i]), 2);
	A[I][j+1] 	= 2.0 * (xdata[i+1] - xdata[i]);
	A[I][j+2] 	= 1.0;
	A[I][j+5] 	= -1.0;

	B.push_back( 0 );
  }
  for (int i = 0; i< n-1; ++i){
	std::size_t const I = 2*n-1+i;
	std::size_t const j = 3*i;
	A[I][j] 	= 3.0 * (xdata[i+1] - xdata[i]);
	A[I][j+1] 	= 1.0;
	A[I][j+4] 	= -1.0;

	B.push_back( 0 );
  }
  A[A.size()-2][1] = 1.0;
  A[A.size()-1][A.size() -3] = 3.0 * (xdata[n] - xdata[n-1]);
  A[A.size()-1][A.size() -2] = 1.0;
  B.push_back( 0 );
  B.push_back( 0 );

  std::cout << "B size: " << B.size() << "\n A size: " << A.size() << " " << A[0].size() << '\n';
  coeff = this->Solve(A, B);
  return coeff;
}

template<typename T, typename V, typename M, typename F>
T func<T,V,M,F>::Spline_eval(T x){
}

template<typename T, typename V, typename M, typename F>
V func<T,V,M,F>::Derivative(int degree){
  T x = 1.0;
  M v = this->Transpose(Vandermonde(xdata));
  V b;
  V w = this->Solve(v, b);
}

template<typename T, typename V, typename M, typename F>
void func<T,V,M,F>::Print_to_file(V f_approx, std::string filename){
    
  std::ofstream file( filename + "_data.dat", std::ios::out);
  if ( file.is_open() ){
    for(int i = 0 ; i < f_approx.size() ; ++i){
      file << xdata[i] << " " << f_approx[i] << '\n';
	}
  }
}


template<typename T, typename V, typename M, typename F>
void func<T,V,M,F>::Abs_error(V f_approx, std::string filename){
  V error;

  for(int i = 0; i < xdata.size(); ++i){
    error.push_back( std::abs( f_approx[i] - data[i])) ;
  }
  std::string error_filename = filename + "_error"; 
  Print_to_file( error, error_filename);
}


template<typename T, typename V, typename M, typename F>
T func<T,V,M,F>::Mean_error(V approx){
  double sum=0, error;

  for(int i = 0; i < xdata.size(); ++i){
	error = std::abs( approx[i] - data[i]);
	sum += error;
  }
  return sum/xdata.size();
}


///////////////////////////////////////////////////
/////////////// I N T E G R A L S /////////////////
///////////////////////////////////////////////////


template<typename T, typename V, typename M, typename F>
integral<T,V,M,F>::integral(V xdata, V data, bool logs, bool verbose):
  xdata{xdata}, data{data}, logs{logs}, verbose{verbose}
{}

template<typename T, typename V, typename M, typename F>
integral<T,V,M,F>::integral(F function, T start, T end, int intervals, bool logs, bool verbose):
  function{function}, step{((end - start)/intervals)}, intervals{intervals}, logs{logs}, verbose{verbose}
{
  //T step = ((end - start)/intervals) ;
  T x_= start;
  for (int i=0; i < intervals; ++i){
	xdata.push_back(x_);
	data.push_back(function(x_));
	x_ += step;
  }

  if (logs){
	std::ofstream file("data.dat", std::ios::out);
	if ( file.is_open() ){
	  for(int i = 0 ; i < intervals ; ++i){
		file << xdata[i] << " " << data[i] << '\n';
	  }
	  file.close();
	}
  }
 //  for (int i=0; i < xdata.size(); ++i){
	// std::cout << xdata[i] << ' ' << data[i] << '\n';
 //  }
  std::cout << '\n';
}


template<typename T, typename V, typename M, typename F>
T integral<T,V,M,F>::Table(){
  T sum = 0;
  sum += data[0]/2;
  for (int i = 1; i < data.size() - 1; ++i){
	sum += data[i];
  }
  sum += data[data.size()-1]/2;
  return step*sum;
};


template<typename T, typename V, typename M, typename F>
T integral<T,V,M,F>::Table_error(T M_){
  T min = *std::min_element(xdata.begin(), xdata.end());
  T max = *std::max_element(xdata.begin(), xdata.end());
  std::cout << "step = " << step << '\n';
  return ( (max - min) * M_ * step * step) / 12.;
}


template<typename T, typename V, typename M, typename F>
T integral<T,V,M,F>::Simpson(){
  T sum1 = 0, sum2 = 0,sum3 = 0;
  T min = *std::min_element(xdata.begin(), xdata.end());
  T max = *std::max_element(xdata.begin(), xdata.end());
  T k = (max - min)/(2*step);
  for (int i = 1; i <= k; ++i){
	sum2 += data[2*i -1];
	if ( i < k){
	  sum3 += data[2*i];
	}
  }
  sum1 = (step/3.) * (data[0] + data[2*k] + 4*sum2 + 2*sum3);
  return sum1;
}


template<typename T, typename V, typename M, typename F>
T integral<T,V,M,F>::Simpson_error(T M_){
  T min = *std::min_element(xdata.begin(), xdata.end());
  T max = *std::max_element(xdata.begin(), xdata.end());
  return ( (max - min) * M_ * step * step * step * step) / 180.;
}


template<typename T, typename V, typename M, typename F>
T integral<T,V,M,F>::Simpson_3over8(){
  T sum1 = 0, sum2 = 0, sum3 = 0;
  T k = xdata.size()/3;
  for (int i = 0; i <= k-1; ++i){
	sum1 += data[3*i+1];
	sum2 += data[3*i+2];
	if ( i < k-1 ){
	  sum3 += data[3*i+3];
	}
  }
  return ((3*step)/8) * (data[0] + 3*sum1 + 3*sum2 + 2*sum3 + data[data.size()-1]);
}

/////////////////////////////////////////////////////////
/////////////// D I F E R R E N T I A L /////////////////
/////////////////////////////////////////////////////////

template<typename T, typename V, typename M, typename F>
diff<T,V,M,F>::diff(F function, T y0, T start, T end, int intervals, bool logs, bool verbose):
  function{function}, y0{y0}, step{(end - start)/intervals}, intervals{intervals}, logs{logs}, verbose{verbose}
{
  T x_ = start;
  for (int i = 0; i < intervals; ++i){
	xdata.push_back(x_);
	x_ += step;
  }
  std::cout << "y : " << y0 << '\n';
 //  for (int i = 0; i < intervals; ++i){
	// std::cout << xdata[i] << '\n';
 //  }
}


template<typename T, typename V, typename M, typename F>
T diff<T,V,M,F>::Euler(){
  T y = y0;
  auto kernel = [this](T x, T y){ return y + step * function(x, y); };
  for ( int i = 0; i < xdata.size(); ++i){
	y = kernel(xdata[i], y);
  }
  return y;
}


template<typename T, typename V, typename M, typename F>
T diff<T,V,M,F>::Heun(){
  T y = y0;
  T b = 1./2.;
  auto k1 = 	[&, this](int i){ return step * function(xdata[i], y); };
  auto k2 = 	[&, this](int i){ return step * function(xdata[i] + step, y + k1(i)); };
  auto kernel = [&, this](int i){ return y + b*(k1(i) + k2(i) ); };
  for ( int i = 0; i < xdata.size(); ++i){
	y = kernel(i);
  }
  return y;
}


template<typename T, typename V, typename M, typename F>
T diff<T,V,M,F>::Ralston(){
  T y = y0;
  T b = 1./4.;
  auto k1 = 	[&, this](int i){ return step * function(xdata[i], y); };
  auto k2 = 	[&, this](int i){ return step * function(xdata[i] + (2*step)/2. , y + (2*k1(i))/3.); };
  auto kernel = [&, this](int i){ return y + b*(k1(i) + 3*k2(i) ); };
  for ( int i = 0; i < xdata.size(); ++i){
	y = kernel(i);
  }
  return y;
}


template<typename T, typename V, typename M, typename F>
T diff<T,V,M,F>::RK4(){
  T y = y0;
  T b = 1./6.;
  auto k1 = 	[&, this](int i){ return step * function(xdata[i], y); };
  auto k2 = 	[&, this](int i){ return step * function(xdata[i] + step/2. , y + k1(i)/2.); };
  auto k3 = 	[&, this](int i){ return step * function(xdata[i] + step/2. , y + k2(i)/2.); };
  auto k4 = 	[&, this](int i){ return step * function(xdata[i] + step, y + k3(i)); };
  auto kernel = [&, this](int i){ return y + b*(k1(i) + 2*k2(i) + 2*k3(i) + k4(i) ); };
  for ( int i = 0; i < xdata.size(); ++i){
	y = kernel(i);
  }
  return y;
}


template<typename T, typename V, typename M, typename F>
T diff<T,V,M,F>::RK_3over8(){
  T y = y0;
  T b = 1./8.;
  auto k1 = 	[&, this](int i){ return step * function(xdata[i], y); };
  auto k2 = 	[&, this](int i){ return step * function(xdata[i] + step/3. , y + k1(i)/3.); };
  auto k3 = 	[&, this](int i){ return step * function(xdata[i] + 2*step/3. , y - k1(i)/3. + k2(i)); };
  auto k4 = 	[&, this](int i){ return step * function(xdata[i] + step, y + k1(i) - k2(i) + k3(i)); };
  auto kernel = [&, this](int i){ return y + b*(k1(i) + 3*k2(i) + 3*k3(i) + k4(i) ); };
  for ( int i = 0; i < xdata.size(); ++i){
	y = kernel(i);
  }
  return y;
}


template<typename T, typename V, typename M, typename F>
T diff<T,V,M,F>::Crank_Nicolson(){
  T y = y0;
  auto kernel = [this](T xi, T xi_, T y){ return y + ((xi_ + xi)/2.) * (function(xi, y) + function(xi_, y));};
  for ( int i = 0; i < xdata.size(); ++i){
	y = kernel(xdata[i], xdata[i+1], y);
  }
  return y;
}
