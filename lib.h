#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <functional>
#include <sstream>
#include <algorithm>


template<typename T, typename V=std::vector<T>, 
		  typename M=std::vector<std::vector<T>>, typename F=std::function<T(T)>>
class linear{
  public:
	// methods
	M Elimination(M, int, int, double);
	M Upper_triangular(M);
	M Gauss_jordan(M);
	M Cramer_replace(M, V , int);
	V Cramer(M, V);
	V Solve(M, V);
	V Solve_diagonally_dominant(M, V);
	V Jacobi(M, V);
	V Gauss_seidel(M, V);
	T Eigenvalue(M, T, T);
	// helper methods
	M Read();
	M Identity_matrix(int);
	M Reduced(M matrix, int, int);
	M Push_matrix(M, M);
	M Inverse(M);
	M Transpose(M);
	M Add_lines(M, int, int);
	M Switch_lines(M, int, int);
	M Multiply_lines(M, int, T);
	M Multiply_matrices(M, M);
	M Pivot(M matrix, int, int);
	M Reverse_pivot(M, int, int);
	V Multiply_matrix_vec(M, V);
	V Line(M, int);
	V Column(M, int);
	void Print_matrix(M);
	void Print_vec(V);
	T Multiply_vec(V, V);
	T Std_det(M);
	T Upper_tri_det(M);
	T Characteristic(M, T);
	bool Is_square(M);
	bool Is_diagonally_dominant(M);
};


template<typename T, typename V=std::vector<T>, 
		  typename M=std::vector<std::vector<T>>, typename F=std::function<T(T)>>
class nonlinear: public linear<T,V,M,F>{
  private: 
	T x0, x1, x2;
	bool const verbose;
	bool Valid();
	int const precision;
  public:
	nonlinear(T x0, T x1, std::function<T(T)> f, bool const verbose, int const precision):
	  x0{x0}, x1{x1}, Function{f}, verbose{verbose}, precision{precision} {}
	nonlinear(T x0, T x1, T x2, std::function<T(T)> f, std::function<T(T)> g, bool const verbose, int const precision):
	  x0{x0}, x1{x1}, x2{x2}, Function{f}, Function2{g}, verbose{verbose}, precision{precision} {}
	nonlinear(nonlinear const & other):
	  x0{other.x0}, x1{other.x1},  verbose{other.verbose} {}
	~nonlinear() = default;

	T Bisection();
	T False_position();
	T Illinois();
	T Secant();
	T Muller();
	T Fixed_point();
	T Newton_raphson();
	std::function<T(T)> Function;
	std::function<T(T)> Function2;
};


template<typename T, typename V=std::vector<T>, 
		  typename M=std::vector<std::vector<T>>, typename F=std::function<T(T)>>
class func: public linear<T,V,M,F>{
  private:
	F function;
	V data, xdata; 		// data: function values
						// xdata: x values
	T start, end;
	int intervals=50;
	bool logs=false, verbose=true;
  public:
	// constuctor
	// extrapolation of known function from generated data
	func(F function, T start, T end , int intervals, bool logs, bool verbose);
	// constructor
	// extrapolation of function from data input
	func(V xdata, V data, bool logs, bool verbose);
	M Vandermonde(V);
	V Polynomial();
	T Poly_eval(V, T);
	V Newton_polynomial();
	T Newton_eval(V, T);
	T Lagrance_sum(int, int, int);
	V Lagrance_polynomial();
	T Lagrance_eval(V, T);
	V Rational_poly(int, int);
	V Vector_B(T);
	V Spline();
	T Spline_eval(T);
	V Derivative(int);
	// helper methods
	void Print_to_file(V, std::string);
	void Abs_error(V, std::string);
	T Mean_error(V);
};


template<typename T, typename V=std::vector<T>, 
		  typename M=std::vector<std::vector<T>>, typename F=std::function<T(T)>>
class integral: public linear<T,V,M,F>{
  private:
	V data, xdata;
	F function;
	T step;
	int intervals;
	bool logs=false, verbose=false;

  public:
	integral(F function, T start, T end, int intervals, bool logs, bool verbose);
	integral(V data, V xdata, bool logs, bool verbose);
	T Table();
	T Table_error(T M_);
	T Simpson();
	T Simpson_error(T M_);
	T Simpson_3over8();
	T Legendre(T);
	T Gauss_Legendre();
	void Desired_acccuracy(int decimal);
};


template<typename T, typename V=std::vector<T>, 
		  typename M=std::vector<std::vector<T>>, typename F=std::function<T(T,T)>>
class diff: public linear<T,V,M,F>{
  private:
	V xdata;
	F function;
	T y0;
	T step;
	int intervals;
	bool logs=false, verbose=false;

  public:
	diff(F function, T yo, T start, T end, int intervals, bool logs, bool verbose);
	T Euler();
	T Heun();
	T Ralston();
	T RK4();
	T RK_3over8();
	T Crank_Nicolson();


};
