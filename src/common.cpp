//==============================================================================
// common.cpp
// Utility functions: approximate equality, formatted printing, and simple
// text/CSV writers; small quadratic least-squares helper via LAPACK.
//==============================================================================

#include "common.hpp"

//------------------------------------------------------------------------------
// Return true if complex numbers are equal within absolute tolerance `tol`.
// Uses component-wise check on real/imag parts to avoid NaN issues.
//------------------------------------------------------------------------------
bool almost_equal(complex_t a, complex_t b, double tol)
{
    return std::abs(a.real() - b.real()) < tol && std::abs(a.imag() - b.imag()) < tol;
}

//------------------------------------------------------------------------------
// Return true if |a - b| < tol (absolute tolerance).
//------------------------------------------------------------------------------
bool almost_equal(double a, double b, double tol)
{
    return std::abs(a - b) < tol;
}

//Check if number is integer
bool isInteger(real_t x)
{
    return std::abs(x - std::round(x)) < 1E-14;
}

//------------------------------------------------------------------------------
// Print real vector with fixed precision; break line after every 10 entries.
//------------------------------------------------------------------------------
void print_vec(const vec_real& vec)
{
    std::cout << std::setprecision(10);
    for (size_t i=0; i<vec.size(); ++i)
    {
        if ((i+1)%10 == 0)
        {
            std::cout << vec[i] << "," << std::endl;
        }
        else if (i == vec.size()-1)
        {
            std::cout << vec[i] << std::endl;
        }
        else
        {
            std::cout << vec[i] << ", ";
        }
    }
}

//------------------------------------------------------------------------------
// Print string vector with fixed precision; break line after every 10 entries.
//------------------------------------------------------------------------------
void print_vec(const std::vector<std::string>& vec)
{
    std::cout << std::setprecision(5);
    for (size_t i=0; i<vec.size(); ++i)
    {
        if ((i+1)%10 == 0)
        {
            std::cout << vec[i] << "," << std::endl;
        }
        else if (i == vec.size()-1)
        {
            std::cout << vec[i] << std::endl;
        }
        else
        {
            std::cout << vec[i] << ", ";
        }
    }
}

//------------------------------------------------------------------------------
// Print complex vector as {real, imag}; wrap lines every 3 entries.
//------------------------------------------------------------------------------
void print_vec(const vec_complex& vec)
{
    std::cout << std::setprecision(16);
    for (size_t i=0; i<vec.size(); ++i)
    {
        if ((i+1)%3 == 0)
        {
            std::cout << "{" << vec[i].real() << ", ";
            std::cout << vec[i].imag() << "}," << std::endl;
        }
        else if (i == vec.size()-1)
        {
            std::cout << "{" << vec[i].real() << ", ";
            std::cout << vec[i].imag() << "}" << std::endl;
        }
        else
        {
            std::cout << "{" << vec[i].real() << ", " << vec[i].imag() << "}, ";
        }
    }
}

//------------------------------------------------------------------------------
// Write real matrix to text file; comma-separated columns, one row per line.
//------------------------------------------------------------------------------
void write_mat(std::string filename, mat_real& mat)
{
    std::ofstream outfile(filename, std::ios::out);

    outfile << std::setprecision(16);
    for (size_t i=0; i<mat.size(); ++i)
    {
        for (size_t j=0; j<mat[i].size()-1; ++j)
        {
            outfile << mat[i][j] << ", ";
        }
        outfile << mat[i][mat[i].size()-1] << std::endl;
    }
}

//------------------------------------------------------------------------------
// Write real vector to text file; one value per line.
//------------------------------------------------------------------------------
void write_vec(std::string filename, vec_real& vec)
{
    std::ofstream outfile(filename, std::ios::out);

    outfile << std::setprecision(16);

    for (size_t i=0; i<vec.size()-1; ++i)
    {
        outfile << vec[i] << std::endl;
    }
    outfile << vec[vec.size()-1] << std::endl;
}

//------------------------------------------------------------------------------
// Write complex vector to text file as "(real, imag)" per line.
//------------------------------------------------------------------------------
void write_vec(std::string filename, vec_complex& vec)
{
    std::ofstream outfile(filename, std::ios::out);

    outfile << std::setprecision(16);

    for (size_t i=0; i<vec.size()-1; ++i)
    {
        outfile <<"(" << vec[i].real() << ", " << vec[i].imag() << "), " << std::endl;
    }
    outfile <<"(" << vec[vec.size()-1].real() << ", " << vec[vec.size()-1].imag() << ")" << std::endl;
}

//------------------------------------------------------------------------------
// Build flat 3x3 quadratic Vandermonde-like design matrix from x1,x2,x3:
// [x1^2 x1 1; x2^2 x2 1; x3^2 x3 1] flattened row-major.
//------------------------------------------------------------------------------
vec_real build_design_matrix(real_t x1, real_t x2, real_t x3)
{
    vec_real design_matrix(3*3, 0);
    
    design_matrix[0] = x1*x1;
    design_matrix[1] = x1;
    design_matrix[2] = 1.0;

    design_matrix[3] = x2*x2;
    design_matrix[4] = x2;
    design_matrix[5] = 1.0;

    design_matrix[6] = x3*x3;
    design_matrix[7] = x3;
    design_matrix[8] = 1.0;

    return design_matrix;
}

//------------------------------------------------------------------------------
// Fit y ≈ a + b x + c x^2 using 3 points via LAPACK least squares (dgels).
// Returns coefficients {a, b, c}. Input y_vals is resized by LAPACK.
//------------------------------------------------------------------------------
vec_real fit_quadratic_least_squares(const vec_real& x_vals, const vec_real& y_vals)
{
    const lapack_int m = 3;  // rows
    const lapack_int n = 3;  // cols
    const lapack_int nrhs = 1;
    const lapack_int lda = n;
    const lapack_int ldb = nrhs;

    vec_real A = build_design_matrix(x_vals[0], x_vals[1], x_vals[2]);
    vec_real b = y_vals;        // LAPACK overwrites b with solution
    b.resize(n, 0.0);

    // Solve A * coeffs = b (row-major). On success, b[0..2] = {a,b,c}.
    lapack_int info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, A.data(), lda, b.data(), ldb);
    if (info != 0)
    {
        throw std::runtime_error("LAPACKE_dgels failed with info = " + std::to_string(info));
    }

    return {b[0], b[1], b[2]};
}

//------------------------------------------------------------------------------
// Extracts elements at indices i where (i % n == 0) or (i == last_index).
// Ensures both the first and last elements are always included (closed interval).
// Returns the subset vector. Input 'a' remains unmodified.
//------------------------------------------------------------------------------
vec_real every_nth_element(const vec_real& a, size_t n)
{
    size_t i = 0;
    size_t last_idx = a.size()-1;
    vec_real tmp_vec;
    std::for_each(a.begin(), a.end(), [&i,&tmp_vec,n,last_idx](auto& ele) { 
        if ((i%n==0) || (i==last_idx))
        {
            tmp_vec.push_back(ele);
        }
        ++i;
    });
    
    return tmp_vec;
}
