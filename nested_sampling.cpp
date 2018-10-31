/* C++ code for thermodynamic formulation of the nested sampling algorithm, applied to
  exploration of continuous potential energy functions */


#include <iostream>
#include <cmath>
#include <array>
using namespace std;


template <size_t S>
using Dbl_arr = std::array<double, S>;

template <class T, unsigned J, unsigned I>
using Matrix = std::array<std::array<T, J>, I>;

constexpr double pi() { return std::atan(1.)*4.; }

/* Main class for performing a nested sampling (NS) simulation */
class Nested_Sampling {

    protected:

    int n_iter;

    public:

    Nested_Sampling(int n);
    ~Nested_Sampling();
};

// Constructor
Nested_Sampling::Nested_Sampling(int n) {
    n_iter = n;
}

// Destructor
Nested_Sampling::~Nested_Sampling() {
    std::cout << "Called destructor\n";
}

/* Damavandi's function. Takes 2D list of args. Domain ((0.,14.),(0.,14.)).
   Hard-to-locate global min (f=0.) at (2.,2.). Local min (f=2.) at (7.,7.).  */
double damavandi (const Dbl_arr<2>& x) {
    try {
        double denom = pow(pi(),2.)*(x[0]-2.)*(x[1]-2.);
        if (denom < 1.0E-20) { throw std::overflow_error("Divide by zero"); }
        return (1. - pow ( abs( (sin(pi()*(x[0]-2.))*sin(pi()*(x[1]-2.)) \
                ) / denom ), 5.) )  \
                *(2. + pow (x[0]-7.,2.) + 2.*pow(x[1]-7.,2.));
    } catch (std::overflow_error err) {
        return 0.; }
    }


int main () {

    Dbl_arr<2> gm {{2.,2.}};
    Dbl_arr<2> local_min {{7.,7.}};

    Matrix<double,2,2> domain {{ {0.,14.}, {0.,14.} }};
    // Matrix<double,1,2> global_min {{2.,2.}};


    std::cout << "Global min of damavandi:" << damavandi(gm) << "\n";
    std::cout << "local min of damavandi:" << damavandi(local_min) << "\n";


    Nested_Sampling ns1(5);

    return 1;
}
