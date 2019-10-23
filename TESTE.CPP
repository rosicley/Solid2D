#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>

using namespace boost::numeric::ublas;

int main()
{

    matrix<double> hammer(7, 3, 0.0);

    vector<double> testando(5000,0.0);

    matrix<double> teste(10000, 10000, 0.0);

    for (int i = 0; i < 10000; i++)
    {
        for (int j = 0; j < 10000; j++)
        {
            teste(i, j)=3.141592;
        }

    }

    double h = 10.0;

    hammer(0, 0) = 1.0 / 3.0;
    hammer(0, 1) = 1.0 / 3.0;
    hammer(0, 2) = 9.0 / 80.0;

    hammer(1, 0) = (9.0 + 2.0 * sqrt(15.0)) / 21.0;
    hammer(1, 1) = (6.0 - sqrt(15.0)) / 21.0;
    hammer(1, 2) = (155.0 - sqrt(15.0)) / 2400.0;

    hammer(2, 0) = (6.0 - sqrt(15.0)) / 21.0;
    hammer(2, 1) = (9.0 + 2.0 * sqrt(15.0)) / 21.0;
    hammer(2, 2) = (155.0 - sqrt(15.0)) / 2400.0;

    hammer(3, 0) = (6.0 - sqrt(15.0)) / 21.0;
    hammer(3, 1) = (6.0 - sqrt(15.0)) / 21.0;
    hammer(3, 2) = (155.0 - sqrt(15.0)) / 2400.0;

    hammer(4, 0) = (6.0 + sqrt(15.0)) / 21.0;
    hammer(4, 1) = (6.0 + sqrt(15.0)) / 21.0;
    hammer(4, 2) = (155.0 + sqrt(15.0)) / 2400.0;

    hammer(5, 0) = (9.0 - 2.0 * sqrt(15.0)) / 21.0;
    hammer(5, 1) = (6.0 + sqrt(15.0)) / 21.0;
    hammer(5, 2) = (155.0 + sqrt(15.0)) / 2400.0;

    hammer(6, 0) = (6.0 + sqrt(15.0)) / 21.0;
    hammer(6, 1) = (9.0 - 2.0 * sqrt(15.0)) / 21.0;
    hammer(6, 2) = (155.0 + sqrt(15.0)) / 2400.0;

    std::cout << sizeof(testando) << " " << sizeof(teste);
    // for(int i = 0; i<7; i++)
    // {
    //     std::cout << hammer(i,0)<<" "<<hammer(i,1)<<" "<<2*hammer(i,2)<<std::endl;
    // }
}