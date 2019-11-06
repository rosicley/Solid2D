static char help[] = "Code to solve static and dynamic nonlinear analysis of 2D elastic solids";

#include "Solid/Solid.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);
    
    boost::posix_time::ptime initial =
        boost::posix_time::microsec_clock::local_time();

    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    Solid *problem = new Solid;

    problem->setAnalysisParameters("EPT", "T6", 12); //("EPD or EPT", "T3, T6 or T10", 7 or 12)
    //problem->setDynamicAnalysisParameters(0.0001, 0.25, 0.5);

    problem->readAnsysInput("cplusplus.txt");
    
    //problem->readFibersInput("cplusplusFibers.txt");

    problem->solveStaticProblem(20, 15, 1.0e-08);
    //solveDynamicProblem(mesmo do estático)

    boost::posix_time::ptime end =
        boost::posix_time::microsec_clock::local_time();

    if (rank == 0)
    {
        boost::posix_time::time_duration diff = end - initial;
        std::cout << "  TOTAL (s) = " << std::fixed
                  << diff.total_milliseconds() / 1000. << std::endl;
    }

    PetscFinalize();

    return 0;
}