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

    //ESTADO PLANO, TIPO DE ELEMENTO, NÚMEROS DE PONTOS DE HAMMER
    problem->setAnalysisParameters("EPT", "T6", 12);

    //LENDO ARQUIVO
    problem->readAnsysInput("cplusplus.txt");

    //LENDO ARQUIVO DE FIBRAS, SE EXISTIR
    FILE *fibras;
    fibras = fopen("cplusplusFibers.txt", "r");
    if (fibras != NULL)
    {
        fclose(fibras);
        problem->readFibersInput("cplusplusFibers.txt");
    }

    //RESOLVENDO O PROBLEMA
    problem->solveStaticProblem(1, 30, 1.0e-07);

    //problem->setDynamicAnalysisParameters(0.0001, 0.25, 0.5);
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
