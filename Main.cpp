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

    //ANTES DE LER OS ELEMENTOS
    problem->setDynamicAnalysisParameters(200.0, 0.25, 0.5);

    //ESTADO PLANO, TIPO DE ELEMENTO, NÚMEROS DE PONTOS DE HAMMER
    problem->setAnalysisParameters("EPT", "T6", 7);

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

    //RESOLVER PROBLEMA ESTÁTICO
    problem->solveStaticProblem(10, 20, 1.0e-06);
   
    //RESOLVER PROBLEMA DINÂMICO
    //problem->solveDynamicProblem(500, 10, 1.0e-07);

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
