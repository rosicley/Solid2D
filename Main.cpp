static char help[] = "Code to solve static and dynamic nonlinear analysis of 2D elastic solids";

#include "Solid/Solid.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char*)0, help);

    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    Solid *problem = new Solid;

    problem->setAnalysisParameters("EPD", "T6", 12); //("EPD or EPT", "T3, T6 or T10", 7 or 12)
    //problem->setDynamicAnalysisParameters(0.0001, 0.25, 0.5);

    problem->readAnsysInput("cplusplus.txt");
    problem->readFibersInput("cplusplusFibers.txt");
    
    //problem->exportToParaview(0);
    problem->solveStaticProblem(1, 30, 1.0e-8);
    //solveDynamicProblem(mesmo do estático)


    PetscFinalize();

    return 0;
}