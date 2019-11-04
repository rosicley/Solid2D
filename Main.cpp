#include "Solid/Solid.h"

int main()
{
    Solid *problem = new Solid;

    problem->setAnalysisParameters("EPD", "T6", 12); //("EPD or EPT", "T3, T6 or T10", 7 or 12)
    //problem->setDynamicAnalysisParameters(0.0001, 0.25, 0.5);

    problem->readAnsysInput("cplusplus.txt");
    problem->readFibersInput("cplusplusFibers.txt");
    
    //problem->exportToParaview(0);
    problem->solveStaticProblem(1, 30, 1.0e-7);
    //solveDynamicProblem(mesmo do estático)



    return 0;
}