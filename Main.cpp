#include "Solid.h"

int main()
{
    Solid *problem = new Solid;

    problem->setAnalysisParameters("EPT", "T6", 7); //("EPD or EPT", "T3, T6 or T10", 7 or 12)
    //problem->setDynamicAnalysisParameters(0.0001, 0.25, 0.5);

    problem->readAnsysInput("cplusplus.txt");
    
    problem->solveStaticProblem(2, 25, 1.0e-9);
    //solveDynamicProblem(mesmo do estático)



    return 0;
}