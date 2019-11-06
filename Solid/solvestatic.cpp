// int Solid::solveStaticProblem(const int &numberOfSteps, const int &maximumOfIteration, const double &tolerance)
// {
//     Mat A;
//     Vec b, x, All;
//     PetscErrorCode ierr;
//     PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
//     KSP ksp;
//     PC pc;
//     VecScatter ctx;
//     PetscScalar val, value;

//     int rank;

//     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

//     vector<double> externalForces(2 * nodes_.size(), 0.0);
//     externalForces = ExternalForces();
//     int n = (order_ + 1) * (order_ + 2) / 2.0;

//     if (rank == 0)
//     {
//         exportToParaview(0);
//     }

//     PetscMalloc1(dirichletConditions_.size(), &dof);
//     for (size_t i = 0; i < dirichletConditions_.size(); i++)
//     {
//         int indexNode = dirichletConditions_[i]->getNode()->getIndex();
//         int direction = dirichletConditions_[i]->getDirection();
//         dof[i] = (2 * indexNode + direction);
//     }

//     for (int loadStep = 1; loadStep <= numberOfSteps; loadStep++)
//     {
//         if (rank == 0)
//         {
//             std::cout << "------------------------- TIME STEP = "
//                       << loadStep << " -------------------------\n";
//         }

//         for (int iteration = 0; iteration < maximumOfIteration; iteration++) //definir o máximo de interações por passo de carga
//         {
//             //Create PETSc sparse parallel matrix
//             ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
//                                 2 * nodes_.size(), 2 * nodes_.size(),
//                                 100, NULL, 300, NULL, &A);
//             CHKERRQ(ierr);

//             ierr = MatGetOwnershipRange(A, &Istart, &Iend);
//             CHKERRQ(ierr);

//             //Create PETSc vectors
//             ierr = VecCreate(PETSC_COMM_WORLD, &b);
//             CHKERRQ(ierr);
//             ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nodes_.size());
//             CHKERRQ(ierr);
//             ierr = VecSetFromOptions(b);
//             CHKERRQ(ierr);
//             ierr = VecDuplicate(b, &x);
//             CHKERRQ(ierr);
//             ierr = VecDuplicate(b, &All);
//             CHKERRQ(ierr);

//             if (rank == 0)
//             {
//                 externalForces = ((1.0 * loadStep) / (1.0 * numberOfSteps)) * externalForces;
//                 for (size_t i = 0; i < nodes_.size(); i++)
//                 {
//                     if (fabs(externalForces(2 * i)) >= 1.0e-15)
//                     {
//                         int dof = 2 * i;
//                         ierr = VecSetValues(b, 1, &dof, &externalForces(2 * i), ADD_VALUES);
//                     }
//                     if (fabs(externalForces(2 * i + 1)) >= 1.0e-15)
//                     {
//                         int dof = 2 * i + 1;
//                         ierr = VecSetValues(b, 1, &dof, &externalForces(2 * i + 1), ADD_VALUES);
//                     }
//                 }
//             }

//             if(iteration==0)
//             {
//                 for(DirichletCondition* con : dirichletConditions_)
//                 {
//                     Node* node = con->getNode();
//                     int direction = con->getDirection();
//                     double value = (con->getValue())/numberOfSteps;
                    
//                     node->incrementCurrentCoordinate(direction, value);
//                 }
//             }

//             boost::posix_time::ptime t1 =
//                 boost::posix_time::microsec_clock::local_time();
//             for (Element *el : elements_)
//             {
//                 if (elementPartition_[el->getIndex()] == rank)
//                 {
//                     std::pair<vector<double>, matrix<double>> elementMatrices;
//                     elementMatrices = el->elementContributions(planeState_, "STATIC", loadStep, numberOfSteps);

//                     for (size_t i = 0; i < el->getConnection().size(); i++)
//                     {
//                         if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
//                         {
//                             int dof = 2 * el->getConnection()[i]->getIndex();
//                             ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
//                         }
//                         if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
//                         {
//                             int dof = 2 * el->getConnection()[i]->getIndex() + 1;
//                             ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
//                         }

//                         for (size_t j = 0; j < el->getConnection().size(); j++)
//                         {
//                             if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
//                             {
//                                 int dof1 = 2 * el->getConnection()[i]->getIndex();
//                                 int dof2 = 2 * el->getConnection()[j]->getIndex();
//                                 ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
//                             }
//                             if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
//                             {
//                                 int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
//                                 int dof2 = 2 * el->getConnection()[j]->getIndex();
//                                 ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
//                             }
//                             if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
//                             {
//                                 int dof1 = 2 * el->getConnection()[i]->getIndex();
//                                 int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
//                                 ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
//                             }
//                             if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
//                             {
//                                 int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
//                                 int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
//                                 ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
//                             }
//                         }
//                     }
//                 }
//             }
//             boost::posix_time::ptime t2 =
//                 boost::posix_time::microsec_clock::local_time();

//             if (rank == 0)
//             {
//                 boost::posix_time::time_duration diff = t2 - t1;
//                 std::cout << "  Time matriz elementos (s) = " << std::fixed
//                           << diff.total_milliseconds() / 1000. << std::endl;
//             }

//             t1 = boost::posix_time::microsec_clock::local_time();

//             if (fiberInsideSolid_.size() >= 1)
//             {
//                 for (FiberElement *fib : fiberInsideSolid_)
//                 {
//                     if (fiberElementPartition_[fib->getIndex()] == rank)
//                     {
//                         std::pair<vector<double>, matrix<double>> elementMatrices;
//                         elementMatrices = fiberContribution("STATIC", fib);
//                         FiberNode *initial = fib->getConnection()[0];
//                         FiberNode *end = fib->getConnection()[1];
//                         int indexSolidInitical = initial->getIndexSolidElement();
//                         int indexSolidEnd = end->getIndexSolidElement();

//                         vector<int> indexNodes(2 * n);
//                         for (int i = 0; i < n; i++)
//                         {
//                             indexNodes(i) = elements_[indexSolidInitical]->getConnection()[i]->getIndex();
//                             indexNodes(i + n) = elements_[indexSolidEnd]->getConnection()[i]->getIndex();
//                         }

//                         for (size_t i = 0; i < n; i++)
//                         {
//                             if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
//                             {
//                                 int dof = 2 * indexNodes(i);
//                                 ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
//                             }
//                             if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
//                             {
//                                 int dof = 2 * indexNodes(i) + 1;
//                                 ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
//                             }

//                             if (fabs(elementMatrices.first(2 * (i + n))) >= 1.0e-15)
//                             {
//                                 int dof = 2 * indexNodes(i + n);
//                                 ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * (i + n)), ADD_VALUES);
//                             }
//                             if (fabs(elementMatrices.first(2 * (i + n) + 1)) >= 1.0e-15)
//                             {
//                                 int dof = 2 * indexNodes(i + n) + 1;
//                                 ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * (i + n) + 1), ADD_VALUES);
//                             }

//                             for (size_t j = 0; j < n; j++)
//                             {
//                                 if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i);
//                                     int dof2 = 2 * indexNodes(j);
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
//                                 }
//                                 if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i);
//                                     int dof2 = 2 * indexNodes(j) + 1;
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
//                                 }
//                                 if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i) + 1;
//                                     int dof2 = 2 * indexNodes(j);
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
//                                 }
//                                 if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i) + 1;
//                                     int dof2 = 2 * indexNodes(j) + 1;
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
//                                 }
//                                 //segundo elemento
//                                 if (fabs(elementMatrices.second(2 * (i + n), 2 * (j + n))) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i + n);
//                                     int dof2 = 2 * indexNodes(j + n);
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * (j + n)), ADD_VALUES);
//                                 }
//                                 if (fabs(elementMatrices.second(2 * (i + n), 2 * (j + n) + 1)) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i + n);
//                                     int dof2 = 2 * indexNodes(j + n) + 1;
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * (j + n) + 1), ADD_VALUES);
//                                 }
//                                 if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * (j + n))) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i + n) + 1;
//                                     int dof2 = 2 * indexNodes(j + n);
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * (j + n)), ADD_VALUES);
//                                 }
//                                 if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * (j + n) + 1)) >= 1.e-15)
//                                 {
//                                     int dof1 = 2 * indexNodes(i + n) + 1;
//                                     int dof2 = 2 * indexNodes(j + n) + 1;
//                                     ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * (j + n) + 1), ADD_VALUES);
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }

//             t2 = boost::posix_time::microsec_clock::local_time();

//             if (rank == 0)
//             {
//                 boost::posix_time::time_duration diff = t2 - t1;
//                 std::cout << "  Time matriz fibras(s) = " << std::fixed
//                           << diff.total_milliseconds() / 1000. << std::endl;
//             }

//             t1 = boost::posix_time::microsec_clock::local_time();

//             //Assemble matrices and vectors
//             ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//             CHKERRQ(ierr);
//             ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
//             CHKERRQ(ierr);

//             ierr = VecAssemblyBegin(b);
//             CHKERRQ(ierr);
//             ierr = VecAssemblyEnd(b);
//             CHKERRQ(ierr);

//             MatZeroRowsColumns(A, dirichletConditions_.size(), dof, 1.0, x, b);

//             //Create KSP context to solve the linear system
//             ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
//             CHKERRQ(ierr);
//             ierr = KSPSetOperators(ksp, A, A);
//             CHKERRQ(ierr);

//             //Solve using MUMPS
// #if defined(PETSC_HAVE_MUMPS)
//             // ierr = KSPSetType(ksp, KSPCG);
//             // ierr = KSPGetPC(ksp, &pc);
//             // ierr = PCSetType(pc, PCJACOBI);
//             ierr = KSPSetType(ksp, KSPPREONLY);
//             ierr = KSPGetPC(ksp, &pc);
//             ierr = PCSetType(pc, PCLU);
// #endif
//             ierr = KSPSetFromOptions(ksp);
//             CHKERRQ(ierr);
//             ierr = KSPSetUp(ksp);

//             //Solve linear system
//             ierr = KSPSolve(ksp, b, x);
//             CHKERRQ(ierr);
//             ierr = KSPGetTotalIterations(ksp, &iterations);

//             //Gathers the solution vector to the master process
//             ierr = VecScatterCreateToAll(x, &ctx, &All);
//             CHKERRQ(ierr);
//             ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
//             ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
//             ierr = VecScatterDestroy(&ctx);
//             CHKERRQ(ierr);

//             t2 = boost::posix_time::microsec_clock::local_time();

//             if (rank == 0)
//             {
//                 boost::posix_time::time_duration diff = t2 - t1;
//                 std::cout << "  Time sistema (s) = " << std::fixed
//                           << diff.total_milliseconds() / 1000. << std::endl;
//             }

//             //Updates nodal variables
//             double norm = 0.0;
//             Ione = 1;

//             for (size_t i = 0; i < nodes_.size(); i++)
//             {
//                 Idof = 2 * i;
//                 ierr = VecGetValues(All, Ione, &Idof, &val);
//                 CHKERRQ(ierr);
//                 norm += val * val;
//                 nodes_[i]->incrementCurrentCoordinate(0, val);

//                 Idof = 2 * i + 1;
//                 ierr = VecGetValues(All, Ione, &Idof, &val);
//                 CHKERRQ(ierr);
//                 norm += val * val;
//                 nodes_[i]->incrementCurrentCoordinate(1, val);
                
//             }

//             for (FiberElement *fiber : fiberInsideSolid_)
//             {
//                 for (FiberNode *fnode : fiber->getConnection())
//                 {
//                     Element *solid = elements_[fnode->getIndexSolidElement()];
//                     bounded_vector<double, 2> dimensionless = fnode->getDimensionlessCoordinates();
//                     vector<double> phi(n, 0.0);
//                     phi = domainShapeFunction(dimensionless(0), dimensionless(1));
//                     vector<double> correction_X1(n, 0.0);
//                     vector<double> correction_X2(n, 0.0);

//                     int auxiliar = 0;
//                     for (Node *solidNode : solid->getConnection())
//                     {
//                         Idof = 2 * solidNode->getIndex();
//                         ierr = VecGetValues(All, Ione, &Idof, &val);
//                         CHKERRQ(ierr);
//                         correction_X1(auxiliar) = val;

//                         Idof = 2 * solidNode->getIndex() + 1;
//                         ierr = VecGetValues(All, Ione, &Idof, &val);
//                         CHKERRQ(ierr);

//                         correction_X2(auxiliar) = val;
//                         auxiliar = auxiliar + 1;
//                     }

//                     bounded_vector<double, 2> currentCoordinate;
//                     currentCoordinate = fnode->getCurrentCoordinate();
//                     currentCoordinate(0) = currentCoordinate(0) + inner_prod(phi, correction_X1);
//                     currentCoordinate(1) = currentCoordinate(1) + inner_prod(phi, correction_X2);

//                     fnode->setCurrentCoordinate(currentCoordinate);
//                 }
//                 fiber->updateNormalForce();
//             }

//             if (rank == 0)
//             {
//                 //boost::posix_time::time_duration time = t2 - t1;

//                 std::cout << "Iteration = " << iteration
//                           << " (" << iterations << ")"
//                           << "   x Norm = " << std::scientific << sqrt(norm)
//                           << std::endl;
//                 //<< "  Time (s) = " << std::fixed
//                 //<< diff.total_milliseconds()/1000. << std::endl;
//             }

//             ierr = KSPDestroy(&ksp);
//             CHKERRQ(ierr);
//             ierr = VecDestroy(&b);
//             CHKERRQ(ierr);
//             ierr = VecDestroy(&x);
//             CHKERRQ(ierr);
//             ierr = VecDestroy(&All);
//             CHKERRQ(ierr);
//             ierr = MatDestroy(&A);
//             CHKERRQ(ierr);

//             if (sqrt(norm) <= tolerance)
//             {
//                 break;
//             }
//         }

//         if (rank == 0)
//         {
//             for (Node *n : nodes_)
//             {
//                 n->setZeroStressState();
//             }
//             for (int i = 0; i < elements_.size(); i++)
//             {
//                 elements_[i]->StressCalculate(planeState_);
//             }
//             exportToParaview(loadStep);
//         }
//     }
//     PetscFree(dof);
// }
