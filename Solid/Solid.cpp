#include "Solid.h"

Solid::Solid()
{
    deltat_ = 1.0;
    gamma_ = 0.5;
    beta_ = 0.25;
    shapeForces_(0) = 0.0;
    shapeForces_(1) = 0.0;
}

Solid::~Solid() {}

Node *Solid::getNode(const int &index)
{
    return nodes_[index];
}

Element *Solid::getElement(const int &index)
{
    return elements_[index];
}

Material *Solid::getMaterial(const int &index)
{
    return materials_[index];
}

void Solid::addMaterial(const int &index, const double &young, const double &poisson, const double &density)
{
    Material *mat = new Material(index, young, poisson, density);
    materials_.push_back(mat);
}

void Solid::addFiberMaterial(const int &index, const double &young, const double &plastStrain, const double &hardeningModulus, const double &density)
{
    FiberMaterial *mat = new FiberMaterial(index, young, plastStrain, hardeningModulus, density);
    fiberMaterials_.push_back(mat);
}

void Solid::addNode(const int &index, const bounded_vector<double, 2> &initialCoordinate)
{
    Node *n = new Node(index, initialCoordinate);
    nodes_.push_back(n);
}

void Solid::addFiberNode(const int &index, const bounded_vector<double, 2> &initialCoordinate)
{
    FiberNode *n = new FiberNode(index, initialCoordinate);
    fiberNodes_.push_back(n);
}

void Solid::addElement(const int &index, const std::vector<int> &nodesIndex, const int &materialIndex, const double &thickness, const std::string &elementType)
{
    std::vector<Node *> nodes;
    //nodes.reserve(nodesIndex.size());
    for (int i = 0; i < nodesIndex.size(); i++)
    {
        nodes.push_back(nodes_[nodesIndex[i]]);
    }
    Element *e = new Element(index, nodes, materials_[materialIndex], thickness, elementType);

    e->setAnalysisParameters(numberOfHammer_, deltat_, beta_, gamma_);

    e->setShapeForce(shapeForces_);

    elements_.push_back(e);
}

void Solid::addFiberElement(const int &index, const std::vector<int> &connection, const int &material, const double &area)
{
    std::vector<FiberNode *> fiberNodes;
    fiberNodes.push_back(fiberNodes_[connection[0]]);
    fiberNodes.push_back(fiberNodes_[connection[1]]);
    FiberElement *el = new FiberElement(index, fiberNodes, fiberMaterials_[material], area);
    fiberElements_.push_back(el);
}

void Solid::addDirichletCondition(const int &index, const int &direction, const double &value)
{
    DirichletCondition *cond = new DirichletCondition(nodes_[index], direction, value);
    dirichletConditions_.push_back(cond);
}

void Solid::addNeumannCondition(const int &index, const int &direction, const double &value)
{
    NeumannCondition *cond = new NeumannCondition(nodes_[index], direction, value);
    neumannConditions_.push_back(cond);
}

void Solid::setAnalysisParameters(const std::string &planeState, const std::string &elementType, const int &numberOfHammer)
{
    elementType_ = elementType;
    numberOfHammer_ = numberOfHammer;
    planeState_ = planeState;

    if (elementType == "T3")
    {
        order_ = 1;
    }
    else if (elementType == "T6")
    {
        order_ = 2;
    }
    else if (elementType == "T10")
    {
        order_ = 3;
    }
}

void Solid::setDynamicAnalysisParameters(const double &deltat, const double &beta, const double &gama)
{
    deltat_ = deltat;
    beta_ = beta;
    gamma_ = gama;
}

std::pair<vector<double>, matrix<double, column_major>> Solid::fiberContribution(const std::string &typeAnalyze, FiberElement *fib)
{
    int n = (order_ + 1) * (order_ + 2) / 2.0;

    std::pair<vector<double>, matrix<double>> local;

    local = fib->fiberLocalContributions(typeAnalyze, deltat_, beta_);
    vector<double> first = local.first;
    matrix<double> second = local.second;

    FiberNode *initial = fib->getConnection()[0];
    FiberNode *end = fib->getConnection()[1];

    // int indexSolidInitical = initial->getIndexSolidElement();
    // int indexSolidEnd = end->getIndexSolidElement();

    bounded_vector<double, 2> dimensionlessInitial = initial->getDimensionlessCoordinates();
    bounded_vector<double, 2> dimensionlessEnd = end->getDimensionlessCoordinates();

    vector<double> phi_inicial(n, 0.0);
    vector<double> phi_end(n, 0.0);

    phi_inicial = domainShapeFunction(dimensionlessInitial(0), dimensionlessInitial(1));
    phi_end = domainShapeFunction(dimensionlessEnd(0), dimensionlessEnd(1));

    matrix<double> mfibra(4 * n, 4, 0.0);

    for (int i = 0; i < n; i++)
    {
        mfibra(2 * i, 0) = phi_inicial(i);
        mfibra(2 * i + 1, 1) = phi_inicial(i);

        mfibra(2 * n + 2 * i, 2) = phi_end(i);
        mfibra(2 * n + 2 * i + 1, 3) = phi_end(i);
    }

    matrix<double> kfc(4 * n, 4 * n, 0.0);
    matrix<double> aux(4 * n, 4, 0.0);
    aux = prod(mfibra, second);
    kfc = prod(aux, trans(mfibra));

    vector<double> fc(4 * n, 0.0);
    fc = prod(mfibra, first);
    fc = -1.0 * fc;

    return std::make_pair(fc, kfc);
}

vector<double> Solid::ExternalForces()
{
    vector<double> externalForces(2 * nodes_.size(), 0.0);

    for (NeumannCondition *cond : neumannConditions_)
    {
        int indexNode = cond->getNode()->getIndex();
        int direction = cond->getDirection();
        double value = cond->getValue();
        externalForces(2 * indexNode + direction) += value;
    }

    return externalForces;
}

int Solid::solveStaticProblem(const int &numberOfSteps, const int &maximumOfIteration, const double &tolerance)
{
    Mat A;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // vector<double> externalForces(2 * nodes_.size(), 0.0);
    // externalForces = ExternalForces();
    int n = (order_ + 1) * (order_ + 2) / 2.0;

    if (rank == 0)
    {
        exportToParaview(0);
    }

    double initialNorm = 0.0;
    for (Node *node : nodes_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }

    PetscMalloc1(dirichletConditions_.size(), &dof);
    for (size_t i = 0; i < dirichletConditions_.size(); i++)
    {
        int indexNode = dirichletConditions_[i]->getNode()->getIndex();
        int direction = dirichletConditions_[i]->getDirection();
        dof[i] = (2 * indexNode + direction);
    }

    for (int loadStep = 1; loadStep <= numberOfSteps; loadStep++)
    {
        boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0)
        {
            std::cout << "------------------------- LOAD STEP = "
                      << loadStep << " -------------------------\n";
        }

        double norm = 100.0;

        for (int iteration = 0; iteration < maximumOfIteration; iteration++) //definir o máximo de interações por passo de carga
        {
            //Create PETSc sparse parallel matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2 * nodes_.size(), 2 * nodes_.size(),
                                100, NULL, 300, NULL, &A);
            CHKERRQ(ierr);

            ierr = MatGetOwnershipRange(A, &Istart, &Iend);
            CHKERRQ(ierr);

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b);
            CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nodes_.size());
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &x);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All);
            CHKERRQ(ierr);

            if (rank == 0)
            {
                for (NeumannCondition *con : neumannConditions_)
                {
                    int ind = con->getNode()->getIndex();
                    int dir = con->getDirection();
                    double val1 = con->getValue() * (1.0 * loadStep / (1.0 * numberOfSteps));
                    int dof = 2 * ind + dir;
                    ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
                }
                // externalForces = ((1.0 * loadStep) / (1.0 * numberOfSteps)) * externalForces;
                // for (size_t i = 0; i < nodes_.size(); i++)
                // {
                //     if (fabs(externalForces(2 * i)) >= 1.0e-15)
                //     {
                //         int dof = 2 * i;
                //         ierr = VecSetValues(b, 1, &dof, &externalForces(2 * i), ADD_VALUES);
                //     }
                //     if (fabs(externalForces(2 * i + 1)) >= 1.0e-15)
                //     {
                //         int dof = 2 * i + 1;
                //         ierr = VecSetValues(b, 1, &dof, &externalForces(2 * i + 1), ADD_VALUES);
                //     }
                // }
            }

            if (iteration == 0)
            {
                for (DirichletCondition *con : dirichletConditions_)
                {
                    Node *node = con->getNode();
                    int dir = con->getDirection();
                    double val1 = (con->getValue()) / (1.0 * numberOfSteps);

                    node->incrementCurrentCoordinate(dir, val1);
                }
            }

            for (Element *el : elements_)
            {
                if (elementPartition_[el->getIndex()] == rank)
                {
                    std::pair<vector<double>, matrix<double>> elementMatrices;
                    elementMatrices = el->elementContributions(planeState_, "STATIC", loadStep, numberOfSteps);

                    for (size_t i = 0; i < el->getConnection().size(); i++)
                    {
                        if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
                        {
                            int dof = 2 * el->getConnection()[i]->getIndex();
                            ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
                        {
                            int dof = 2 * el->getConnection()[i]->getIndex() + 1;
                            ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
                        }

                        for (size_t j = 0; j < el->getConnection().size(); j++)
                        {
                            if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex();
                                int dof2 = 2 * el->getConnection()[j]->getIndex();
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
                            }
                            if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
                                int dof2 = 2 * el->getConnection()[j]->getIndex();
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
                            }
                            if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex();
                                int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
                            }
                            if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
                                int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
                            }
                        }
                    }
                }
            }

            if (fiberElements_.size() >= 1)
            {
                for (FiberElement *fib : fiberElements_)
                {
                    if (fiberElementPartition_[fib->getIndex()] == rank)
                    {
                        FiberNode *initial = fib->getConnection()[0];
                        FiberNode *end = fib->getConnection()[1];

                        int indexSolidInitial = initial->getIndexSolidElement();
                        int indexSolidEnd = end->getIndexSolidElement();

                        if (indexSolidInitial >= 0 && indexSolidEnd >= 0)
                        {

                            std::pair<vector<double>, matrix<double>> elementMatrices;
                            elementMatrices = fiberContribution("STATIC", fib);

                            vector<int> indexNodes(4 * n, 0.0);
                            for (int i = 0; i < n; i++)
                            {
                                indexNodes(2 * i) = (elements_[indexSolidInitial]->getConnection()[i]->getIndex()) * 2;
                                indexNodes(2 * i + 1) = (elements_[indexSolidInitial]->getConnection()[i]->getIndex()) * 2 + 1;
                                indexNodes(2 * (i + n)) = (elements_[indexSolidEnd]->getConnection()[i]->getIndex()) * 2;
                                indexNodes(2 * (i + n) + 1) = (elements_[indexSolidEnd]->getConnection()[i]->getIndex()) * 2 + 1;
                                // indexNodes(2 * i) = (elements_[indexSolidEnd]->getConnection()[i]->getIndex()) * 2;
                                // indexNodes(2 * i + 1) = (elements_[indexSolidEnd]->getConnection()[i]->getIndex()) * 2 + 1;
                                // indexNodes(2 * (i + n)) = (elements_[indexSolidInitial]->getConnection()[i]->getIndex()) * 2;
                                // indexNodes(2 * (i + n) + 1) = (elements_[indexSolidInitial]->getConnection()[i]->getIndex()) * 2 + 1;
                            }

                            for (int i = 0; i < (4 * n); i++)
                            {
                                if (fabs(elementMatrices.first(i)) >= 1.0e-15)
                                {
                                    int dof = indexNodes(i);
                                    ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(i), ADD_VALUES);
                                }
                                for (int j = 0; j < (4 * n); j++)
                                {
                                    if (fabs(elementMatrices.second(i, j)) >= 1.e-15)
                                    {
                                        int dof1 = indexNodes(i);
                                        int dof2 = indexNodes(j);
                                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(i, j), ADD_VALUES);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);

            ierr = VecAssemblyBegin(b);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);
            CHKERRQ(ierr);

            MatZeroRowsColumns(A, dirichletConditions_.size(), dof, 1.0, x, b);

            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
            CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp, A, A);
            CHKERRQ(ierr);

            //Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp, KSPPREONLY);
            ierr = KSPGetPC(ksp, &pc);
            ierr = PCSetType(pc, PCLU);
#endif
            ierr = KSPSetFromOptions(ksp);
            CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);

            //Solve linear system
            ierr = KSPSolve(ksp, b, x);
            CHKERRQ(ierr);
            ierr = KSPGetTotalIterations(ksp, &iterations);

            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(x, &ctx, &All);
            CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterDestroy(&ctx);
            CHKERRQ(ierr);

            //Updates nodal variables
            norm = 0.0;
            Ione = 1;

            for (size_t i = 0; i < nodes_.size(); i++)
            {
                Idof = 2 * i;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                nodes_[i]->incrementCurrentCoordinate(0, val);

                Idof = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                nodes_[i]->incrementCurrentCoordinate(1, val);
            }

            // if (fiberElements_.size() >= 1)
            // {
            //     for (FiberElement *fiber : fiberElements_)
            //     {
            //         for (FiberNode *fnode : fiber->getConnection())
            //         {
            //             Element *solid = elements_[fnode->getIndexSolidElement()];
            //             bounded_vector<double, 2> dimensionless = fnode->getDimensionlessCoordinates();
            //             vector<double> phi(n, 0.0);
            //             phi = domainShapeFunction(dimensionless(0), dimensionless(1));
            //             vector<double> correction_X1(n, 0.0);
            //             vector<double> correction_X2(n, 0.0);

            //             int auxiliar = 0;
            //             for (Node *solidNode : solid->getConnection())
            //             {
            //                 Idof = 2 * solidNode->getIndex();
            //                 ierr = VecGetValues(All, Ione, &Idof, &val);
            //                 CHKERRQ(ierr);

            //                 correction_X1(auxiliar) = val;

            //                 Idof = 2 * solidNode->getIndex() + 1;
            //                 ierr = VecGetValues(All, Ione, &Idof, &val);
            //                 CHKERRQ(ierr);

            //                 correction_X2(auxiliar) = val;

            //                 auxiliar = auxiliar + 1;
            //             }

            //             bounded_vector<double, 2> currentCoordinate;
            //             currentCoordinate = fnode->getCurrentCoordinate();
            //             currentCoordinate(0) = currentCoordinate(0) + inner_prod(phi, correction_X1);
            //             currentCoordinate(1) = currentCoordinate(1) + inner_prod(phi, correction_X2);

            //             fnode->setCurrentCoordinate(currentCoordinate);
            //         }
            //         fiber->updateNormalForce();
            //     }
            // }

            if (fiberElements_.size() >= 1)
            {
                for (FiberNode *fnode : fiberNodes_)
                {

                    int indexSolidElement = fnode->getIndexSolidElement();
                    if (indexSolidElement != -1)
                    {
                        Element *solid = elements_[indexSolidElement];
                        bounded_vector<double, 2> dimensionless = fnode->getDimensionlessCoordinates();
                        vector<double> phi(n, 0.0);
                        phi = domainShapeFunction(dimensionless(0), dimensionless(1));
                        vector<double> correction_X1(n, 0.0);
                        vector<double> correction_X2(n, 0.0);

                        int auxiliar = 0;
                        for (Node *solidNode : solid->getConnection())
                        {
                            Idof = 2 * solidNode->getIndex();
                            ierr = VecGetValues(All, Ione, &Idof, &val);
                            CHKERRQ(ierr);

                            correction_X1(auxiliar) = val;

                            Idof = 2 * solidNode->getIndex() + 1;
                            ierr = VecGetValues(All, Ione, &Idof, &val);
                            CHKERRQ(ierr);

                            correction_X2(auxiliar) = val;

                            auxiliar = auxiliar + 1;
                        }

                        bounded_vector<double, 2> currentCoordinate;
                        currentCoordinate = fnode->getCurrentCoordinate();
                        currentCoordinate(0) = currentCoordinate(0) + inner_prod(phi, correction_X1);
                        currentCoordinate(1) = currentCoordinate(1) + inner_prod(phi, correction_X2);

                        fnode->setCurrentCoordinate(currentCoordinate);
                    }
                }
            }

            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();

            if (rank == 0)
            {
                boost::posix_time::time_duration diff = t2 - t1;
                std::cout << "Iteration = " << iteration
                          << " (" << loadStep << ")"
                          << "   x Norm = " << std::scientific << sqrt(norm / initialNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds() / 1000. << std::endl;
            }

            ierr = KSPDestroy(&ksp);
            CHKERRQ(ierr);
            ierr = VecDestroy(&b);
            CHKERRQ(ierr);
            ierr = VecDestroy(&x);
            CHKERRQ(ierr);
            ierr = VecDestroy(&All);
            CHKERRQ(ierr);
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);

            if (sqrt(norm / initialNorm) <= tolerance)
            {
                break;
            }
        }

        if (rank == 0)
        {
            for (FiberElement *fib : fiberInsideSolid_)
            {
                fib->updateNormalForce();
            }

            for (Node *n : nodes_)
            {
                n->setZeroStressState();
            }

            for (int i = 0; i < elements_.size(); i++)
            {
                elements_[i]->StressCalculate(planeState_);
            }
            exportToParaview(loadStep);
        }
    }
    PetscFree(dof);
    return 0;
}

void Solid::exportToParaview(const int &loadstep)
{
    std::stringstream text;
    text << "output" << loadstep << ".vtu";
    std::ofstream file(text.str());

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << (nodes_.size() + 2 * fiberInsideSolid_.size())
         << "\"  NumberOfCells=\"" << (elements_.size() + fiberInsideSolid_.size())
         << "\">"
         << "\n";
    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        file << n->getCurrentCoordinate()(0) << " " << n->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            file << fiber->getConnection()[0]->getCurrentCoordinate()(0) << " " << fiber->getConnection()[0]->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
            file << fiber->getConnection()[1]->getCurrentCoordinate()(0) << " " << fiber->getConnection()[1]->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
        }
    }
    file << "      </DataArray>"
         << "\n"
         << "    </Points>"
         << "\n";
    //element connectivity
    file << "    <Cells>"
         << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">"
         << "\n";

    int n = (order_ + 1) * (order_ + 2) / 2.0;

    for (Element *e : elements_)
    {
        std::vector<Node *> conec;
        conec.reserve(n);
        conec = e->getConnection();

        for (int i = 0; i < n; i++)
        {
            file << conec[i]->getIndex() << " ";
        }
        file << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        int initial, end;
        initial = nodes_.size() - 1;
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            initial = initial + 1;
            end = initial + 1;

            file << initial << " " << end << "\n";

            initial = end;
        }
    }
    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    int aux = 0;
    for (Element *e : elements_)
    {
        int n = e->getConnection().size();
        aux += n;
        file << aux << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            aux += 2;
            file << aux << "\n";
        }
    }
    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (Element *e : elements_)
    {
        file << 69 << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            file << 3 << "\n";
        }
    }
    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";
    //nodal results
    file << "    <PointData>"
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Displacement\" format=\"ascii\">"
         << "\n";

    for (Node *n : nodes_)
    {
        bounded_vector<double, 2> initial = n->getInitialCoordinate();
        bounded_vector<double, 2> current = n->getCurrentCoordinate();

        file << current(0) - initial(0) << " " << current(1) - initial(1) << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            bounded_vector<double, 2> initial = fiber->getConnection()[0]->getInitialCoordinate();
            bounded_vector<double, 2> current = fiber->getConnection()[0]->getCurrentCoordinate();
            file << current(0) - initial(0) << " " << current(1) - initial(1) << "\n";

            initial = fiber->getConnection()[1]->getInitialCoordinate();
            current = fiber->getConnection()[1]->getCurrentCoordinate();
            file << current(0) - initial(0) << " " << current(1) - initial(1) << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Velocity\" format=\"ascii\">"
         << "\n";

    for (Node *n : nodes_)
    {
        bounded_vector<double, 2> currentVelocity = n->getCurrentVelocity();
        file << currentVelocity(0) << " " << currentVelocity(1) << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            bounded_vector<double, 2> currentVelocity = fiber->getConnection()[0]->getCurrentVelocity();
            file << currentVelocity(0) << " " << currentVelocity(1) << "\n";

            currentVelocity = fiber->getConnection()[1]->getCurrentVelocity();
            file << currentVelocity(0) << " " << currentVelocity(1) << "\n";
        }
    }

    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Acceleration\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        file << n->getCurrentAcceleration()(0) << " " << n->getCurrentAcceleration()(1) << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            bounded_vector<double, 2> current = fiber->getConnection()[0]->getCurrentAcceleration();
            file << current(0) << " " << current(1) << "\n";

            current = fiber->getConnection()[1]->getCurrentAcceleration();
            file << current(0) << " " << current(1) << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"CauchyNormalStress\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        double cont = n->getStressState()(3);
        double aux1 = n->getStressState()(0);
        double aux2 = n->getStressState()(1);
        file << aux1 / cont << " " << aux2 / cont << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            file << 0.0 << " " << 0.0 << "\n";
            file << 0.0 << " " << 0.0 << "\n";
        }
    }

    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"CauchyShearStress\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        double cont = n->getStressState()(3);
        double aux3 = n->getStressState()(2);
        file << aux3 / cont << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            file << 0.0 << "\n";
            file << 0.0 << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"NormalFibers\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        file << 0.0 << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            file << fiber->getConnection()[0]->getNormalForce() << "\n";
            file << fiber->getConnection()[1]->getNormalForce() << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "    </PointData>"
         << "\n";
    //elemental results
    file << "    <CellData>"
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Process\" format=\"ascii\">" << std::endl;

    for (Element *el : elements_)
    {
        file << elementPartition_[el->getIndex()] << "\n";
    }
    if (fiberInsideSolid_.size() >= 1)
    {
        for (FiberElement *fiber : fiberInsideSolid_)
        {
            file << fiberElementPartition_[fiber->getIndex()] << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";
    file << "    </CellData>"
         << "\n";
    //footnote
    file << "  </Piece>"
         << "\n"
         << "  </UnstructuredGrid>"
         << "\n"
         << "</VTKFile>"
         << "\n";
}

void Solid::readAnsysInput(const std::string &read)
{
    std::ifstream file(read);
    std::string line;

    std::getline(file, line);

    int nmaterial, nnode, nelement, nneumann, ndirichlet;

    file >> nmaterial;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nmaterial; i++)
    {
        int index;
        double young, poisson, density;

        file >> index >> young >> poisson >> density;

        addMaterial(index, young, poisson, density);

        std::getline(file, line);
    }

    std::getline(file, line);

    file >> nnode;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nnode; i++)
    {
        int index;
        double x1, x2;
        bounded_vector<double, 2> coordinate;

        file >> index >> x1 >> x2;

        coordinate(0) = x1;
        coordinate(1) = x2;

        addNode(index, coordinate);

        std::getline(file, line);
    }

    std::getline(file, line);

    file >> nelement;

    std::getline(file, line);
    std::getline(file, line);

    int n = (order_ + 1) * (order_ + 2) / 2.0;

    for (int i = 0; i < nelement; i++)
    {
        std::vector<int> nodesconec;
        int index, material, aux;
        double thickness, b1, b2;

        file >> index;

        for (int j = 0; j < n; j++)
        {
            file >> aux;
            nodesconec.push_back(aux);
        }

        file >> material >> thickness >> b1 >> b2; //shapeForce_(0) >> shapeForce_(1);
        shapeForces_(0) = b1;
        shapeForces_(1) = b2;

        addElement(index, nodesconec, material, thickness, elementType_);

        std::getline(file, line);
    }

    std::getline(file, line);

    file >> nneumann;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nneumann; i++)
    {
        int index, direction;
        double value;

        file >> index >> direction >> value;

        addNeumannCondition(index, direction, value);

        std::getline(file, line);
    }

    std::getline(file, line);

    file >> ndirichlet;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < ndirichlet; i++)
    {
        int index, direction;
        double value;

        file >> index >> direction >> value;

        addDirichletCondition(index, direction, value);

        std::getline(file, line);
    }
    domainDecompositionMETIS(elementType_);
}

void Solid::readFibersInput(const std::string &read)
{
    std::ifstream file(read);
    std::string line;

    std::getline(file, line);

    int nmaterial, nnode, nelement;

    file >> nmaterial;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nmaterial; i++)
    {
        int index;
        double young, plaststrain, hardening, density;

        file >> index >> young >> plaststrain >> hardening >> density;

        addFiberMaterial(index, young, plaststrain, hardening, density);

        std::getline(file, line);
    }

    std::getline(file, line);

    file >> nnode;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nnode; i++)
    {
        int index;
        double x1, x2;
        bounded_vector<double, 2> coordinate;

        file >> index >> x1 >> x2;

        coordinate(0) = x1;
        coordinate(1) = x2;

        addFiberNode(index, coordinate);

        std::getline(file, line);
    }

    std::getline(file, line);

    file >> nelement;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nelement; i++)
    {

        int index, material, initial, end;
        double area;

        file >> index >> initial >> end >> material >> area;

        addFiberElement(index, {initial, end}, material, area);

        std::getline(file, line);
    }

    fibersDecompositionMETIS();

    incidenceOfFibers();
}

void Solid::incidenceOfFibers()
{
    boost::posix_time::ptime t1 =
        boost::posix_time::microsec_clock::local_time();

    matrix<double> dataelems(elements_.size(), 3, 0.0);

    int n = (order_ + 1) * (order_ + 2) / 2.0;

    for (Element *el : elements_)
    {
        int cont = el->getIndex();
        bounded_vector<double, 2> coordinate_0 = el->getConnection()[0]->getInitialCoordinate();
        bounded_vector<double, 2> coordinate_1 = el->getConnection()[1]->getInitialCoordinate();
        bounded_vector<double, 2> coordinate_2 = el->getConnection()[2]->getInitialCoordinate();
        bounded_vector<double, 2> center = (coordinate_0 + coordinate_1 + coordinate_2) / 3.0;
        double r0 = norm_2(center - coordinate_0);
        double r1 = norm_2(center - coordinate_1);
        double r2 = norm_2(center - coordinate_2);
        double rmax = r0;
        if (r1 > rmax)
        {
            rmax = r1;
        }
        if (r2 > rmax)
        {
            rmax = r2;
        }
        dataelems(cont, 0) = center(0);
        dataelems(cont, 1) = center(1);
        dataelems(cont, 2) = rmax;
    }
    for (FiberNode *fnode : fiberNodes_)
    {
        for (Element *el : elements_)
        {
            int cont = el->getIndex();
            bounded_vector<double, 2> coord = fnode->getInitialCoordinate();
            bounded_vector<double, 2> normf;
            normf(0) = coord(0) - dataelems(cont, 0);
            normf(1) = coord(1) - dataelems(cont, 1);

            if (norm_2(normf) <= 2.0 * dataelems(cont, 2))
            {
                vector<double> nos_X1(n, 0.0);
                vector<double> nos_X2(n, 0.0);

                for (int j = 0; j < n; j++)
                {
                    nos_X1(j) = el->getConnection()[j]->getInitialCoordinate()(0);
                    nos_X2(j) = el->getConnection()[j]->getInitialCoordinate()(1);
                }

                bounded_vector<double, 2> p3 = el->getConnection()[2]->getInitialCoordinate();
                bounded_vector<double, 2> p1 = el->getConnection()[0]->getInitialCoordinate();
                bounded_vector<double, 2> p2 = el->getConnection()[1]->getInitialCoordinate();

                bounded_matrix<double, 2, 2> m1t, m;
                bounded_vector<double, 2> auxsol;
                bounded_vector<double, 2> ponto0;
                bounded_vector<double, 2> dsol;

                m1t(0, 0) = p1(0) - p3(0);
                m1t(0, 1) = p2(0) - p3(0);
                m1t(1, 0) = p1(1) - p3(1);
                m1t(1, 1) = p2(1) - p3(1);

                m = inverseMatrix(m1t);

                auxsol = prod(m, (coord - p3));

                double erro = 1.0;
                int contador = 0;

                while ((erro >= 0.000001) && (contador <= 20))
                {
                    double ksi = auxsol(0);
                    double eta = auxsol(1);

                    ponto0(0) = inner_prod(domainShapeFunction(ksi, eta), nos_X1);
                    ponto0(1) = inner_prod(domainShapeFunction(ksi, eta), nos_X2);

                    m1t = el->referenceJacobianMatrix(ksi, eta);
                    m = inverseMatrix(m1t);
                    dsol = prod(m, (coord - ponto0));
                    auxsol = auxsol + dsol;

                    erro = norm_2(dsol);
                    contador = contador + 1;
                }

                bounded_vector<double, 3> res;
                res(0) = auxsol(0);
                res(1) = auxsol(1);
                res(2) = 1.0 - auxsol(0) - auxsol(1);

                if (res(0) >= 0.0 && res(0) <= 1.0 && res(1) >= 0.0 && res(1) <= 1.0 && res(2) >= 0.0 && res(2) <= 1.0)
                {
                    bounded_vector<double, 2> dimensionless;
                    dimensionless(0) = res(0);
                    dimensionless(1) = res(1);

                    fnode->setIncidenceInSolid(cont, dimensionless);

                    break;
                }
            }
        }
    }

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
        for (FiberElement *fib : fiberElements_)
        {
            FiberNode *initial = fib->getConnection()[0];
            FiberNode *end = fib->getConnection()[1];

            int indexSolidInitical = initial->getIndexSolidElement();
            int indexSolidEnd = end->getIndexSolidElement();

            if (indexSolidInitical >= 0 && indexSolidEnd >= 0)
            {
                fiberInsideSolid_.push_back(fib);
            }
        }
        boost::posix_time::ptime t2 =
            boost::posix_time::microsec_clock::local_time();
        boost::posix_time::time_duration diff = t2 - t1;
        std::cout << "INCIDENCIA DAS FIBRAS NO SÓLIDO EM " << std::fixed
                  << diff.total_milliseconds() / 1000. << "SEGUNDOS " << std::endl;
        std::cout << "O NÚMERO DE FIBRAS DENTRO DO SÓLIDO É " << fiberInsideSolid_.size() << std::endl;
    }
}

bounded_matrix<double, 2, 2> Solid::inverseMatrix(const bounded_matrix<double, 2, 2> &matrix)
{
    bounded_matrix<double, 2, 2> inverse;
    double detinv = 1.0 / (matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0));

    inverse(0, 0) = detinv * matrix(1, 1);
    inverse(1, 0) = -detinv * matrix(1, 0);
    inverse(0, 1) = -detinv * matrix(0, 1);
    inverse(1, 1) = detinv * matrix(0, 0);

    return inverse;
}

vector<double> Solid::domainShapeFunction(const double &xsi1, const double &xsi2)
{
    int n = (order_ + 1) * (order_ + 2) / 2.0; //number of nodes per element
    vector<double> phi(n, 0.0);

    if (order_ == 1)
    {
        phi(0) = xsi1;
        phi(1) = xsi2;
        phi(2) = 1.0 - xsi1 - xsi2;
    }
    else if (order_ == 2)
    {
        phi(0) = xsi1 * (2.0 * xsi1 - 1.0);
        phi(1) = xsi2 * (2.0 * xsi2 - 1.0);
        phi(2) = (xsi2 + xsi1 - 1.0) * (2.0 * xsi2 + 2.0 * xsi1 - 1.0);
        phi(3) = 4.0 * xsi1 * xsi2;
        phi(4) = -4.0 * xsi2 * (xsi2 + xsi1 - 1.0);
        phi(5) = -4.0 * xsi1 * (xsi2 + xsi1 - 1.0);
    }
    else if (order_ == 3)
    {
        phi(0) = (xsi1 * (3.0 * xsi1 - 2.0) * (3.0 * xsi1 - 1.0)) / 2.0;
        phi(1) = (xsi2 * (3.0 * xsi2 - 2.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(2) = -((xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0) * (3.0 * xsi2 + 3.0 * xsi1 - 1.0)) / 2.0;
        phi(3) = (9.0 * xsi1 * xsi2 * (3.0 * xsi1 - 1.0)) / 2.0;
        phi(4) = (9.0 * xsi1 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(5) = -(9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(6) = (9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi(7) = (9.0 * xsi1 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi(8) = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0) * (xsi2 + xsi1 - 1.0)) / 2.0;
        phi(9) = -27.0 * xsi1 * xsi2 * (xsi2 + xsi1 - 1.0);
    }

    return phi;
}

void Solid::domainDecompositionMETIS(const std::string &elementType)
{
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());

    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = elements_.size();
    idx_t numNd = nodes_.size();
    idx_t ssize = size;
    idx_t one = 1;
    idx_t n;
    if (elementType == "T3")
        n = 3;
    else if (elementType == "T6")
        n = 6;
    else if (elementType == "T10")
        n = 10;
    idx_t elem_start[numEl + 1], elem_connec[n * numEl];
    elementPartition_ = new idx_t[numEl];
    nodePartition_ = new idx_t[numNd];

    for (idx_t i = 0; i < numEl + 1; i++)
    {
        elem_start[i] = n * i;
    }
    for (idx_t jel = 0; jel < numEl; jel++)
    {

        for (idx_t i = 0; i < elements_[jel]->getConnection().size(); i++)
        {
            int nodeIndex = elements_[jel]->getConnection()[i]->getIndex();
            elem_connec[n * jel + i] = nodeIndex;
        }
    }

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
                       NULL, NULL, &one, &ssize, NULL, NULL,
                       &objval, elementPartition_, nodePartition_);

    mirrorData << std::endl
               << "DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    for (int i = 0; i < elements_.size(); i++)
    {
        mirrorData << "process = " << elementPartition_[i]
                   << ", element = " << i << std::endl;
    }

    mirrorData << std::endl
               << "DOMAIN DECOMPOSITION - NODES" << std::endl;
    for (int i = 0; i < nodes_.size(); i++)
    {
        mirrorData << "process = " << nodePartition_[i]
                   << ", node = " << i << std::endl;
    }
}

void Solid::fibersDecompositionMETIS()
{
    std::string mirror2;
    mirror2 = "fiber_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());

    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = fiberElements_.size();
    idx_t numNd = 2 * fiberElements_.size();
    idx_t ssize = size;
    idx_t one = 1;
    idx_t n;

    n = 2;

    idx_t elem_start[numEl + 1], elem_connec[n * numEl];
    fiberElementPartition_ = new idx_t[numEl];
    fiberNodePartition_ = new idx_t[numNd];

    for (idx_t i = 0; i < numEl + 1; i++)
    {
        elem_start[i] = n * i;
    }
    for (idx_t jel = 0; jel < numEl; jel++)
    {
        for (idx_t i = 0; i < n; i++)
        {
            int nodeIndex = fiberElements_[jel]->getConnection()[i]->getIndex();
            elem_connec[n * jel + i] = nodeIndex;
        }
    }

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
                       NULL, NULL, &one, &ssize, NULL, NULL,
                       &objval, fiberElementPartition_, fiberNodePartition_);

    mirrorData << std::endl
               << "FIBER DECOMPOSITION - ELEMENTS" << std::endl;
    for (int i = 0; i < fiberElements_.size(); i++)
    {
        mirrorData << "process = " << fiberElementPartition_[i]
                   << ", element = " << i << std::endl;
    }

    mirrorData << std::endl
               << "FIBER DECOMPOSITION - NODES" << std::endl;
    for (int i = 0; i < 2 * fiberElements_.size(); i++)
    {
        mirrorData << "process = " << fiberNodePartition_[i]
                   << ", node = " << i << std::endl;
    }
}

int Solid::firstAccelerationCalculation()
{
    Mat A;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int n = (order_ + 1) * (order_ + 2) / 2.0;

    //Create PETSc sparse parallel matrix
    PetscMalloc1(dirichletConditions_.size(), &dof);
    for (size_t i = 0; i < dirichletConditions_.size(); i++)
    {
        int indexNode = dirichletConditions_[i]->getNode()->getIndex();
        int direction = dirichletConditions_[i]->getDirection();
        dof[i] = (2 * indexNode + direction);
    }
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        2 * nodes_.size(), 2 * nodes_.size(),
                        100, NULL, 300, NULL, &A);
    CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nodes_.size());
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &All);
    CHKERRQ(ierr);

    if (rank == 0)
    {
        for (NeumannCondition *con : neumannConditions_)
        {
            int ind = con->getNode()->getIndex();
            int dir = con->getDirection();
            double val1 = con->getValue();
            int dof = 2 * ind + dir;
            ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
        }
    }

    for (DirichletCondition *con : dirichletConditions_)
    {
        Node *node = con->getNode();
        int dir = con->getDirection();
        double val1 = con->getValue();

        node->incrementCurrentCoordinate(dir, val1);
    }

    for (Element *el : elements_)
    {
        if (elementPartition_[el->getIndex()] == rank)
        {
            std::pair<vector<double>, matrix<double>> elementMatrices;
            elementMatrices = el->elementContributions(planeState_, "STATIC", 1, 1);
            matrix<double> massLocal(2 * n, 2 * n, 0.0);
            massLocal = el->massMatrix();

            for (size_t i = 0; i < el->getConnection().size(); i++)
            {
                if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
                {
                    int dof = 2 * el->getConnection()[i]->getIndex();
                    ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
                }
                if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
                {
                    int dof = 2 * el->getConnection()[i]->getIndex() + 1;
                    ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
                }

                for (size_t j = 0; j < el->getConnection().size(); j++)
                {
                    if (fabs(massLocal(2 * i, 2 * j)) >= 1.e-15)
                    {
                        int dof1 = 2 * el->getConnection()[i]->getIndex();
                        int dof2 = 2 * el->getConnection()[j]->getIndex();
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j), ADD_VALUES);
                    }
                    if (fabs(massLocal(2 * i + 1, 2 * j)) >= 1.e-15)
                    {
                        int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
                        int dof2 = 2 * el->getConnection()[j]->getIndex();
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j), ADD_VALUES);
                    }
                    if (fabs(massLocal(2 * i, 2 * j + 1)) >= 1.e-15)
                    {
                        int dof1 = 2 * el->getConnection()[i]->getIndex();
                        int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j + 1), ADD_VALUES);
                    }
                    if (fabs(massLocal(2 * i + 1, 2 * j + 1)) >= 1.e-15)
                    {
                        int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
                        int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j + 1), ADD_VALUES);
                    }
                }
            }
        }
    }

    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ;

    ierr = VecAssemblyBegin(b);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);
    CHKERRQ(ierr);

    MatZeroRowsColumns(A, dirichletConditions_.size(), dof, 1.0, x, b);

    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);

    //Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp, KSPPREONLY);
    ierr = KSPGetPC(ksp, &pc);
    ierr = PCSetType(pc, PCLU);
#endif
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    //Solve linear system
    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &iterations);

    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(x, &ctx, &All);
    CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    CHKERRQ(ierr);

    Ione = 1;

    for (size_t i = 0; i < nodes_.size(); i++)
    {
        bounded_vector<double, 2> firstAccel;
        Idof = 2 * i;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(0) = val;

        Idof = 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(1) = val;
        nodes_[i]->setCurrentAcceleration(firstAccel);
        nodes_[i]->setPastAcceleration(firstAccel);
    }

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    PetscFree(dof);
    if (rank == 0)
    {
        std::cout << "ACELERAÇÕES NO PRIMEIRO PASSO DE TEMPO CALCULADAS." << std::endl;
    }
    return 0;
}

int Solid::solveDynamicProblem(const int &numberOfTimes, const int &maximumOfIteration, const double &tolerance)
{
    firstAccelerationCalculation();

    Mat A;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int n = (order_ + 1) * (order_ + 2) / 2.0;

    if (rank == 0)
    {
        exportToParaview(0);
    }

    double initialNorm = 0.0;
    for (Node *node : nodes_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }

    PetscMalloc1(dirichletConditions_.size(), &dof);
    for (size_t i = 0; i < dirichletConditions_.size(); i++)
    {
        int indexNode = dirichletConditions_[i]->getNode()->getIndex();
        int direction = dirichletConditions_[i]->getDirection();
        dof[i] = (2 * indexNode + direction);
    }

    for (int timeStep = 1; timeStep <= numberOfTimes; timeStep++)
    {
        boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0)
        {
            std::cout << "------------------------- TIME STEP = "
                      << timeStep << " -------------------------\n";
        }

        double norm = 100.0;

        for (int iteration = 0; iteration < maximumOfIteration; iteration++) //definir o máximo de interações por passo de carga
        {
            //Create PETSc sparse parallel matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2 * nodes_.size(), 2 * nodes_.size(),
                                100, NULL, 300, NULL, &A);
            CHKERRQ(ierr);

            ierr = MatGetOwnershipRange(A, &Istart, &Iend);
            CHKERRQ(ierr);

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b);
            CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nodes_.size());
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &x);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All);
            CHKERRQ(ierr);

            if (rank == 0)
            {
                for (NeumannCondition *con : neumannConditions_)
                {
                    int ind = con->getNode()->getIndex();
                    int dir = con->getDirection();
                    double val1 = con->getValue(); //AS FORÇAS APLICADAS PODEM VARIAR AO LONGO DO TEMPO
                    int dof = 2 * ind + dir;
                    ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
                }
            }

            for (Element *el : elements_)
            {
                if (elementPartition_[el->getIndex()] == rank)
                {
                    std::pair<vector<double>, matrix<double>> elementMatrices;
                    elementMatrices = el->elementContributions(planeState_, "DYNAMIC", 1, 1); //COM 1 E 1 AS FORÇAS DE DOMINIO PERMANCEM CONSTATEM AO LONGO DO TEMPO

                    for (size_t i = 0; i < el->getConnection().size(); i++)
                    {
                        if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
                        {
                            int dof = 2 * el->getConnection()[i]->getIndex();
                            ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
                        {
                            int dof = 2 * el->getConnection()[i]->getIndex() + 1;
                            ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
                        }

                        for (size_t j = 0; j < el->getConnection().size(); j++)
                        {
                            if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex();
                                int dof2 = 2 * el->getConnection()[j]->getIndex();
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
                            }
                            if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
                                int dof2 = 2 * el->getConnection()[j]->getIndex();
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
                            }
                            if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex();
                                int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
                            }
                            if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
                            {
                                int dof1 = 2 * el->getConnection()[i]->getIndex() + 1;
                                int dof2 = 2 * el->getConnection()[j]->getIndex() + 1;
                                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
                            }
                        }
                    }
                }
            }

            // if (fiberElements_.size() >= 1)
            // {
            //     for (FiberElement *fib : fiberElements_)
            //     {
            //         if (fiberElementPartition_[fib->getIndex()] == rank)
            //         {
            //             FiberNode *initial = fib->getConnection()[0];
            //             FiberNode *end = fib->getConnection()[1];

            //             int indexSolidInitical = initial->getIndexSolidElement();
            //             int indexSolidEnd = end->getIndexSolidElement();

            //             if (indexSolidInitical >= 0 && indexSolidEnd >= 0)
            //             {

            //                 std::pair<vector<double>, matrix<double>> elementMatrices;
            //                 elementMatrices = fiberContribution("STATIC", fib);
            //                 FiberNode *initial = fib->getConnection()[0];
            //                 FiberNode *end = fib->getConnection()[1];
            //                 int indexSolidInitical = initial->getIndexSolidElement();
            //                 int indexSolidEnd = end->getIndexSolidElement();

            //                 vector<int> indexNodes(2 * n);
            //                 for (int i = 0; i < n; i++)
            //                 {
            //                     indexNodes(i) = elements_[indexSolidInitical]->getConnection()[i]->getIndex();
            //                     indexNodes(i + n) = elements_[indexSolidEnd]->getConnection()[i]->getIndex();
            //                 }

            //                 for (size_t i = 0; i < n; i++)
            //                 {
            //                     if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
            //                     {
            //                         int dof = 2 * indexNodes(i);
            //                         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
            //                     }
            //                     if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
            //                     {
            //                         int dof = 2 * indexNodes(i) + 1;
            //                         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
            //                     }

            //                     if (fabs(elementMatrices.first(2 * (i + n))) >= 1.0e-15)
            //                     {
            //                         int dof = 2 * indexNodes(i + n);
            //                         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * (i + n)), ADD_VALUES);
            //                     }
            //                     if (fabs(elementMatrices.first(2 * (i + n) + 1)) >= 1.0e-15)
            //                     {
            //                         int dof = 2 * indexNodes(i + n) + 1;
            //                         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * (i + n) + 1), ADD_VALUES);
            //                     }

            //                     for (size_t j = 0; j < n; j++)
            //                     {
            //                         if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i);
            //                             int dof2 = 2 * indexNodes(j);
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
            //                         }
            //                         if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i);
            //                             int dof2 = 2 * indexNodes(j) + 1;
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
            //                         }
            //                         if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i) + 1;
            //                             int dof2 = 2 * indexNodes(j);
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
            //                         }
            //                         if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i) + 1;
            //                             int dof2 = 2 * indexNodes(j) + 1;
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
            //                         }
            //                         //segundo elemento
            //                         if (fabs(elementMatrices.second(2 * (i + n), 2 * (j + n))) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i + n);
            //                             int dof2 = 2 * indexNodes(j + n);
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * (j + n)), ADD_VALUES);
            //                         }
            //                         if (fabs(elementMatrices.second(2 * (i + n), 2 * (j + n) + 1)) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i + n);
            //                             int dof2 = 2 * indexNodes(j + n) + 1;
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * (j + n) + 1), ADD_VALUES);
            //                         }
            //                         if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * (j + n))) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i + n) + 1;
            //                             int dof2 = 2 * indexNodes(j + n);
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * (j + n)), ADD_VALUES);
            //                         }
            //                         if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * (j + n) + 1)) >= 1.e-15)
            //                         {
            //                             int dof1 = 2 * indexNodes(i + n) + 1;
            //                             int dof2 = 2 * indexNodes(j + n) + 1;
            //                             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * (j + n) + 1), ADD_VALUES);
            //                         }
            //                     }
            //                 }
            //             }
            //         }
            //     }
            // }

            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);

            ierr = VecAssemblyBegin(b);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);
            CHKERRQ(ierr);

            MatZeroRowsColumns(A, dirichletConditions_.size(), dof, 1.0, x, b);

            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
            CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp, A, A);
            CHKERRQ(ierr);

            //Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp, KSPPREONLY);
            ierr = KSPGetPC(ksp, &pc);
            ierr = PCSetType(pc, PCLU);
#endif
            ierr = KSPSetFromOptions(ksp);
            CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);

            //Solve linear system
            ierr = KSPSolve(ksp, b, x);
            CHKERRQ(ierr);
            ierr = KSPGetTotalIterations(ksp, &iterations);

            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(x, &ctx, &All);
            CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterDestroy(&ctx);
            CHKERRQ(ierr);

            // VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            //Updates nodal variables
            norm = 0.0;
            Ione = 1;

            for (Node *node : nodes_)
            {
                int i = node->getIndex();

                Idof = 2 * i;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                nodes_[i]->incrementCurrentCoordinate(0, val);

                Idof = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;

                nodes_[i]->incrementCurrentCoordinate(1, val);

                bounded_vector<double, 2> vel, accel;
                accel = node->getCurrentCoordinate() / (beta_ * deltat_ * deltat_) - node->getPastCoordinate() / (beta_ * deltat_ * deltat_) -
                        node->getPastVelocity() / (beta_ * deltat_) - node->getPastAcceleration() * (0.5 / beta_ - 1.0);
                node->setCurrentAcceleration(accel);

                vel = gamma_ * deltat_ * node->getCurrentAcceleration() + node->getPastVelocity() + deltat_ * (1.0 - gamma_) * node->getPastAcceleration();
                node->setCurrentVelocity(vel);
            }

            // if (fiberInsideSolid_.size() >= 1)
            // {
            //     for (FiberElement *fiber : fiberInsideSolid_)
            //     {
            //         for (FiberNode *fnode : fiber->getConnection())
            //         {
            //             Element *solid = elements_[fnode->getIndexSolidElement()];
            //             bounded_vector<double, 2> dimensionless = fnode->getDimensionlessCoordinates();
            //             vector<double> phi(n, 0.0);
            //             phi = domainShapeFunction(dimensionless(0), dimensionless(1));
            //             vector<double> correction_X1(n, 0.0);
            //             vector<double> correction_X2(n, 0.0);

            //             int auxiliar = 0;
            //             for (Node *solidNode : solid->getConnection())
            //             {
            //                 Idof = 2 * solidNode->getIndex();
            //                 ierr = VecGetValues(All, Ione, &Idof, &val);
            //                 CHKERRQ(ierr);
            //                 correction_X1(auxiliar) = val;

            //                 Idof = 2 * solidNode->getIndex() + 1;
            //                 ierr = VecGetValues(All, Ione, &Idof, &val);
            //                 CHKERRQ(ierr);

            //                 correction_X2(auxiliar) = val;
            //                 auxiliar = auxiliar + 1;
            //             }

            //             bounded_vector<double, 2> currentCoordinate;
            //             currentCoordinate = fnode->getCurrentCoordinate();
            //             currentCoordinate(0) = currentCoordinate(0) + inner_prod(phi, correction_X1);
            //             currentCoordinate(1) = currentCoordinate(1) + inner_prod(phi, correction_X2);

            //             fnode->setCurrentCoordinate(currentCoordinate);
            //         }
            //         fiber->updateNormalForce();
            //     }
            // }

            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();

            if (rank == 0)
            {
                boost::posix_time::time_duration diff = t2 - t1;
                std::cout << "Iteration = " << iteration
                          << " (" << timeStep << ")"
                          << "   x Norm = " << std::scientific << sqrt(norm / initialNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds() / 1000. << std::endl;
            }

            ierr = KSPDestroy(&ksp);
            CHKERRQ(ierr);
            ierr = VecDestroy(&b);
            CHKERRQ(ierr);
            ierr = VecDestroy(&x);
            CHKERRQ(ierr);
            ierr = VecDestroy(&All);
            CHKERRQ(ierr);
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);

            if (sqrt(norm / initialNorm) <= tolerance)
            {
                break;
            }
        }

        if (rank == 0)
        {
            for (Node *n : nodes_)
            {
                n->setZeroStressState();
            }

            for (int i = 0; i < elements_.size(); i++)
            {
                elements_[i]->StressCalculate(planeState_);
            }
            exportToParaview(timeStep);
        }

        for (Node *node : nodes_)
        {
            bounded_vector<double, 2> coordinate = node->getCurrentCoordinate();
            node->setPastCoordinate(coordinate);
            bounded_vector<double, 2> vel = node->getCurrentVelocity();
            node->setPastVelocity(vel);
            bounded_vector<double, 2> accel = node->getCurrentAcceleration();
            node->setPastAcceleration(accel);
        }

        for (Node *node : nodes_)
        {
            bounded_vector<double, 2> vel, accel;
            accel = node->getCurrentCoordinate() / (beta_ * deltat_ * deltat_) - node->getPastCoordinate() / (beta_ * deltat_ * deltat_) -
                    node->getPastVelocity() / (beta_ * deltat_) - node->getPastAcceleration() * (0.5 / beta_ - 1.0);
            node->setCurrentAcceleration(accel);
            vel = gamma_ * deltat_ * node->getCurrentAcceleration() + node->getPastVelocity() + deltat_ * (1.0 - gamma_) * node->getPastAcceleration();
            node->setCurrentVelocity(vel);
        }
    }
    PetscFree(dof);
    return 0;
}