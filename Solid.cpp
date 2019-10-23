#include "Solid.h"

std::vector<Node *> Truss::getNodes()
{
    return nodes_;
}

vector<int> Truss::BoundaryCondition()
{
    vector<double> boundary(3 * nodes_.size(), 0.0);

    for (int i = 0; i < 3 * nodes_.size(); i++)
    {
        boundary[i] = boundaryConditions_[i];
    }

    return boundary;
}

std::vector<Element *> Truss::getElements()
{
    return elements_;
}

std::vector<Material *> Truss::getMaterial()
{
    return materials_;
}

void Truss::addMaterial(const int &index,
                        const double &young,
                        const double &plastStrain,
                        const double &hardeningModulus,
                        const double &density,
                        const double &expansionCoef)
{
    Material *mat = new Material(index, young, plastStrain, hardeningModulus, density, expansionCoef);
    materials_.push_back(mat);
}

void Truss::addNode(const int &index, const std::vector<double> &initialCoordinate)
{
    Node *n = new Node(index, initialCoordinate);
    nodes_.push_back(n);
}

void Truss::addElement(const int &index,
                       const std::vector<int> &connection,
                       const int &material,
                       const double &area)
{
    std::vector<Node *> nodes;
    nodes.push_back(nodes_[connection[0]]);
    nodes.push_back(nodes_[connection[1]]);
    Element *el = new Element(index, nodes, materials_[material], area);
    elements_.push_back(el);
}

vector<double> Truss::InternalForces()
{
    int cont = elements_.size();
    vector<double> internalForces_(3 * nodes_.size(), 0.0);

    for (int i = 0; i < cont; i++)
    {
        std::vector<double> force = elements_[i]->InternalForce();

        int indexInitialNode = elements_[i]->getConnection()[0]->getIndex();
        int indexEndNode = elements_[i]->getConnection()[1]->getIndex();

        for (int id = 0; id < 3; id++)
        {
            internalForces_[3 * indexInitialNode + id] += force[id];
        }

        for (int id = 0; id < 3; id++)
        {
            internalForces_[3 * indexEndNode + id] += force[3 + id];
        }
    }
    return internalForces_;
}

vector<double> Truss::TemperatureForces(const int &numberOfSteps, const int &currentStep)
{
    int cont = elements_.size();
    vector<double> tempForces_(3 * nodes_.size(), 0.0);

    for (int i = 0; i < cont; i++)
    {
        std::vector<double> force = elements_[i]->TemperatureForce(numberOfSteps, currentStep);

        int indexInitialNode = elements_[i]->getConnection()[0]->getIndex();
        int indexEndNode = elements_[i]->getConnection()[1]->getIndex();

        for (int id = 0; id < 3; id++)
        {
            tempForces_[3 * indexInitialNode + id] += force[id];
        }

        for (int id = 0; id < 3; id++)
        {
            tempForces_[3 * indexEndNode + id] += force[3 + id];
        }
    }
    return tempForces_;
}

matrix<double> Truss::MassMatrix()
{
    matrix<double> Mass(3 * nodes_.size(), 3 * nodes_.size(), 0.0);

    for (int i = 0; i < elements_.size(); i++)
    {
        bounded_matrix<double, 6, 6> localMass = elements_[i]->localMassMatrix();
        int indexInitialNode = elements_[i]->getConnection()[0]->getIndex();
        int indexEndNode = elements_[i]->getConnection()[1]->getIndex();
        for (int ij = 0; ij < 3; ij++)
        {
            for (int ik = 0; ik < 3; ik++)
            {
                Mass(indexInitialNode * 3 + ij, indexInitialNode * 3 + ik) += localMass(ij, ik);

                Mass(indexInitialNode * 3 + ij, indexEndNode * 3 + ik) += localMass(ij, ik + 3);

                Mass(indexEndNode * 3 + ik, indexInitialNode * 3 + ij) += localMass(ij + 3, ik);

                Mass(indexEndNode * 3 + ij, indexEndNode * 3 + ik) += localMass(ij + 3, ik + 3);
            }
        }
    }
    return Mass;
}

vector<double> Truss::InertialForces(const double &beta, const double &gamma, const double &deltat)
{
    vector<double> currentCoordinate(3 * nodes_.size(), 0.0), pastCoordinate(3 * nodes_.size(), 0.0),
        currentVelocity(3 * nodes_.size(), 0.0), pastVelocity(3 * nodes_.size(), 0.0),
        currentAcceleration(3 * nodes_.size(), 0.0), pastAcceleration(3 * nodes_.size(), 0.0),
        inertialForce(3 * nodes_.size(), 0.0), Q(3 * nodes_.size(), 0.0), R(3 * nodes_.size(), 0.0);
    matrix<double> massMatrix(3 * nodes_.size(), 3 * nodes_.size(), 0.0), amortecimento(3 * nodes_.size(), 3 * nodes_.size(), 0.0);

    massMatrix = MassMatrix();

    amortecimento = 0.05 * massMatrix;

    for (Node *node : nodes_)
    {
        int index = node->getIndex();

        pastCoordinate(3 * index) = node->getPastCoordinate()[0];
        pastCoordinate(3 * index + 1) = node->getPastCoordinate()[1];
        pastCoordinate(3 * index + 2) = node->getPastCoordinate()[2];

        pastVelocity(3 * index) = node->getPastVelocity()[0];
        pastVelocity(3 * index + 1) = node->getPastVelocity()[1];
        pastVelocity(3 * index + 2) = node->getPastVelocity()[2];

        pastAcceleration(3 * index) = node->getPastAcceleration()[0];
        pastAcceleration(3 * index + 1) = node->getPastAcceleration()[1];
        pastAcceleration(3 * index + 2) = node->getPastAcceleration()[2];

        currentCoordinate(3 * index) = node->getCurrentCoordinate()[0];
        currentCoordinate(3 * index + 1) = node->getCurrentCoordinate()[1];
        currentCoordinate(3 * index + 2) = node->getCurrentCoordinate()[2];

        currentVelocity(3 * index) = node->getCurrentVelocity()[0];
        currentVelocity(3 * index + 1) = node->getCurrentVelocity()[1];
        currentVelocity(3 * index + 2) = node->getCurrentVelocity()[2];

        currentAcceleration(3 * index) = node->getCurrentAcceleration()[0];
        currentAcceleration(3 * index + 1) = node->getCurrentAcceleration()[1];
        currentAcceleration(3 * index + 2) = node->getCurrentAcceleration()[2];
    }

    Q = pastCoordinate / (beta * deltat * deltat) + pastVelocity / (beta * deltat) + (1 / (2 * beta) - 1) * pastAcceleration;
    R = pastVelocity + deltat * (1 - gamma) * pastAcceleration;

    inertialForce = prod(massMatrix, currentCoordinate) / (beta * deltat * deltat) - prod(massMatrix, Q) +
                    prod(amortecimento, currentCoordinate) * gamma / (beta * deltat) + prod(amortecimento, R) -
                    gamma * deltat * prod(amortecimento, Q);

    return inertialForce;
}

vector<double> Truss::ExternalForces()
{
    vector<double> externalForces(3 * nodes_.size(), 0.0);

    for (int i = 0; i < 3 * nodes_.size(); i++)
    {
        externalForces[i] = externalForces_[i];
    }

    return externalForces;
}

matrix<double> Truss::Hessian()
{
    matrix<double> hessian(3 * nodes_.size(), 3 * nodes_.size(), 0.0);

    for (int i = 0; i < elements_.size(); i++)
    {
        bounded_matrix<double, 6, 6> localHessian = elements_[i]->localHessian();
        int indexInitialNode = elements_[i]->getConnection()[0]->getIndex();
        int indexEndNode = elements_[i]->getConnection()[1]->getIndex();
        for (int ij = 0; ij < 3; ij++)
        {
            for (int ik = 0; ik < 3; ik++)
            {
                hessian(indexInitialNode * 3 + ij, indexInitialNode * 3 + ik) += localHessian(ij, ik);

                hessian(indexInitialNode * 3 + ij, indexEndNode * 3 + ik) += localHessian(ij, ik + 3);

                hessian(indexEndNode * 3 + ik, indexInitialNode * 3 + ij) += localHessian(ij + 3, ik);

                hessian(indexEndNode * 3 + ij, indexEndNode * 3 + ik) += localHessian(ij + 3, ik + 3);
            }
        }
    }

    return hessian;
}

matrix<double> Truss::TemperatureHessian(const int &numberOfSteps, const int &currentStep)
{
    matrix<double> hessian(3 * nodes_.size(), 3 * nodes_.size(), 0.0);

    for (int i = 0; i < elements_.size(); i++)
    {
        bounded_matrix<double, 6, 6> localHessian = elements_[i]->localTemperatureHessian(numberOfSteps, currentStep);
        int indexInitialNode = elements_[i]->getConnection()[0]->getIndex();
        int indexEndNode = elements_[i]->getConnection()[1]->getIndex();
        for (int ij = 0; ij < 3; ij++)
        {
            for (int ik = 0; ik < 3; ik++)
            {
                hessian(indexInitialNode * 3 + ij, indexInitialNode * 3 + ik) += localHessian(ij, ik);

                hessian(indexInitialNode * 3 + ij, indexEndNode * 3 + ik) += localHessian(ij, ik + 3);

                hessian(indexEndNode * 3 + ik, indexInitialNode * 3 + ij) += localHessian(ij + 3, ik);

                hessian(indexEndNode * 3 + ij, indexEndNode * 3 + ik) += localHessian(ij + 3, ik + 3);
            }
        }
    }
    return hessian;
}

int Truss::solveStaticProblem(const int &numberOfSteps, const double &tolerance)
{
    double normInitialCoordinate = 0.0;

    std::stringstream text;
    text << name_ << "-forçaXdeslocamento.txt";
    std::ofstream file(text.str());

    for (int i = 0; i < nodes_.size(); i++)
    {
        std::vector<double> initialCoordinate = nodes_[i]->getInitialCoordinate();
        normInitialCoordinate += initialCoordinate[0] * initialCoordinate[0] + initialCoordinate[1] * initialCoordinate[1] + initialCoordinate[2] * initialCoordinate[2];
    }

    //vector<double> auxiliarg(3 * nodes_.size(), 0.0);
    for (int loadStep = 0; loadStep <= numberOfSteps; loadStep++)
    {

        std::cout << "------------------------- LOAD STEP = "
                  << loadStep << " -------------------------\n";

        vector<double> dexternalForces = (loadStep)*ExternalForces() / numberOfSteps;

        for (int interation = 0; interation < 15; interation++) //definir o máximo de interações por passo de carga
        {
            vector<int> c(3 * nodes_.size(), 0.0);
            vector<double> g(3 * nodes_.size(), 0.0), deltaY(3 * nodes_.size(), 0.0);
            g = InternalForces() - TemperatureForces(numberOfSteps, loadStep) - dexternalForces;
            matrix<double, column_major> hessian = Hessian() + TemperatureHessian(numberOfSteps, loadStep);

            for (int i = 0; i < 3 * nodes_.size(); i++)
            {
                if (boundaryConditions_[i] == 1) //quando =1 é porque o deslocamento naquela direção está sendo
                {
                    for (int k = 0; k < 3 * nodes_.size(); k++)
                    {
                        hessian(i, k) = 0.0;
                        hessian(k, i) = 0.0;
                    }
                    hessian(i, i) = 1.0;
                    g[i] = 0.0;
                }
            }
            // vector<int> ipiv(3 * nodes_.size());
            // boost::numeric::bindings::lapack::getrf(hessian, ipiv);
            // boost::numeric::bindings::lapack::getri(hessian, ipiv, boost::numeric::bindings::lapack::optimal_workspace());
            // deltaY = g;
            // deltaY[7]=-2.8/7000;
            // g = -deltaY;
            // boost::numeric::bindings::lapack::gesv(hessian, c, g);

            deltaY = -g;

            boost::numeric::bindings::lapack::gesv(hessian, c, deltaY);

            double normDeltaY = 0.0;

            for (int ih = 0; ih < nodes_.size(); ih++) //loop para atualizar as coordenadas dos nós
            {
                int index = nodes_[ih]->getIndex();
                std::vector<double> currentCoordinate = nodes_[ih]->getCurrentCoordinate();

                normDeltaY += deltaY[3 * index] * deltaY[3 * index] + deltaY[3 * index + 1] * deltaY[3 * index + 1] + deltaY[3 * index + 2] * deltaY[3 * index + 2];

                currentCoordinate[0] += deltaY[3 * index];
                currentCoordinate[1] += deltaY[3 * index + 1];
                currentCoordinate[2] += deltaY[3 * index + 2];

                nodes_[ih]->setCurrentCoordinate(currentCoordinate);
            }

            //auxiliarg=InternalForces()-g;
            std::cout << "Iteration = " << interation
                      << "   x Norm = " << std::scientific << sqrt(normDeltaY / normInitialCoordinate)
                      << std::endl;

            if (sqrt(normDeltaY / normInitialCoordinate) <= tolerance)
                break;
        }

        exportToParaview(loadStep);

        // file //<< nodes_[0]->getCurrentCoordinate()[0] - nodes_[0]->getInitialCoordinate()[0] << " "
        //      << (nodes_[2]->getCurrentCoordinate()[1] - nodes_[2]->getInitialCoordinate()[1])*(-1.0) << " "
        //      << auxiliarg[7]*(-1.0) << std::endl;

        file //<< nodes_[0]->getCurrentCoordinate()[0] - nodes_[0]->getInitialCoordinate()[0] << " "
            << -1.0*(nodes_[2]->getCurrentCoordinate()[1] - nodes_[2]->getInitialCoordinate()[1]) << " "
            << -1.0*dexternalForces[7] << std::endl;
    }
}

int Truss::solveDynamicProblem(const int &numberOfTimes, const double &tolerance)
{
    double beta = 0.25, gamma = 0.5;
    double time = 0.0;

    matrix<double, column_major> hessian0(3 * nodes_.size(), 3 * nodes_.size());
    hessian0 = Hessian();
    matrix<double, column_major> mass(3 * nodes_.size(), 3 * nodes_.size());
    mass = MassMatrix();

    for (int i = 0; i < 3 * nodes_.size(); i++)
    {
        if (boundaryConditions_[i] == 1) //quando =1 é porque o deslocamento naquela direção está sendo
        {
            for (int k = 0; k < 3 * nodes_.size(); k++)
            {
                hessian0(i, k) = 0.0;
                hessian0(k, i) = 0.0;
                mass(i, k) = 0.0;
                mass(k, i) = 0.0;
            }
            hessian0(i, i) = 1.0;
            mass(i, i) = 1.0;
        }
    }

    vector<double> alphar(3 * nodes_.size()), alphai(3 * nodes_.size()), betinha(3 * nodes_.size());
    matrix<double, column_major> vl(1, 3 * nodes_.size());
    matrix<double, column_major> vr(3 * nodes_.size(), 3 * nodes_.size());
    double work[10 * (3 * nodes_.size())];

    boost::numeric::bindings::lapack::ggev('N', 'V', hessian0, mass, alphar, alphai, betinha, vl, vr, boost::numeric::bindings::lapack::optimal_workspace());

    vector<double> eigenvalues(3 * nodes_.size());
    for (int i = 0; i < 3 * nodes_.size(); i++)
    {
        eigenvalues(i) = alphar(i) / betinha(i);
    }

    for (int i = 0; i < 3 * nodes_.size(); i++)
    {
        eigenvalues(i) = sqrt(eigenvalues(i));
    }

    int auxiliar = int(3 * nodes_.size() / 20) + 2;
    vector<double> vecAuxiliar(auxiliar, 1000.0);
    vecAuxiliar(0) = -1000.0;
    //matrix<double> defmodo(auxiliar - 1, 3 * nodes_.size());

    for (int i = 0; i < (auxiliar - 1); i++)
    {
        for (int j = 0; j < 3 * nodes_.size(); j++)
        {
            if (eigenvalues(j) < vecAuxiliar(i + 1) and eigenvalues(j) > vecAuxiliar(i))
            {
                vecAuxiliar(i + 1) = eigenvalues(j);
                // for (int h = 0; h < 3 * nodes_.size(); h++)
                // {
                //     defmodo(i, h) = vr(h, j);
                // }
            }
        }
    }

    std::stringstream text2;
    text2 << name_ << "-ModosDeVibracao.txt";
    std::ofstream file2(text2.str());

    // for (int i = 0; i < (auxiliar - 1); i++)
    // {
    //     for (Node *n : nodes_)
    //     {
    //         int index = n->getIndex();
    //         std::vector<double> coord = n->getInitialCoordinate();
    //         coord[0] = coord[0] + defmodo(i, 3 * index)/10;
    //         coord[1] = coord[1] + defmodo(i, 3 * index + 1)/10;
    //         coord[2] = coord[2] + defmodo(i, 3 * index + 2)/10;
    //         n->setCurrentCoordinate(coord);
    //     }
    //     exportToParaview(i);
    // }
    // for (Node *n : nodes_)
    // {
    //     std::vector<double> coord = n->getInitialCoordinate();
    //     n->setCurrentCoordinate(coord);
    // }

    for (int i = 0; i < (vecAuxiliar.size() - 1); i++)
    {
        file2 << "w" << i << " = " << vecAuxiliar(i + 1) << " rad/(ms)" << std::endl;
    }

    const double deltat = 2.0 * 3.1415926535 / (10.0 * vecAuxiliar(auxiliar - 1)); //ntp=10
    file2 << std::endl;
    file2 << "deltat = " << deltat << " ms" << std::endl;

    vector<int> c(3 * nodes_.size(), 0);
    matrix<double, column_major> massMatrix = MassMatrix();
    vector<double> initialAcceleration(3 * nodes_.size(), 0.0);
    initialAcceleration = ExternalForces();

    for (int i = 0; i < 3 * nodes_.size(); i++)
    {
        if (boundaryConditions_[i] == 1) //quando =1 é porque o deslocamento naquela direção está sendo
        {
            for (int k = 0; k < 3 * nodes_.size(); k++)
            {
                massMatrix(i, k) = 0.0;
                massMatrix(k, i) = 0.0;
            }
            massMatrix(i, i) = 1.0;
            initialAcceleration(i) = 0.0;
        }
    }

    boost::numeric::bindings::lapack::gesv(massMatrix, c, initialAcceleration);

    for (Node *node : nodes_)
    {
        int index = node->getIndex();
        node->setCurrentAcceleration({initialAcceleration(3 * index), initialAcceleration(3 * index + 1), initialAcceleration(3 * index + 2)});
        node->setPastAcceleration({initialAcceleration(3 * index), initialAcceleration(3 * index + 1), initialAcceleration(3 * index + 2)});
    }

    std::stringstream text1;
    text1 << name_ << "-DeslocamentoxTempo.txt";
    std::ofstream file1(text1.str());

    double normInitialCoordinate = 0.0;
    for (int i = 0; i < nodes_.size(); i++)
    {
        std::vector<double> initialCoordinate = nodes_[i]->getInitialCoordinate();
        normInitialCoordinate += initialCoordinate[0] * initialCoordinate[0] + initialCoordinate[1] * initialCoordinate[1] + initialCoordinate[2] * initialCoordinate[2];
    }

    vector<double> dexternalForces = ExternalForces();

    for (int timeStep = 0; timeStep <= numberOfTimes; timeStep++)
    {
        time += deltat;

        std::cout << "------------------------- TIME STEP = "
                  << timeStep << " -------------------------\n";

        for (int interation = 0; interation < 20; interation++) //definir o máximo de interações por passo de carga
        {
            vector<int> c(3 * nodes_.size(), 0);
            vector<double> g(3 * nodes_.size(), 0.0), deltaY(3 * nodes_.size(), 0.0);

            g = InternalForces() + InertialForces(beta, gamma, deltat) - dexternalForces;
            //amortemciento
            matrix<double, column_major> hessian = Hessian() + MassMatrix() / (beta * deltat * deltat) + 0.05 * MassMatrix() * gamma / (beta * deltat);

            for (int i = 0; i < 3 * nodes_.size(); i++)
            {
                if (boundaryConditions_[i] == 1) //quando =1 é porque o deslocamento naquela direção está sendo
                {
                    for (int k = 0; k < 3 * nodes_.size(); k++)
                    {
                        hessian(i, k) = 0.0;
                        hessian(k, i) = 0.0;
                    }
                    hessian(i, i) = 1.0;
                    g[i] = 0.0;
                }
            }

            deltaY = -g;

            boost::numeric::bindings::lapack::gesv(hessian, c, deltaY);

            double normDeltaY = 0.0;

            for (int h = 0; h < nodes_.size(); h++) //loop para atualizar as coordenadas dos nós
            {
                int index = nodes_[h]->getIndex();
                std::vector<double> currentCoordinate = nodes_[h]->getCurrentCoordinate();

                normDeltaY += deltaY[3 * index] * deltaY[3 * index] + deltaY[3 * index + 1] * deltaY[3 * index + 1] + deltaY[3 * index + 2] * deltaY[3 * index + 2];

                currentCoordinate[0] += deltaY[3 * index];
                currentCoordinate[1] += deltaY[3 * index + 1];
                currentCoordinate[2] += deltaY[3 * index + 2];
                nodes_[h]->setCurrentCoordinate(currentCoordinate);

                std::vector<double> accel;
                accel.reserve(3);

                for (int i = 0; i < 3; i++)
                {
                    accel[i] = (nodes_[h]->getCurrentCoordinate()[i]) / (beta * deltat * deltat) - (nodes_[h]->getPastCoordinate()[i]) / (beta * deltat * deltat) -
                               (nodes_[h]->getPastVelocity()[i]) / (beta * deltat) - (nodes_[h]->getPastAcceleration()[i]) * (0.5 / beta - 1.0);
                }
                nodes_[h]->setCurrentAcceleration({accel[0], accel[1], accel[2]});

                std::vector<double> vel;
                vel.reserve(3);

                for (int i = 0; i < 3; i++)
                {
                    vel[i] = (nodes_[h]->getCurrentAcceleration()[i]) * gamma * deltat + (nodes_[h]->getPastVelocity()[i]) * 1.0 +
                             (nodes_[h]->getPastAcceleration()[i]) * deltat * (1.0 - gamma);
                }
                nodes_[h]->setCurrentVelocity({vel[0], vel[1], vel[2]});
            }

            std::cout << "Iteration = " << interation
                      << "   x Norm = " << std::scientific << sqrt(normDeltaY / normInitialCoordinate)
                      << std::endl;

            if (sqrt(normDeltaY / normInitialCoordinate) <= tolerance)
                break;
        }
        exportToParaview(timeStep);

        for (Node *node : nodes_)
        {
            std::vector<double> updatingCoordinate = node->getCurrentCoordinate();
            node->setPastCoordinate(updatingCoordinate);
            std::vector<double> updatingVel = node->getCurrentVelocity();
            node->setPastVelocity(updatingVel);
            std::vector<double> updatingAccel = node->getCurrentAcceleration();
            node->setPastAcceleration(updatingAccel);
        }

        file1 << nodes_[23]->getCurrentCoordinate()[0] - nodes_[23]->getInitialCoordinate()[0] << " " << time << std::endl;
    }
}

void Truss::exportToParaview(const int &loadstep)
{
    std::stringstream text;
    text << name_ << loadstep << ".vtu";
    std::ofstream file(text.str());

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << nodes_.size()
         << "\"  NumberOfCells=\"" << elements_.size()
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
        file << n->getCurrentCoordinate()[0] << " " << n->getCurrentCoordinate()[1] << " " << n->getCurrentCoordinate()[2] << "\n";
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

    for (Element *el : elements_)
    {

        file << el->getConnection()[0]->getIndex() << " " << el->getConnection()[1]->getIndex() << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    int aux = 0;

    for (Element *el : elements_)
    {
        int n = el->getConnection().size();
        aux += n;
        file << aux << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (Element *el : elements_)
    {
        file << 4 << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";
    //nodal results
    file << "    <PointData>"
         << "\n";
    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         << "Name=\"Displacement\" format=\"ascii\">"
         << "\n";

    for (Node *n : nodes_)
    {
        file << n->getCurrentCoordinate()[0] - n->getInitialCoordinate()[0] << " "
             << n->getCurrentCoordinate()[1] - n->getInitialCoordinate()[1] << " "
             << n->getCurrentCoordinate()[2] - n->getInitialCoordinate()[2] << "\n";
    }

    file << "      </DataArray> "
         << "\n";
    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         << "Name=\"Velocity\" format=\"ascii\">"
         << "\n";

    for (Node *n : nodes_)
    {
        file << n->getCurrentVelocity()[0] << " " << n->getCurrentVelocity()[1] << " " << n->getCurrentVelocity()[2] << "\n";
    }
    file << "      </DataArray> "
         << "\n";
    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         << "Name=\"Acceleration\" format=\"ascii\">"
         << "\n";

    for (Node *n : nodes_)
    {
        file << n->getCurrentAcceleration()[0] << " " << n->getCurrentAcceleration()[1] << " " << n->getCurrentAcceleration()[2] << "\n";
    }
    file << "      </DataArray> "
         << "\n";
    file << "    </PointData>"
         << "\n";
    //elemental results
    file << "    <CellData>"
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

void Truss::readInput(const std::string &read, const std::string &typeAnalyze)
{
    std::ifstream file(read);
    std::string line;

    file >> name_;

    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);

    int nmaterial, nnode, nelement, nforce, ntemp, nboundary, numberOfTimesOrSteps_;
    double deltat_, tolerance_;

    file >> numberOfTimesOrSteps_ >> tolerance_;

    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);

    file >> nmaterial;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nmaterial; i++)
    {
        int index;
        double young, plastStrain, hardeningModulus, density, expansion;

        file >> index >> young >> plastStrain >> hardeningModulus >> density >> expansion;

        addMaterial(index, young, plastStrain, hardeningModulus, density, expansion);

        std::getline(file, line);
    }

    std::getline(file, line);
    std::getline(file, line);

    file >> nnode;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nnode; i++)
    {
        int index;
        double x1, x2, x3;

        file >> index >> x1 >> x2 >> x3;

        addNode(index, {x1, x2, x3});

        std::getline(file, line);
    }

    std::getline(file, line);
    std::getline(file, line);

    file >> nelement;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nelement; i++)
    {
        int index, initialNode, endNode, material;
        double area;

        file >> index >> initialNode >> endNode >> material >> area;

        addElement(index, {initialNode, endNode}, material, area);

        std::getline(file, line);
    }

    std::getline(file, line);
    std::getline(file, line);

    file >> nforce;

    std::getline(file, line);
    std::getline(file, line);

    int cont = nodes_.size();
    boundaryConditions_.reserve(3 * cont);
    externalForces_.reserve(3 * cont);

    for (int i = 0; i < (3 * cont); i++)
    {
        externalForces_[i] = 0.0;
        boundaryConditions_[i] = 0;
    }

    for (int i = 0; i < nforce; i++)
    {
        int index;
        double x1, x2, x3;

        file >> index >> x1 >> x2 >> x3;

        externalForces_[3 * index] = x1;
        externalForces_[3 * index + 1] = x2;
        externalForces_[3 * index + 2] = x3;

        std::getline(file, line);
    }

    std::getline(file, line);
    std::getline(file, line);

    file >> ntemp;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < ntemp; i++)
    {
        int index;
        double deltaTemp;

        file >> index >> deltaTemp;

        elements_[index]->setDeltaTemp(deltaTemp);

        std::getline(file, line);
    }

    std::getline(file, line);
    std::getline(file, line);

    file >> nboundary;

    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < nboundary; i++)
    {
        int index;
        double x1, x2, x3;

        file >> index >> x1 >> x2 >> x3;

        boundaryConditions_[3 * index] = x1;
        boundaryConditions_[3 * index + 1] = x2;
        boundaryConditions_[3 * index + 2] = x3;

        std::getline(file, line);
    }

    if (typeAnalyze == "Static")
    {
        solveStaticProblem(numberOfTimesOrSteps_, tolerance_);
    }

    if (typeAnalyze == "Dynamic")
    {
        solveDynamicProblem(numberOfTimesOrSteps_, tolerance_);
    }
}
