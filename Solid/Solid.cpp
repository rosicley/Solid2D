#include "Solid.h"

Solid::Solid()
{
    deltat_=1.0;
    gamma_=0.5;
    beta_=0.25;
    shapeForces_(0)=0.0;
    shapeForces_(1)=0.0;
}

Solid::~Solid() {}

void Solid::addMaterial(const int &index,
                        const double &young,
                        const double &poisson,
                        const double &density)
{
    Material *mat = new Material(index, young, poisson, density);
    materials_.push_back(mat);
}

void Solid::addFiberMaterial(const int &index,
                             const double &young,
                             const double &plastStrain,
                             const double &hardeningModulus,
                             const double &density)
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

void Solid::addElement(const int &index,
                       const std::vector<int> &nodesIndex,
                       const int &materialIndex,
                       const double &thickness,
                       const std::string &elementType)
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

void Solid::addFiberElement(const int &index,
                            const std::vector<int> &connection,
                            const int &material,
                            const double &area)
{
    std::vector<FiberNode *> fiberNodes;
    fiberNodes.push_back(fiberNodes_[connection[0]]);
    fiberNodes.push_back(fiberNodes_[connection[1]]);
    FiberElement *el = new FiberElement(index, fiberNodes, fiberMaterials_[material], area);
    fiberElements_.push_back(el);
}

void Solid::addDirichletCondition(const int &index,
                                  const int &direction,
                                  const double &value)
{
    DirichletCondition *cond = new DirichletCondition(nodes_[index], direction, value);
    dirichletConditions_.push_back(cond);
}

void Solid::addNeumannCondition(const int &index,
                                const int &direction,
                                const double &value)
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

std::pair<vector<double>, matrix<double, column_major>> Solid::globalSolid(const std::string &typeAnalyze, const int &step, const int &numberOfStep)
{
    vector<double> first(2 * nodes_.size(), 0.0);
    matrix<double, column_major> second(2 * nodes_.size(), 2 * nodes_.size(), 0.0);
    int n = (order_ + 1) * (order_ + 2) / 2.0;

    for (Element *el : elements_)
    {
        std::pair<vector<double>, matrix<double>> local;

        local = el->elementContributions(planeState_, typeAnalyze, step, numberOfStep);

        vector<int> indexNodes(n);

        for (int i = 0; i < n; i++)
        {
            indexNodes(i) = el->getConnection()[i]->getIndex();
        }

        for (int in = 0; in < n; in++)
        {

            first(2 * indexNodes(in)) += local.first(2 * in);
            first(2 * indexNodes(in) + 1) += local.first(2 * in + 1);

            for (int jn = 0; jn < n; jn++)
            {
                second(2 * indexNodes(in), 2 * indexNodes(jn)) += local.second(2 * in, 2 * jn);
                second(2 * indexNodes(in), 2 * indexNodes(jn) + 1) += local.second(2 * in, 2 * jn + 1);
                second(2 * indexNodes(in) + 1, 2 * indexNodes(jn)) += local.second(2 * in + 1, 2 * jn);
                second(2 * indexNodes(in) + 1, 2 * indexNodes(jn) + 1) += local.second(2 * in + 1, 2 * jn + 1);
            }
        }
    }

    if (fiberElements_.size() >= 1)
    {
        for (FiberElement *fib : fiberInsideSolid_)
        {
            std::pair<vector<double>, matrix<double>> local;
            local = fib->fiberContributions(typeAnalyze, deltat_, beta_);

            FiberNode *initial = fib->getConnection()[0];
            FiberNode *end = fib->getConnection()[1];

            int indexSolidInitical = initial->getIndexSolidElement();
            int indexSolidEnd = end->getIndexSolidElement();

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
                mfibra((2 * n + 1) + 2 * i, 3) = phi_end(i);
            }

            matrix<double> kfc(4 * n, 4 * n, 0.0);
            matrix<double> aux(4 * n, 4, 0.0);
            aux = prod(mfibra, local.second);
            kfc = prod(aux, trans(mfibra));

            vector<double> fc(4 * n, 0.0);
            fc = prod(mfibra, local.first);

            vector<int> indexNodes(2 * n);

            for (int i = 0; i < n; i++)
            {
                indexNodes(i) = elements_[indexSolidInitical]->getConnection()[i]->getIndex();
                indexNodes(i + n) = elements_[indexSolidEnd]->getConnection()[i]->getIndex();
            }

            for (int in = 0; in < n; in++)
            {
                first(2 * indexNodes(in)) += fc(2 * in);
                first(2 * indexNodes(in) + 1) += fc(2 * in + 1);

                first(2 * indexNodes(in + n)) += fc(2 * (in + n));
                first(2 * indexNodes(in + n) + 1) += fc(2 * (in + n) + 1);

                for (int jn = 0; jn < n; jn++)
                {
                    second(2 * indexNodes(in), 2 * indexNodes(jn)) += kfc(2 * in, 2 * jn);
                    second(2 * indexNodes(in), 2 * indexNodes(jn) + 1) += kfc(2 * in, 2 * jn + 1);
                    second(2 * indexNodes(in) + 1, 2 * indexNodes(jn)) += kfc(2 * in + 1, 2 * jn);
                    second(2 * indexNodes(in) + 1, 2 * indexNodes(jn) + 1) += kfc(2 * in + 1, 2 * jn + 1);

                    second(2 * indexNodes(in + n), 2 * indexNodes(jn + n)) += kfc(2 * (in + n), 2 * (jn + n));
                    second(2 * indexNodes(in + n), 2 * indexNodes(jn + n) + 1) += kfc(2 * (in + n), 2 * (jn + n) + 1);
                    second(2 * indexNodes(in + n) + 1, 2 * indexNodes(jn + n)) += kfc(2 * (in + n) + 1, 2 * (jn + n));
                    second(2 * indexNodes(in + n) + 1, 2 * indexNodes(jn + n) + 1) += kfc(2 * (in + n) + 1, 2 * (jn + n) + 1);
                }
            }
        }
    }

    return std::make_pair(first, second);
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

int Solid::solveStaticProblem(const int &numberOfSteps, const int &maximumOfInteration, const double &tolerance)
{
    // double normInitialCoordinate = 0.0;

    // for (int i = 0; i < nodes_.size(); i++)
    // {
    //     bounded_vector<double, 2> initialCoordinate = nodes_[i]->getInitialCoordinate();
    //     normInitialCoordinate += initialCoordinate(0) * initialCoordinate(0) + initialCoordinate(1) * initialCoordinate(1);
    // }

    //normInitialCoordinate = sqrt(normInitialCoordinate);

    vector<double> externalForces(2 * nodes_.size(), 0.0);
    externalForces = ExternalForces();
    int n = (order_ + 1) * (order_ + 2) / 2.0;

    exportToParaview(0);

    for (int loadStep = 1; loadStep <= numberOfSteps; loadStep++)
    {

        std::cout << "------------------------- LOAD STEP = "
                  << loadStep << " -------------------------\n";

        for (int interation = 0; interation < maximumOfInteration; interation++) //definir o máximo de interações por passo de carga
        {
            std::pair<vector<double>, matrix<double, column_major>> variationOfPI = globalSolid("STATIC", loadStep, numberOfSteps);

            variationOfPI.first = variationOfPI.first + ((1.0 * loadStep) / (1.0 * numberOfSteps)) * externalForces;

            for (DirichletCondition *cond : dirichletConditions_)
            {
                int indexNode = (cond->getNode())->getIndex();
                int direction = cond->getDirection();
                double value;
                if (interation == 0)
                {
                    value = (cond->getValue()) / numberOfSteps;
                }
                else
                {
                    value = 0.0;
                }

                for (int j = 0; j < 2 * nodes_.size(); j++)
                {
                    variationOfPI.first(j) += -value * (variationOfPI.second(j, 2 * indexNode + direction));
                    variationOfPI.second(j, 2 * indexNode + direction) = 0.0;
                    variationOfPI.second(2 * indexNode + direction, j) = 0.0;
                    variationOfPI.second(2 * indexNode + direction, 2 * indexNode + direction) = 1.0;
                    variationOfPI.first(2 * indexNode + direction) = value;
                }
            }

            //boost::numeric::bindings::lapack::posv(variationOfPI.second, variationOfPI.first);
            vector<int> c(2 * nodes_.size(), 0.0);
            boost::numeric::bindings::lapack::gesv(variationOfPI.second, c, variationOfPI.first);

            for (int ih = 0; ih < nodes_.size(); ih++) //loop para atualizar as coordenadas dos nós
            {
                int index = nodes_[ih]->getIndex();
                bounded_vector<double, 2> currentCoordinate = nodes_[ih]->getCurrentCoordinate();

                currentCoordinate(0) += variationOfPI.first(2 * index);
                currentCoordinate(1) += variationOfPI.first(2 * index + 1);

                nodes_[ih]->setCurrentCoordinate(currentCoordinate);
            }

            for (FiberElement *fiber : fiberInsideSolid_)
            {
                for (FiberNode *fnode : fiber->getConnection())
                {
                    Element *solid = elements_[fnode->getIndexSolidElement()];
                    bounded_vector<double, 2> dimensionless = fnode->getDimensionlessCoordinates();
                    vector<double> phi(n, 0.0);
                    phi = domainShapeFunction(dimensionless(0), dimensionless(1));
                    vector<double> correction_X1(n, 0.0);
                    vector<double> correction_X2(n, 0.0);

                    int auxiliar = 0;
                    for (Node *solidNode : solid->getConnection())
                    {
                        int index = solidNode->getIndex();
                        correction_X1(auxiliar) = variationOfPI.first(2 * index);
                        correction_X2(auxiliar) = variationOfPI.first(2 * index + 1);
                        auxiliar = auxiliar + 1;
                    }

                    bounded_vector<double, 2> currentCoordinate;
                    currentCoordinate = fnode->getCurrentCoordinate();
                    currentCoordinate(0) = currentCoordinate(0) + inner_prod(phi, correction_X1);
                    currentCoordinate(1) = currentCoordinate(1) + inner_prod(phi, correction_X2);

                    fnode->setCurrentCoordinate(currentCoordinate);
                }
                fiber->updateNormalForce();
            }

            double error = (norm_2(variationOfPI.first)); /// normInitialCoordinate;

            std::cout << "Iteration = " << interation
                      << "   x Norm = " << std::scientific << error
                      << std::endl;

            if (error <= tolerance)
                break;
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        file << fiber->getConnection()[0]->getCurrentCoordinate()(0) << " " << fiber->getConnection()[0]->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
        file << fiber->getConnection()[1]->getCurrentCoordinate()(0) << " " << fiber->getConnection()[1]->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
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

    int initial, end;
    initial = nodes_.size() - 1;
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        initial = initial + 1;
        end = initial + 1;

        file << initial << " " << end << "\n";

        initial = end;
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        aux += 2;
        file << aux << "\n";
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        file << 3 << "\n";
    }
    //ATÉ AQUI ESTÁ PARA FIBRA!!!
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        bounded_vector<double, 2> initial = fiber->getConnection()[0]->getInitialCoordinate();
        bounded_vector<double, 2> current = fiber->getConnection()[0]->getCurrentCoordinate();
        file << current(0) - initial(0) << " " << current(1) - initial(1) << "\n";

        initial = fiber->getConnection()[1]->getInitialCoordinate();
        current = fiber->getConnection()[1]->getCurrentCoordinate();
        file << current(0) - initial(0) << " " << current(1) - initial(1) << "\n";
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        bounded_vector<double, 2> currentVelocity = fiber->getConnection()[0]->getCurrentVelocity();
        file << currentVelocity(0) << " " << currentVelocity(1) << "\n";

        currentVelocity = fiber->getConnection()[1]->getCurrentVelocity();
        file << currentVelocity(0) << " " << currentVelocity(1) << "\n";
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        bounded_vector<double, 2> current = fiber->getConnection()[0]->getCurrentAcceleration();
        file << current(0) << " " << current(1) << "\n";

        current = fiber->getConnection()[1]->getCurrentAcceleration();
        file << current(0) << " " << current(1) << "\n";
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        file << 0.0 << " " << 0.0 << "\n";
        file << 0.0 << " " << 0.0 << "\n";
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        file << 0.0 << "\n";
        file << 0.0 << "\n";
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
    for (FiberElement *fiber : fiberInsideSolid_)
    {
        file << fiber->getConnection()[0]->getNormalForce() << "\n";
        file << fiber->getConnection()[1]->getNormalForce() << "\n";
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
        shapeForces_(0)=b1;
        shapeForces_(1)=b2;

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

    // int teste = elements_[255]->getConnection()[4]->getIndex();
    // //bounded_vector<double, 2> teste = nodes_[1000]->getCurrentCoordinate();
    // std::cout<<elements_[147]-;
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
    incidenceOfFibers();
}

void Solid::incidenceOfFibers()
{
    matrix<double> dataelems(elements_.size(), 3, 0.0);

    int n = (order_ + 1) * (order_ + 2) / 2.0;

    int cont = 0;
    for (Element *el : elements_)
    {
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
        cont = cont + 1;
    }
    int no = 0;
    for (FiberNode *fnode : fiberNodes_)
    {
        cont = 0;
        for (Element *el : elements_)
        {
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

                    ponto0(0) = inner_prod((el->domainShapeFunction(ksi, eta)), nos_X1);
                    ponto0(1) = inner_prod((el->domainShapeFunction(ksi, eta)), nos_X2);

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
            cont = cont + 1;
        }
    }

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