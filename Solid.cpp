#include "Solid.h"

// Solid::Solid()
// {
//     deltat_=1.0;
//     gamma_=0.5;
//     beta_=0.25;
//     shapeForce_=0.0;
// }

// Solid::~Solid() {}

void Solid::addMaterial(const int &index,
                        const double &young,
                        const double &poisson,
                        const double &density)
{
    Material *mat = new Material(index, young, poisson, density);
    materials_.push_back(mat);
}

void Solid::addNode(const int &index, const bounded_vector<double, 2> &initialCoordinate)
{
    Node *n = new Node(index, initialCoordinate);
    nodes_.push_back(n);
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

    e->setShapeForce({shapeForces_[0], shapeForces_[1]});

    elements_.push_back(e);
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

        //ERRRO AQUI:
        local = el->elementContributions(planeState_, typeAnalyze, step, numberOfStep);

        vector<int> indexNodes(n);

        for (int i = 0; i < n; i++)
        {
            indexNodes(i) = el->getConnection()[i]->getIndex();
        }

        for (int in = 0; in < n; in++)
        {
            for (int id = 0; id <= 1; id++)
            {
                first(2 * indexNodes(in) + id) += local.first(2 * in + id);

                for (int jn = 0; jn < n; jn++)
                {
                    for (int jd = 0; jd <= 1; jd++)
                    {
                        second(2 * indexNodes(in) + id, 2 * indexNodes(jn) + jd) += local.second(2 * in + id, 2 * jn + jd);
                    }
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

            //boost::numeric::bindings::lapack::posv(variationOfPI.second, variationOfPI.first);

            for (int ih = 0; ih < nodes_.size(); ih++) //loop para atualizar as coordenadas dos nós
            {
                int index = nodes_[ih]->getIndex();
                bounded_vector<double, 2> currentCoordinate = nodes_[ih]->getCurrentCoordinate();

                currentCoordinate(0) += variationOfPI.first(2 * index);
                currentCoordinate(1) += variationOfPI.first(2 * index + 1);

                nodes_[ih]->setCurrentCoordinate(currentCoordinate);
            }

            double error = (norm_2(variationOfPI.first)); /// normInitialCoordinate;

            std::cout << "Iteration = " << interation
                      << "   x Norm = " << std::scientific << error
                      << std::endl;

            if (error <= tolerance)
                break;
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
        file << n->getCurrentCoordinate()(0) << " " << n->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
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
    // for (int i=0; i<nodes_.size(); i++)
    // {
    //     file << nodes_[i]->getCurrentCoordinate()(0) - nodes_[i]->getInitialCoordinate()(0) << " "
    //          << nodes_[i]->getCurrentCoordinate()(1) - nodes_[i]->getInitialCoordinate()(1) << "\n";
    // }
    // {
    //     file << n->getCurrentCoordinate()(0) - n->getInitialCoordinate()(0) << " "
    //          << n->getCurrentCoordinate()(1) - n->getInitialCoordinate()(1) << "\n";
    // }
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
    file << "      </DataArray> "
         << "\n";
    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Acceleration\" format=\"ascii\">"
         << "\n";

    for (Node *n : nodes_)
    {
        file << n->getCurrentAcceleration()(0) << " " << n->getCurrentAcceleration()(1) << "\n";
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

    shapeForces_.reserve(2);

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

        file >> material >> thickness >> shapeForces_[0] >> shapeForces_[1]; //shapeForce_(0) >> shapeForce_(1);

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

    // std::cout<<"teste";
    // int teste = elements_[255]->getConnection()[4]->getIndex();
    // //bounded_vector<double, 2> teste = nodes_[1000]->getCurrentCoordinate();
    // std::cout<<elements_[147]-;
}
