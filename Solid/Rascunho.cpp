// //vector<double> deltaY(2 * nodes_.size(), 0.0);

// // deltaY = variationOfPI.first + (loadStep / numberOfSteps) * externalForces;

// // for (DirichletCondition *cond : dirichletConditions_)
// // {
// //     int indexNode = cond->getNode()->getIndex();
// //     int direction = cond->getDirection();
// //     double value = (cond->getValue()) / numberOfSteps;

// //     for (int j = 0; j < 2 * nodes_.size(); j++)
// //     {
// //         deltaY(j) += -value*(variationOfPI.second(j,2*indexNode+direction));
// //         variationOfPI.second(j,2*indexNode+direction)=0.0;
// //         variationOfPI.second(2*indexNode+direction,j)=0.0;
// //         variationOfPI.second(2*indexNode+direction,2*indexNode+direction)=1.0;
// //         deltaY(2*indexNode+direction)=value;
// //     }
// // }

// // double normDeltaY = norm_2(deltaY);

// // for (int ih = 0; ih < nodes_.size(); ih++) //loop para atualizar as coordenadas dos nós
// // {
// //     int index = nodes_[ih]->getIndex();
// //     bounded_vector<double, 2> currentCoordinate = nodes_[ih]->getCurrentCoordinate();

// //     currentCoordinate[0] += deltaY[3 * index];
// //     currentCoordinate[1] += deltaY[3 * index + 1];
// //     currentCoordinate[2] += deltaY[3 * index + 2];

// //     nodes_[ih]->setCurrentCoordinate(currentCoordinate);
// // }

// //auxiliarg=InternalForces()-g;



// // int Solid::solveDynamicProblem(const int &numberOfTimes, const double &tolerance)
// // {
// //     double beta = 0.25, gamma = 0.5;
// //     double time = 0.0;

// //     matrix<double, column_major> hessian0(3 * nodes_.size(), 3 * nodes_.size());
// //     hessian0 = Hessian();
// //     matrix<double, column_major> mass(3 * nodes_.size(), 3 * nodes_.size());
// //     mass = MassMatrix();

// //     for (int i = 0; i < 3 * nodes_.size(); i++)
// //     {
// //         if (boundaryConditions_[i] == 1) //quando =1 é porque o deslocamento naquela direção está sendo
// //         {
// //             for (int k = 0; k < 3 * nodes_.size(); k++)
// //             {
// //                 hessian0(i, k) = 0.0;
// //                 hessian0(k, i) = 0.0;
// //                 mass(i, k) = 0.0;
// //                 mass(k, i) = 0.0;
// //             }
// //             hessian0(i, i) = 1.0;
// //             mass(i, i) = 1.0;
// //         }
// //     }

// //     vector<double> alphar(3 * nodes_.size()), alphai(3 * nodes_.size()), betinha(3 * nodes_.size());
// //     matrix<double, column_major> vl(1, 3 * nodes_.size());
// //     matrix<double, column_major> vr(3 * nodes_.size(), 3 * nodes_.size());
// //     double work[10 * (3 * nodes_.size())];

// //     boost::numeric::bindings::lapack::ggev('N', 'V', hessian0, mass, alphar, alphai, betinha, vl, vr, boost::numeric::bindings::lapack::optimal_workspace());

// //     vector<double> eigenvalues(3 * nodes_.size());
// //     for (int i = 0; i < 3 * nodes_.size(); i++)
// //     {
// //         eigenvalues(i) = alphar(i) / betinha(i);
// //     }

// //     for (int i = 0; i < 3 * nodes_.size(); i++)
// //     {
// //         eigenvalues(i) = sqrt(eigenvalues(i));
// //     }

// //     int auxiliar = int(3 * nodes_.size() / 20) + 2;
// //     vector<double> vecAuxiliar(auxiliar, 1000.0);
// //     vecAuxiliar(0) = -1000.0;
// //     //matrix<double> defmodo(auxiliar - 1, 3 * nodes_.size());

// //     for (int i = 0; i < (auxiliar - 1); i++)
// //     {
// //         for (int j = 0; j < 3 * nodes_.size(); j++)
// //         {
// //             if (eigenvalues(j) < vecAuxiliar(i + 1) and eigenvalues(j) > vecAuxiliar(i))
// //             {
// //                 vecAuxiliar(i + 1) = eigenvalues(j);
// //                 // for (int h = 0; h < 3 * nodes_.size(); h++)
// //                 // {
// //                 //     defmodo(i, h) = vr(h, j);
// //                 // }
// //             }
// //         }
// //     }

// //     std::stringstream text2;
// //     text2 << name_ << "-ModosDeVibracao.txt";
// //     std::ofstream file2(text2.str());

// //     // for (int i = 0; i < (auxiliar - 1); i++)
// //     // {
// //     //     for (Node *n : nodes_)
// //     //     {
// //     //         int index = n->getIndex();
// //     //         std::vector<double> coord = n->getInitialCoordinate();
// //     //         coord[0] = coord[0] + defmodo(i, 3 * index)/10;
// //     //         coord[1] = coord[1] + defmodo(i, 3 * index + 1)/10;
// //     //         coord[2] = coord[2] + defmodo(i, 3 * index + 2)/10;
// //     //         n->setCurrentCoordinate(coord);
// //     //     }
// //     //     exportToParaview(i);
// //     // }
// //     // for (Node *n : nodes_)
// //     // {
// //     //     std::vector<double> coord = n->getInitialCoordinate();
// //     //     n->setCurrentCoordinate(coord);
// //     // }

// //     for (int i = 0; i < (vecAuxiliar.size() - 1); i++)
// //     {
// //         file2 << "w" << i << " = " << vecAuxiliar(i + 1) << " rad/(ms)" << std::endl;
// //     }

// //     const double deltat = 2.0 * 3.1415926535 / (10.0 * vecAuxiliar(auxiliar - 1)); //ntp=10
// //     file2 << std::endl;
// //     file2 << "deltat = " << deltat << " ms" << std::endl;

// //     vector<int> c(3 * nodes_.size(), 0);
// //     matrix<double, column_major> massMatrix = MassMatrix();
// //     vector<double> initialAcceleration(3 * nodes_.size(), 0.0);
// //     initialAcceleration = ExternalForces();

// //     for (int i = 0; i < 3 * nodes_.size(); i++)
// //     {
// //         if (boundaryConditions_[i] == 1) //quando =1 é porque o deslocamento naquela direção está sendo
// //         {
// //             for (int k = 0; k < 3 * nodes_.size(); k++)
// //             {
// //                 massMatrix(i, k) = 0.0;
// //                 massMatrix(k, i) = 0.0;
// //             }
// //             massMatrix(i, i) = 1.0;
// //             initialAcceleration(i) = 0.0;
// //         }
// //     }

// //     boost::numeric::bindings::lapack::gesv(massMatrix, c, initialAcceleration);

// //     for (Node *node : nodes_)
// //     {
// //         int index = node->getIndex();
// //         node->setCurrentAcceleration({initialAcceleration(3 * index), initialAcceleration(3 * index + 1), initialAcceleration(3 * index + 2)});
// //         node->setPastAcceleration({initialAcceleration(3 * index), initialAcceleration(3 * index + 1), initialAcceleration(3 * index + 2)});
// //     }

// //     std::stringstream text1;
// //     text1 << name_ << "-DeslocamentoxTempo.txt";
// //     std::ofstream file1(text1.str());

// //     double normInitialCoordinate = 0.0;
// //     for (int i = 0; i < nodes_.size(); i++)
// //     {
// //         std::vector<double> initialCoordinate = nodes_[i]->getInitialCoordinate();
// //         normInitialCoordinate += initialCoordinate[0] * initialCoordinate[0] + initialCoordinate[1] * initialCoordinate[1] + initialCoordinate[2] * initialCoordinate[2];
// //     }

// //     vector<double> dexternalForces = ExternalForces();

// //     for (int timeStep = 0; timeStep <= numberOfTimes; timeStep++)
// //     {
// //         time += deltat;

// //         std::cout << "------------------------- TIME STEP = "
// //                   << timeStep << " -------------------------\n";

// //         for (int interation = 0; interation < 20; interation++) //definir o máximo de interações por passo de carga
// //         {
// //             vector<int> c(3 * nodes_.size(), 0);
// //             vector<double> g(3 * nodes_.size(), 0.0), deltaY(3 * nodes_.size(), 0.0);

// //             g = InternalForces() + InertialForces(beta, gamma, deltat) - dexternalForces;
// //             //amortemciento
// //             matrix<double, column_major> hessian = Hessian() + MassMatrix() / (beta * deltat * deltat) + 0.05 * MassMatrix() * gamma / (beta * deltat);

// //             for (int i = 0; i < 3 * nodes_.size(); i++)
// //             {
// //                 if (boundaryConditions_[i] == 1) //quando =1 é porque o deslocamento naquela direção está sendo
// //                 {
// //                     for (int k = 0; k < 3 * nodes_.size(); k++)
// //                     {
// //                         hessian(i, k) = 0.0;
// //                         hessian(k, i) = 0.0;
// //                     }
// //                     hessian(i, i) = 1.0;
// //                     g[i] = 0.0;
// //                 }
// //             }

// //             deltaY = -g;

// //             boost::numeric::bindings::lapack::gesv(hessian, c, deltaY);

// //             double normDeltaY = 0.0;

// //             for (int h = 0; h < nodes_.size(); h++) //loop para atualizar as coordenadas dos nós
// //             {
// //                 int index = nodes_[h]->getIndex();
// //                 std::vector<double> currentCoordinate = nodes_[h]->getCurrentCoordinate();

// //                 normDeltaY += deltaY[3 * index] * deltaY[3 * index] + deltaY[3 * index + 1] * deltaY[3 * index + 1] + deltaY[3 * index + 2] * deltaY[3 * index + 2];

// //                 currentCoordinate[0] += deltaY[3 * index];
// //                 currentCoordinate[1] += deltaY[3 * index + 1];
// //                 currentCoordinate[2] += deltaY[3 * index + 2];
// //                 nodes_[h]->setCurrentCoordinate(currentCoordinate);

// //                 std::vector<double> accel;
// //                 accel.reserve(3);

// //                 for (int i = 0; i < 3; i++)
// //                 {
// //                     accel[i] = (nodes_[h]->getCurrentCoordinate()[i]) / (beta * deltat * deltat) - (nodes_[h]->getPastCoordinate()[i]) / (beta * deltat * deltat) -
// //                                (nodes_[h]->getPastVelocity()[i]) / (beta * deltat) - (nodes_[h]->getPastAcceleration()[i]) * (0.5 / beta - 1.0);
// //                 }
// //                 nodes_[h]->setCurrentAcceleration({accel[0], accel[1], accel[2]});

// //                 std::vector<double> vel;
// //                 vel.reserve(3);

// //                 for (int i = 0; i < 3; i++)
// //                 {
// //                     vel[i] = (nodes_[h]->getCurrentAcceleration()[i]) * gamma * deltat + (nodes_[h]->getPastVelocity()[i]) * 1.0 +
// //                              (nodes_[h]->getPastAcceleration()[i]) * deltat * (1.0 - gamma);
// //                 }
// //                 nodes_[h]->setCurrentVelocity({vel[0], vel[1], vel[2]});
// //             }

// //             std::cout << "Iteration = " << interation
// //                       << "   x Norm = " << std::scientific << sqrt(normDeltaY / normInitialCoordinate)
// //                       << std::endl;

// //             if (sqrt(normDeltaY / normInitialCoordinate) <= tolerance)
// //                 break;
// //         }
// //         exportToParaview(timeStep);

// //         for (Node *node : nodes_)
// //         {
// //             std::vector<double> updatingCoordinate = node->getCurrentCoordinate();
// //             node->setPastCoordinate(updatingCoordinate);
// //             std::vector<double> updatingVel = node->getCurrentVelocity();
// //             node->setPastVelocity(updatingVel);
// //             std::vector<double> updatingAccel = node->getCurrentAcceleration();
// //             node->setPastAcceleration(updatingAccel);
// //         }

// //         file1 << nodes_[23]->getCurrentCoordinate()[0] - nodes_[23]->getInitialCoordinate()[0] << " " << time << std::endl;
// //     }
// // }



// // void Solid::exportToParaview(const int &loadstep)
// // {
// //     std::stringstream text;
// //     text << "output" << loadstep << ".vtu";
// //     std::ofstream file(text.str());

// //     //header
// //     file << "<?xml version=\"1.0\"?>"
// //          << "\n"
// //          << "<VTKFile type=\"UnstructuredGrid\">"
// //          << "\n"
// //          << "  <UnstructuredGrid>"
// //          << "\n"
// //          << "  <Piece NumberOfPoints=\"" << nodes_.size()
// //          << "\"  NumberOfCells=\"" << elements_.size()
// //          << "\">"
// //          << "\n";
// //     //nodal coordinates
// //     file << "    <Points>"
// //          << "\n"
// //          << "      <DataArray type=\"Float64\" "
// //          << "NumberOfComponents=\"3\" format=\"ascii\">"
// //          << "\n";
// //     for (Node *n : nodes_)
// //     {
// //         file << n->getCurrentCoordinate()(0) << " " << n->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
// //     }
// //     file << "      </DataArray>"
// //          << "\n"
// //          << "    </Points>"
// //          << "\n";
// //     //element connectivity
// //     file << "    <Cells>"
// //          << "\n"
// //          << "      <DataArray type=\"Int32\" "
// //          << "Name=\"connectivity\" format=\"ascii\">"
// //          << "\n";

// //     int n = (order_ + 1) * (order_ + 2) / 2.0;

// //     for (Element *e : elements_)
// //     {
// //         std::vector<Node *> conec;
// //         conec.reserve(n);
// //         conec = e->getConnection();

// //         for (int i = 0; i < n; i++)
// //         {
// //             file << conec[i]->getIndex() << " ";
// //         }
// //         file << "\n";
// //     }
// //     file << "      </DataArray>"
// //          << "\n";
// //     //offsets
// //     file << "      <DataArray type=\"Int32\""
// //          << " Name=\"offsets\" format=\"ascii\">"
// //          << "\n";
// //     int aux = 0;
// //     for (Element *e : elements_)
// //     {
// //         int n = e->getConnection().size();
// //         aux += n;
// //         file << aux << "\n";
// //     }
// //     file << "      </DataArray>"
// //          << "\n";
// //     //elements type
// //     file << "      <DataArray type=\"UInt8\" Name=\"types\" "
// //          << "format=\"ascii\">"
// //          << "\n";

// //     for (Element *e : elements_)
// //     {
// //         file << 69 << "\n";
// //     }
// //     file << "      </DataArray>"
// //          << "\n"
// //          << "    </Cells>"
// //          << "\n";
// //     //nodal results
// //     file << "    <PointData>"
// //          << "\n";

// //     file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
// //          << "Name=\"Displacement\" format=\"ascii\">"
// //          << "\n";
// //     for (Node *n : nodes_)
// //     {
// //         bounded_vector<double, 2> initial = n->getInitialCoordinate();
// //         bounded_vector<double, 2> current = n->getCurrentCoordinate();

// //         file << current(0) - initial(0) << " " << current(1) - initial(1) << "\n";
// //     }
// //     file << "      </DataArray> "
// //          << "\n";

// //     file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
// //          << "Name=\"Velocity\" format=\"ascii\">"
// //          << "\n";
// //     for (Node *n : nodes_)
// //     {
// //         bounded_vector<double, 2> currentVelocity = n->getCurrentVelocity();
// //         file << currentVelocity(0) << " " << currentVelocity(1) << "\n";
// //     }
// //     file << "      </DataArray> "
// //          << "\n";

// //     file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
// //          << "Name=\"Acceleration\" format=\"ascii\">"
// //          << "\n";
// //     for (Node *n : nodes_)
// //     {
// //         file << n->getCurrentAcceleration()(0) << " " << n->getCurrentAcceleration()(1) << "\n";
// //     }
// //     file << "      </DataArray> "
// //          << "\n";

// //     file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
// //          << "Name=\"CauchyNormalStress\" format=\"ascii\">"
// //          << "\n";
// //     for (Node *n : nodes_)
// //     {
// //         double cont = n->getStressState()(3);
// //         double aux1 = n->getStressState()(0);
// //         double aux2 = n->getStressState()(1);
// //         file << aux1 / cont << " " << aux2 / cont << "\n";
// //     }
// //     file << "      </DataArray> "
// //          << "\n";

// //     file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
// //          << "Name=\"CauchyShearStress\" format=\"ascii\">"
// //          << "\n";
// //     for (Node *n : nodes_)
// //     {
// //         double cont = n->getStressState()(3);
// //         double aux3 = n->getStressState()(2);
// //         file << aux3 / cont << "\n";
// //     }
// //     file << "      </DataArray> "
// //          << "\n";

// //     file << "    </PointData>"
// //          << "\n";
// //     //elemental results
// //     file << "    <CellData>"
// //          << "\n";

// //     file << "    </CellData>"
// //          << "\n";
// //     //footnote
// //     file << "  </Piece>"
// //          << "\n"
// //          << "  </UnstructuredGrid>"
// //          << "\n"
// //          << "</VTKFile>"
// //          << "\n";
// // }

//                 FiberNode *initial = fiber->getConnection()[0];
//                 FiberNode *end = fiber->getConnection()[1];

//                 Element *initialSolid = elements_[initial->getIndexSolidElement()];
//                 Element *endSolid = elements_[end->getIndexSolidElement()];

//                 bounded_vector<double, 2> dimensionlessInitial = initial->getDimensionlessCoordinates();
//                 bounded_vector<double, 2> dimensionlessEnd = end->getDimensionlessCoordinates();

//                 vector<double> phi_initial(n, 0.0);
//                 vector<double> phi_end(n, 0.0);
//                 phi_initial = domainShapeFunction(dimensionlessInitial(0), dimensionlessInitial(1));
//                 phi_end = domainShapeFunction(dimensionlessEnd(0), dimensionlessEnd(1));

//                 vector<double> correction_X1ElementInitial(n, 0.0);
//                 vector<double> correction_X2ElementInitial(n, 0.0);

//                 vector<double> correction_X1ElementEnd(n, 0.0);
//                 vector<double> correction_X2ElementEnd(n, 0.0);

//                 int auxiliar = 0;
//                 for (Node *solidNode : initialSolid->getConnection())
//                 {
//                     int index = solidNode->getIndex();
//                     correction_X1ElementInitial(auxiliar) = variationOfPI.first(2 * index);
//                     correction_X2ElementInitial(auxiliar) = variationOfPI.first(2 * index + 1);
//                     auxiliar = auxiliar + 1;
//                 }

//                 bounded_vector<double, 2> currentCoordinate;
//                 currentCoordinate=initial->getCurrentCoordinate();
//                 currentCoordinate(0)=currentCoordinate(0)+inner_prod(phi_initial, correction_X1ElementInitial);
//                 currentCoordinate(1)=currentCoordinate(1)+inner_prod(phi_initial,correction_X2ElementInitial);



// vector<int> indexNodes(2 * n, 0.0);
                            // for (int i = 0; i < n; i++)
                            // {
                            //     indexNodes(i) = elements_[indexSolidInitical]->getConnection()[i]->getIndex();
                            //     indexNodes(i + n) = elements_[indexSolidEnd]->getConnection()[i]->getIndex();
                            // }

                            // for (size_t i = 0; i < n; i++)
                            // {
                            //     if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
                            //     {
                            //         int dof = 2 * indexNodes(i);
                            //         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
                            //     }
                            //     if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
                            //     {
                            //         int dof = 2 * indexNodes(i) + 1;
                            //         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
                            //     }

                            //     if (fabs(elementMatrices.first(2 * (i + n))) >= 1.0e-15)
                            //     {
                            //         int dof = 2 * indexNodes(i + n);
                            //         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * (i + n)), ADD_VALUES);
                            //     }
                            //     if (fabs(elementMatrices.first(2 * (i + n) + 1)) >= 1.0e-15)
                            //     {
                            //         int dof = 2 * indexNodes(i + n) + 1;
                            //         ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * (i + n) + 1), ADD_VALUES);
                            //     }

                            //     for (size_t j = 0; j < n; j++)
                            //     {
                            //         //Primeiro quadrante
                            //         if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i);
                            //             int dof2 = 2 * indexNodes(j);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i);
                            //             int dof2 = 2 * indexNodes(j) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i) + 1;
                            //             int dof2 = 2 * indexNodes(j);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i) + 1;
                            //             int dof2 = 2 * indexNodes(j) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
                            //         }
                            //         //Segundo quadrante
                            //         if (fabs(elementMatrices.second(2 * i, 2 * (j + n)) >= 1.e-15))
                            //         {
                            //             int dof1 = 2 * indexNodes(i);
                            //             int dof2 = 2 * indexNodes(j + n);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * (j + n)), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * i, 2 * (j + n) + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i);
                            //             int dof2 = 2 * indexNodes(j + n) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * (j + n) + 1), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * i + 1, 2 * (j + n))) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i) + 1;
                            //             int dof2 = 2 * indexNodes(j + n);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * (j + n)), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * i + 1, 2 * (j + n) + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i) + 1;
                            //             int dof2 = 2 * indexNodes(j + n) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * (j + n) + 1), ADD_VALUES);
                            //         }
                            //         //Terceiro quadrante
                            //         if (fabs(elementMatrices.second(2 * (i + n), 2 * j)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n);
                            //             int dof2 = 2 * indexNodes(j);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * j), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * (i + n), 2 * j + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n);
                            //             int dof2 = 2 * indexNodes(j) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * j + 1), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * j)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n) + 1;
                            //             int dof2 = 2 * indexNodes(j);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * j), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * j + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n) + 1;
                            //             int dof2 = 2 * indexNodes(j) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * j + 1), ADD_VALUES);
                            //         }
                            //         //Quarto quadrante
                            //         if (fabs(elementMatrices.second(2 * (i + n), 2 * (j + n))) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n);
                            //             int dof2 = 2 * indexNodes(j + n);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * (j + n)), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * (i + n), 2 * (j + n) + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n);
                            //             int dof2 = 2 * indexNodes(j + n) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n), 2 * (j + n) + 1), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * (j + n))) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n) + 1;
                            //             int dof2 = 2 * indexNodes(j + n);
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * (j + n)), ADD_VALUES);
                            //         }
                            //         if (fabs(elementMatrices.second(2 * (i + n) + 1, 2 * (j + n) + 1)) >= 1.e-15)
                            //         {
                            //             int dof1 = 2 * indexNodes(i + n) + 1;
                            //             int dof2 = 2 * indexNodes(j + n) + 1;
                            //             ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * (i + n) + 1, 2 * (j + n) + 1), ADD_VALUES);
                            //         }
                            //     }
                            // }