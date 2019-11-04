#pragma once
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
//#include "../Element.h"

using namespace boost::numeric::ublas;

class FiberNode
{
public:
	FiberNode();

	FiberNode(const int &index,
			  const bounded_vector<double, 2> &initialCoordinate); //(nº do nó, {x1, x2})

	~FiberNode();

	int getIndex();

	bounded_vector<double, 2> getInitialCoordinate();

	bounded_vector<double, 2> getPastCoordinate();

	bounded_vector<double, 2> getPastVelocity();

	bounded_vector<double, 2> getPastAcceleration();

	bounded_vector<double, 2> getCurrentCoordinate();

	bounded_vector<double, 2> getCurrentVelocity();

	bounded_vector<double, 2> getCurrentAcceleration();

	bounded_vector<double, 2> getDimensionlessCoordinates();

	int getIndexSolidElement();

	double getNormalForce();

	void setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate);

	void setPastVelocity(const bounded_vector<double, 2> &pastVelocity);

	void setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration);

	void setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate);

	void setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity);

	void setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration);

	void setIncidenceInSolid(const int &indexElement, const bounded_vector<double, 2> &incidence);

	void setNormalForce(const double &force);


private:
	int index_;

	bounded_vector<double, 2> initialCoordinate_;

	bounded_vector<double, 2> pastCoordinate_;

	bounded_vector<double, 2> pastVelocity_;

	bounded_vector<double, 2> pastAcceleration_;

	bounded_vector<double, 2> currentCoordinate_;

	bounded_vector<double, 2> currentVelocity_;

	bounded_vector<double, 2> currentAcceleration_;

	double normalForce_;

	int indexElement_;

	bounded_vector<double, 2> incidence_; // xsi1, xsi2
};