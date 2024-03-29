#pragma once
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

using namespace boost::numeric::ublas;

class Node
{
public:
	Node();

	Node(const int &index,
		 const bounded_vector<double, 2> &initialCoordinate); //(nº do nó, {x1, x2})

	~Node();

	int getIndex();

	bounded_vector<double, 2> getInitialCoordinate();

	bounded_vector<double, 2> getPastCoordinate();

	bounded_vector<double, 2> getPastVelocity();

	bounded_vector<double, 2> getPastAcceleration();

	bounded_vector<double, 2> getCurrentCoordinate();

	bounded_vector<double, 2> getCurrentVelocity();

	bounded_vector<double, 2> getCurrentAcceleration();

	bounded_vector<double, 4> getStressState();

	void setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate);

	void setPastVelocity(const bounded_vector<double, 2> &pastVelocity);

	void setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration);

	void setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate);

	void setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity);

	void setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration);

	void setStressState(const bounded_vector<double, 3> &stressState);

	void setZeroStressState();

	void incrementCurrentCoordinate(const int& direction, const double& value);


private:
	int index_;

	bounded_vector<double, 2> initialCoordinate_;

	bounded_vector<double, 2> pastCoordinate_;

	bounded_vector<double, 2> pastVelocity_;

	bounded_vector<double, 2> pastAcceleration_;

	bounded_vector<double, 2> currentCoordinate_;

	bounded_vector<double, 2> currentVelocity_;

	bounded_vector<double, 2> currentAcceleration_;

	bounded_vector<double, 4> stressState_; //{SigmaX1, SigmaX2, TalX1X2, contador}
};