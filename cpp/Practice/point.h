#ifndef POINT_H
#define POINT_H
#include <iostream>
using namespace std;

class Point {

private:
	// this variables/objects are only accesible by the methods in this class
	double x; // object
	double y; // object
	// no private methods are declared.

public:
	Point(); // first constructor
	Point(double xx, double yy); // second constructor
	double getX() const; 		// method to learn x coordiante.
	double getY() const; 		// method to learn y coordinate.
	void setX(double xx); 		// method to modify x coordinate.
	void setY(double yy); 		// method to modify y coordinate.
	double getR() const;  		// method to learn the radious.
	double getA() const;		// method to learn the angle.
	void setR(double r);		// method to modify the radious.
	void setA(double theta);	// method to modify the angle.
	void rotate(double theta);	// method/procedure intended to rotate point about an origin.
	bool operator==(const Point& Q) const; // operator declaration for c++
	bool operator!=(const Point& Q) const; // operator declaration for c++
	// becuase c++ doesn't know how to do Q==P if this are type Points!

};

// this methods/ prodedures cannot access the objects inside the class definition!
double dist(Point P, Point Q); 	// procedure to compute distance between 2 points.
Point midpoint(Point P, Point Q); // procedure to compute the midpoint between 2 points.
ostream& operator << ostream& os, const Point& P);

#endif