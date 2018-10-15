#include "lattpoint.h"
#include <stdexcept>

LattPoint::LattPoint(int* initX, int initL) {
	L = initL;
	for(int i = 0; i < 4; i++)
		x[i] = initX[i] % L;
}

LattPoint::LattPoint(const LattPoint& x2) {
	L = x2.L;
	for(int i = 0; i < 4; i++)
		x[i] = x2.x[i];
}

LattPoint::LattPoint(const LattPoint&& x2) {
	L = x2.L;
	for(int i = 0; i < 4; i++)
		x[i] = x2.x[i];
}

void LattPoint::setX(const int* newX) {
	for(int i = 0; i < 4; i++)
		x[i] = newX[i] % L;
}

void LattPoint::getX(int* targetX) const {
	for(int i = 0; i < 4; i++)
		targetX[i] = x[i];
}

int LattPoint::getL() const {
	return L;
}

LattPoint LattPoint::scale(const int s) {
	int newX[4];
	for(int i = 0; i < 4; i++)
		newX[i] = (x[i] * s) % L;

	LattPoint scaled(newX, L);
	return scaled;
}

string LattPoint::toString() const {
	return "[" + std::to_string(x[0]) + "," + std::to_string(x[1]) + "," + std::to_string(x[2]) + "," + std::to_string(x[3]) + "]";
}

LattPoint LattPoint::operator=(const LattPoint& x2) {
	if (this != &x2) {
		L = x2.L;
		for(int i = 0; i < 4; i++)
			x[i] = x2.x[i];
	}
	return *this;
}

LattPoint LattPoint::operator=(const LattPoint&& x2) {
	if (this != &x2)
		for(int i = 0; i < 4; i++)
			x[i] = x2.x[i];
	return *this;
}

LattPoint LattPoint::operator+=(const LattPoint& x2) {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");

	int x2x[4];
	x2.getX(x2x);

	for(int i = 0; i < 4; i++)
		x[i] = (x[i] + x2x[i]) % L;
	return *this;
}

LattPoint LattPoint::operator+(const LattPoint& x2) const {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");
	
	int sumX[4];

	int x2x[4];
	x2.getX(x2x);

	for(int i = 0; i < 4; i++)
		sumX[i] = (x[i] + x2x[i]) % L;
	LattPoint sum(sumX, L);
	return sum;
}

LattPoint LattPoint::operator-=(const LattPoint& x2) {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");

	int x2x[4];
	x2.getX(x2x);

	for(int i = 0; i < 4; i++) {
		x[i] = x[i] - x2x[i];
		while(x[i] < 0)
			x[i] += L;
		x[i] %= L;
	}
	return *this;
}

LattPoint LattPoint::operator-(const LattPoint& x2) const {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");
	
	int diffX[4];

	int x2x[4];
	x2.getX(x2x);

	for(int i = 0; i < 4; i++) {
		diffX[i] = x[i] - x2x[i];
		while(diffX[i] < 0)
			diffX[i] += L;
		diffX[i] %= L;
	}
	LattPoint diff(diffX, L);
	return diff;
}

bool LattPoint::operator<(const LattPoint& x2) const {
	int x2x[4];
	x2.getX(x2x);

	for(int i = 0; i < 4; i++) {
		if(x[i] < x2x[i]) return true;
		else if(x[i] > x2x[i]) return false;
	}
	return false;
}