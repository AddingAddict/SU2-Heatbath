#include "link.h"

LinkAhhh::LinkAhhh() {
	a[0] = 1;
	for(int i = 1; i < 4; i++)
		a[i] = 0;
}

LinkAhhh::LinkAhhh(double* initA) {
	for(int i = 0; i < 4; i++)
		a[i] = initA[i];
}

LinkAhhh::LinkAhhh(const LinkAhhh& U2) {
	for(int i = 0; i < 4; i++)
		a[i] = U2.a[i];
}

LinkAhhh::LinkAhhh(const LinkAhhh&& U2) {
	for(int i = 0; i < 4; i++)
		a[i] = U2.a[i];
}

LinkAhhh LinkAhhh::inv() const {
	double norm2 = this->det();
	double invA[4] = {a[0]/norm2, -a[1]/norm2, -a[2]/norm2, -a[3]/norm2};

	LinkAhhh inv(invA);
	return inv;
}

double LinkAhhh::det() const {
	return a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
}

double LinkAhhh::tr() const {
	return 2*a[0];
}

LinkAhhh LinkAhhh::scale(const double s) {
	double newA[4];
	for(int i = 0; i < 4; i++)
		newA[i] = a[i] * s;

	LinkAhhh scaled(newA);
	return scaled;
}

void LinkAhhh::setA(const double* newA) {
	for(int i = 0; i < 4; i++)
		a[i] = newA[i];
}

void LinkAhhh::getA(double* targetA) const {
	for(int i = 0; i < 4; i++)
		targetA[i] = a[i];
}

LinkAhhh LinkAhhh::operator=(const LinkAhhh& U2) {
	if (this != &U2)
		for(int i = 0; i < 4; i++)
			a[i] = U2.a[i];
	return *this;
}

LinkAhhh LinkAhhh::operator=(const LinkAhhh&& U2) {
	if (this != &U2)
		for(int i = 0; i < 4; i++)
			a[i] = U2.a[i];
	return *this;
}

LinkAhhh LinkAhhh::operator+=(const LinkAhhh& U2) {
	double U2a[4];
	U2.getA(U2a);

	for(int i = 0; i < 4; i++)
		a[i] += U2a[i];
	return *this;
}

LinkAhhh LinkAhhh::operator+(const LinkAhhh& U2) const{
	double sumA[4];

	double U2a[4];
	U2.getA(U2a);

	for(int i = 0; i < 4; i++)
		sumA[i] = a[i] + U2a[i];
	LinkAhhh sum(sumA);
	return sum;
}

LinkAhhh LinkAhhh::operator-=(const LinkAhhh& U2) {
	double U2a[4];
	U2.getA(U2a);

	for(int i = 0; i < 4; i++)
		a[i] -= U2a[i];
	return *this;
}

LinkAhhh LinkAhhh::operator-(const LinkAhhh& U2) const{
	double diffA[4];

	double U2a[4];
	U2.getA(U2a);

	for(int i = 0; i < 4; i++)
		diffA[i] = a[i] - U2a[i];
	LinkAhhh diff(diffA);
	return diff;
}

LinkAhhh LinkAhhh::operator*=(const LinkAhhh& U2) {
	double a1 = a[0];
	double b1 = a[3];
	double c1 = a[2];
	double d1 = a[1];

	double U2a[4];
	U2.getA(U2a);

	double a2 = U2a[0];
	double b2 = U2a[3];
	double c2 = U2a[2];
	double d2 = U2a[1];

	double prodA[4];

	prodA[0] = a1*a2 - b1*b2 - c1*c2 - d1*d2;
	prodA[3] = a1*b2 + b1*a2 + c1*d2 - d1*c2;
	prodA[2] = a1*c2 - b1*d2 + c1*a2 + d1*b2;
	prodA[1] = a1*d2 + b1*c2 - c1*b2 + d1*a2;

	this->setA(prodA);
	return *this;
}

LinkAhhh LinkAhhh::operator*(const LinkAhhh& U2) const{
	double prodA[4];

	double a1 = a[0];
	double b1 = a[3];
	double c1 = a[2];
	double d1 = a[1];

	double U2a[4];
	U2.getA(U2a);

	double a2 = U2a[0];
	double b2 = U2a[3];
	double c2 = U2a[2];
	double d2 = U2a[1];

	prodA[0] = a1*a2 - b1*b2 - c1*c2 - d1*d2;
	prodA[3] = a1*b2 + b1*a2 + c1*d2 - d1*c2;
	prodA[2] = a1*c2 - b1*d2 + c1*a2 + d1*b2;
	prodA[1] = a1*d2 + b1*c2 - c1*b2 + d1*a2;

	LinkAhhh prod(prodA);
	return prod;
}