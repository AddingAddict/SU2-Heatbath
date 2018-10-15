#include <string>
using namespace std;

class LattPoint {
	private:
		int x[4];						// 4-vector representation
		int L;							// lattice size
	public:
		LattPoint(int* initX, int initL);// constructor
		LattPoint(const LattPoint& x);	// copy constructor
		LattPoint(const LattPoint&& x);	// move constructor
		void setX(const int* newX);		// setter function
		void getX(int* targetX) const;	// getter function
		int getL() const;
		LattPoint scale(const int s);	// scale by integer
		string toString() const;
		LattPoint operator=(const LattPoint& U2);
		LattPoint operator=(const LattPoint&& U2);
		LattPoint operator+=(const LattPoint& U2);
		LattPoint operator+(const LattPoint& U2) const;
		LattPoint operator-=(const LattPoint& U2);
		LattPoint operator-(const LattPoint& U2) const;
		bool operator<(const LattPoint& U2) const;
};