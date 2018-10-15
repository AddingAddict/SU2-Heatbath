class Link {
	private:
		double a[4];					// 4-vector representation
	public:
		Link();							// constructor
		Link(double* initA);
		Link(const Link& U2);			// copy constructor
		Link(const Link&& U2);			// move constructor
		Link inv() const;				// multiplicative inverse
		double det() const;				// determinant of representation
		double tr() const;				// trace of representation
		Link scale(const double s);		// returns scaled version of link variable 
		void setA(const double* newA);	// setter function
		void getA(double* targetA) const;	// getter function
		Link operator=(const Link& U2);
		Link operator=(const Link&& U2);
		Link operator+=(const Link& U2);
		Link operator+(const Link& U2) const;
		Link operator-=(const Link& U2);
		Link operator-(const Link& U2) const;
		Link operator*=(const Link& U2);
		Link operator*(const Link& U2) const;
};