//Emily Maxey
//Hydrogen Wave Function Generator
//4/29/2021

/*	
	This file produces the normalized wave functions for hydrogen.
	You will likely not need to interfere in this script at all.  Refer to takeIntegral.h to create
	new perturbations and order them to be used in the main function in generateWaveFunctions.C.
	
	It follows the Generating functions defined as:
	
	Radial:
	R(r) = -[(2/na) * ((n-l-1)!/2n((n+l)!))]^.5 * e^(na) L(2/na)
		where L is the Associated Laguerre polynomials.  We use p= 2l+1, q = n+l and x = 2/(na)
		L(x) = (d/dx)^pL(x)	with the Laguerre polynomials: L(x)= e^x (d/dx)^q(x^qe^-x)
		
	Spherical:
	Y(𝜭,𝞍) = (-1)^m [((2l+1)(l-m)!/(4PI(l+m)!]^.5 e^im𝞍 * P(cos𝜭)
		Where P is the associated Legendre polynomial,
		P(x) = (-1)^m sqrt((1-x^2)^m) * (d/dx)^m * P(x) and Legendre polynomials P(x) = (-1)^l/(2^l*l!) * (d/dx)^l(1-x^2)^l

	To actually generate these solutions, we defined two primary classes,
	Terms and Expressions.  A term is comprised of a coefficient, a base and an exponent.
	Since some exponents require a variable, the exponent is broken down into an integer
	and a string component.  All necessary get and set functions were also created.
	
	Expressions consist of several specific terms as well as some additional coefficients.
	Since there are some constants inside of square roots as we generate the wave functions,
	Expressions have specific variables for the integer and string components in both the 
	numerator and denominator.  Once again any and all necessary get and set functions were 
	defined inside the class.
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

/*	class declarations */

class Term
{
	private:
		int coeff = 1;
		string baseStr = "";
		double power = 1;
		string powerStr = "";
		
	public:
		double getPower() { return power; };
		int getCoeff() { return coeff; };
		string getBaseStr() { return baseStr; };
		string getPowerStr() { return powerStr; };
		
		void setPower(double newPower) { power = newPower; };
		void setCoeff(int newCoeff) { coeff = newCoeff; };
		void setBaseStr(string newBase) { baseStr = newBase; };
		void setTerm(Term myTerm);
		void setPowerStr(string newPower) { powerStr = newPower; };
		
		void phiSphericalTerm(int m);
};

class Expression
{
	private:
		int numerator = 1;
		int numInSquare = 1;
		int denominator = 1;
		int denomInSquare = 1;

		string numeratorStr = "";
		string denominatorStr = "";
		string numInSquareStr = "";
		string denomInSquareStr = "";
		
		Term termNumerator;		//this is the cos term for spherical, r term for radial
		Term termNum2;			//this is the sin term
		Term expTerm;			//this is the e^() term
		Term denomTerm;			//This is to take care of a in the denominator
		Term sphTerm;		//This is to take care of the e^(imΦ) term.
		Term soidalPhiTerm;
		
	public:
	
		Expression();

		int getNumerator() { return numerator; };
		int getDenominator() { return denominator; };
		int getSquareNum() { return numInSquare; };
		int getDenomInSquare() { return denomInSquare; };
		
		string getNumeratorStr() { return numeratorStr; };
		string getDenominatorStr() { return denominatorStr; };
		string getSquareNumStr() { return numInSquareStr; };
		string getDenomInSquareStr() { return denomInSquareStr; };
		
		Term getTermNumerator() { return termNumerator; };
		Term getTerm2() { return termNum2; };
		Term getExpTerm() { return expTerm; };
		Term getDenomTerm() { return denomTerm; };
		Term getSphTerm() { return sphTerm; };
		Term getSoidalPhiTerm() {return soidalPhiTerm;}
		
		void setNumerator(int newNum) { numerator = newNum; };
		void setDenominator(int newDenom) { denominator = newDenom; };
		void setSquareNum(int newRootNum) { numInSquare = newRootNum; };
		void setDenomInSquare(int newRootDenom) { denomInSquare = newRootDenom; };
		void setNumeratorStr(string newNum) { numeratorStr = newNum; };
		void setDenominatorStr(string newDenom) { denominatorStr = newDenom; };
		void setSquareNumStr(string newRootNum) { numInSquareStr = newRootNum; };
		void setDenomInSquareStr(string newRootDenom) { denomInSquareStr = newRootDenom; };
		void setTermNumerator(Term newTerm);
		void setDenomTerm(Term newTerm);
		void termToExpression(Term myTerm);
		void setTermNum2(Term newTerm2);
		void term2ToExpression(Term myTerm2);
		void setExpTerm(Term myExpTerm);
		void setSoidalPhiTerm(Term newPhiTerm);
		void multiplyExpressions(Expression multiplier);
		void simplifyRadical();
		void simplifyFraction();
		void takeSquareRoot();
		void expressionCube();
		void setSphTerm(Term newTerm);
};

class Radial
{
	private:
		Expression primary;
		vector<Expression> laguerrePolynomials;
		
	public:
		Expression getPrimary() { return primary; };
		vector<Expression> getLaguerrePolynomials() { return laguerrePolynomials; };
		
		void setPrimary(Expression primary);
		void setLaguerrePolynomials(vector<Expression> laguerreSolutions);
};

/*	Function Prototypes	*/
//Function explanations are provided with each individual definition.

vector<Expression> legendrePolynomial(int l);

vector<Expression> associatedLegendrePolynomial(int l, int m);

Expression cubedRadialPart(int n);

vector<Term> polynomialExpansion(int l);

void takeDerivative(vector<Term> &terms);

Term takeDerivative(Term);

Expression insideSphericalBrackets(int l, int m);

vector<Expression> LaguerreSolutions(int n, int l);

int nChooseK(int n, int k);

int factorial(int n);

int sphericalSign(int m);

vector<Expression> totalSpherical(int sign, Expression insideSphericalBrackets, vector<Expression> associatedLegendres, Term sphTerm);

vector<Expression> totalSpherical(int n, int l, int m);

int gcd(int u, int v);

Expression getRadialSqrRootPart(int n, int l);

vector<Expression> sphericalPortion(int l, int m);

void setExpression1To2(Expression& exp1, Expression exp2);

ostream& outputExpression(ostream& out, Expression expression);

ostream& denominatorOutput(ostream& out, Expression expression);

ostream& numeratorOutput(ostream& out, Expression expression);

ostream& polynomialOutput(ostream& out, vector<Term> terms);

ostream& termOutput(ostream& out, Term term);

Expression powerToTheLRadial(int l, int n);

Term exponentialRadialTerm(int n);

bool termIsEmpty(Term term);

vector<Expression> radialTotal(Expression inSqrRoot, Term exponentialTerm, Expression powerToL, vector<Expression> laguerreSolutions);

vector<Expression> createRadial(int n, int l);

/*
	Basically sets up the process for the other totalSpherical overloaded function
	to actually create the spherical.  Creates the initial expressions and signs that
	will be needed from only the quantum numbers.
*/
vector<Expression> totalSpherical(int n, int l, int m)
{
	Term sphTerm;
	int sphSign = sphericalSign(m);

	Expression insideBrackets;
	setExpression1To2(insideBrackets, insideSphericalBrackets(l,m));
				
	vector<Expression> assLegendres = associatedLegendrePolynomial(l,m);

				
	sphTerm.phiSphericalTerm(m); 
	
	return totalSpherical(sphSign, insideBrackets, assLegendres, sphTerm);
}

/*
	Call to output the wave function, only takes a single expression.  If you wish
	to output a wave function with multiple expression terms, you need to do a for
	loop looping through each index of the vector<expressions> calling this function
	for every index.  Operates by outputting the numerator and denominator seperately.
*/
ostream& outputExpression(ostream& out, Expression expression)
{
	if(expression.getNumerator()==0)			//If the numerator is 0, no point in outputting the rest of the wave function.
		return out;

	numeratorOutput(out, expression);			
	out << endl << "----------------" << endl;
	denominatorOutput(out, expression);
	out << endl;
	
	return out;
}

/*	
	Outputs the denominator of our expression	
*/
ostream& denominatorOutput(ostream& out, Expression expression)
{
	out << expression.getDenominator() << expression.getDenominatorStr();

	termOutput(cout , expression.getDenomTerm());
	if(expression.getDenomInSquare() != 1)
		out << "sqrt(" << expression.getDenomInSquare() << expression.getDenomInSquareStr() << ")";
	else if(expression.getDenomInSquareStr() != "")
		out << "sqrt(" << expression.getDenomInSquareStr() << ")";
	return out;
}

/*	
	Used during development to make sure factoring worked as expected
	for the Legendre Polynomials, is not obsolete.						
*/
ostream& polynomialOutput(ostream& out, vector<Term> terms)
{
	for(int i = 0 ; i < terms.size() ; i++)
	{
		termOutput(out, terms[i]);
		out << " + " ;
	}
	
	return out;
}

/*	
	Outputs an indiividual term with its base and power
	Note that it still might need to output the powerStr
	even if the power is just 1(our case, so that it 
	does print out the iPhi).							
*/
ostream& termOutput(ostream& out, Term term)
{
	if(term.getCoeff() != 1)
	{
		out << term.getCoeff();
	}
	if(term.getPower() == 1)
	{
		out << term.getBaseStr();
		if(term.getPowerStr() != "")
		{
			out << "^(" << term.getPowerStr() << ")";
		}
	}
	else if(term.getPower() != 0)
	{
		out << term.getBaseStr() << "^(" << term.getPower() << term.getPowerStr() << ")";
	}
	return out;
}

/*	
	Numerator output for expression needs to include the numerator, numInSquare
	term1 which holds our cos, term2 which holds sign and expTerm which holds
	the e^imphi.  Uses termOutput function to do this.  Look at those rules
	for more specifics.														
*/
ostream& numeratorOutput(ostream& out, Expression expression)
{
	out << expression.getNumerator() << expression.getSquareNumStr() ;
	if(expression.getSquareNum() != 1)
		out << "sqrt(" << expression.getSquareNum() << expression.getSquareNumStr() << ")";
	else if(expression.getSquareNumStr() != "")
		out << "sqrt(" << expression.getSquareNumStr() << ")";
		
	termOutput(out, expression.getTermNumerator());
	termOutput(out, expression.getTerm2());
	termOutput(out, expression.getExpTerm());
	termOutput(out, expression.getSphTerm());
	termOutput(out, expression.getSoidalPhiTerm());
	
	return out;
}

/*	
	Computes the factorial of the given integer: Necessary to calculate the Spherical portion.
	Returns n!									
*/
int factorial(int n)
{
	if(n == 0)
		return 1;
	
	if(n == 1)
		return 1;
		
		
	int fact = 1;
	for(int i = n ; i >1 ; i--)
	{
		fact *= i;
	}
	
	return fact;
}

/*	returns (-1)^m		*/
int sphericalSign(int m)
{
	return pow(-1, m);
}

/*	This is used for factoring purposes.  The formula comes from 
	pascals triangle.  It is the coeff of that specific term in 
	the expanded polynomial					
	Returns 	n!
			---------
			 k!(n-k)!					*/
int nChooseK(int n, int k)
{
	int numFact = factorial(n);
	int denomFact1  = factorial(k);
	int denomFact2 = factorial(n-k);
	
	int coeff = numFact / denomFact1;
	coeff /= denomFact2;
	
	return coeff;
}

/*	Gives us the greatest Common Denominator.
	Reused code from past project 1337.57				*/
int gcd(int u, int v)
{
	if(u == 0 && v ==0)
		return 0;
	if(v == 0)
		return abs(u);
	if(u == 0)
		return abs(v);
	
	int divisor = min(abs(u), abs(v));

	
	while((u%divisor)!=0 || (v%divisor)!=0 )
	{
		divisor--;
	}
	
	
	return divisor;
}

/*Returns true if the term contains no information, false otherwise	*/
bool termIsEmpty(Term term)
{
	if(term.getCoeff() != 1 && term.getBaseStr() != "" && term.getPower() != 1 && term.getPowerStr() != "")
		return true;
	return false;
}

/*	Constructor.  The exponenet power should start at 0
	if we have an e^x factor, the power will be changed 
	later, we dont want one to start.					*/
Expression :: Expression()
{
	expTerm.setPower(0);
	denomTerm.setPower(0);
	
}

/*	
	Used to take care of the 1/2 power of the square brackets
	in calculating the spherical portion of our wave function.
	Essentially sets the normal number to being inside the sqrt.						
*/
void Expression :: takeSquareRoot()
{
	if(getDenomInSquare() != 1 || getDenomInSquareStr() != "")
		cout << endl << "---------------BIG PROBLEM, TAKING SQRT OF SOMETHING IN SQRT !!!!!-----------------------\n Please reference takeSquareRoot in generateWaveFunctions.h" << endl;	
		//This is a precautionary statement.  If you encounter this, you will need to create your own function and if cases to handle what is essentailly a fourth root.
		//It has not occurred yet and as such is not accounted for currently.  
		
	setSquareNum(getNumerator());
	setSquareNumStr(getNumeratorStr());
	setDenomInSquare(getDenominator());
	setDenomInSquareStr(getDenominatorStr());
	
	setNumerator(1);
	setNumeratorStr("");
	setDenominator(1);
	setDenominatorStr("");
	
	Term newTerm;
	newTerm = getDenomTerm();
	newTerm.setPower(newTerm.getPower() / 2);
	setDenomTerm(newTerm);
	
	getDenomTerm().setPower(getDenomTerm().getPower() / 2);
	getTermNumerator().setPower(getTermNumerator().getPower() / 2);
	getTerm2().setPower(getTerm2().getPower() / 2);
	getExpTerm().setPower(getExpTerm().getPower() / 2);
}

/*	
	Pretty self explanatory, cubes the expression which calls it
*/
void Expression :: expressionCube()
{
	setNumerator(pow(getNumerator(), 3));
	setDenominator(pow(getDenominator() , 3));
	setNumerator(getNumerator() * getSquareNum());
	if(getSquareNumStr() != "")
		setNumeratorStr(getNumeratorStr() + " * " + getSquareNumStr());
	setDenominator(getDenominator() * getDenomInSquare());
	if(getDenomInSquareStr() != "")
		setDenominatorStr(getDenominatorStr() + " * " + getDenomInSquareStr());
	
	Term newTerm;
	newTerm.setTerm(getDenomTerm());
	
	newTerm.setPower(3);
	
	setDenomTerm(newTerm);
}

/*	Set function for the denom term which carries the a term in the denominator.	*/
void Expression :: setDenomTerm(Term newTerm)
{
	denomTerm.setPower(newTerm.getPower());
	setDenominator(getDenominator() * newTerm.getCoeff());
	denomTerm.setCoeff(1);
	denomTerm.setBaseStr(newTerm.getBaseStr());
}

/*	Set function for the soidalPhiTerm	*/
void Expression :: setSoidalPhiTerm(Term newTerm)
{
	soidalPhiTerm.setPower(newTerm.getPower());
	setNumerator(getNumerator() * newTerm.getCoeff());
	soidalPhiTerm.setCoeff(1);
	soidalPhiTerm.setBaseStr(newTerm.getBaseStr());
}

/*	Set function for the sphTerm which holds the e^(imphi) term		*/
void Expression :: setSphTerm(Term newTerm)
{
	sphTerm.setPower(newTerm.getPower());
	setNumerator(getNumerator() * newTerm.getCoeff());
	sphTerm.setCoeff(1);
	sphTerm.setBaseStr(newTerm.getBaseStr());
	sphTerm.setPowerStr(newTerm.getPowerStr());
}

/*	Name explains, basically finds the greatest perfect
	square that it can take out of the radical. Does the 
	same process for both the numerator and denom. For ex:
	sqrt(24).  Greatest is 4, so it becomes 2sqrt(6)	*/
void Expression :: simplifyRadical()
{
	int inSquareNum = getSquareNum();
	int max = pow(inSquareNum, .5);
	int iSquared = 1;
	
	for(int i = max; i > 0 ; i--)
	{
		
		iSquared = pow(i,2);
		if(inSquareNum % iSquared == 0)
		{
			setNumerator(getNumerator() * i);
			setSquareNum(getSquareNum() / pow(i,2));
			inSquareNum = getSquareNum();
		}
	}
	
	
	int inSquareDenom = getDenomInSquare();
	max = pow(inSquareDenom, .5);
	iSquared = 1;
	
	for(int i = max; i > 0 ; i--)
	{
		
		iSquared = pow(i,2);
		if(inSquareDenom % iSquared == 0)
		{
			setDenominator(getDenominator() * i);
			setDenomInSquare(getDenomInSquare() / pow(i,2));
			inSquareDenom = getDenomInSquare();
		}
	}
}

//Initializes termNum2 which holds the sin term
void Expression :: setTermNum2(Term newTerm2)
{
	termNum2.setPower(newTerm2.getPower());
	termNum2.setCoeff(newTerm2.getCoeff());
	termNum2.setBaseStr(newTerm2.getBaseStr());
}

/*	Initializes expTerm.  Note that since our e^(stuff)
	term requires a power string, this sets that part 
	as well.										*/
void Expression :: setExpTerm(Term myExpTerm)
{
	expTerm.setPower(myExpTerm.getPower());
	expTerm.setCoeff(myExpTerm.getCoeff());
	expTerm.setBaseStr(myExpTerm.getBaseStr());
	expTerm.setPowerStr(myExpTerm.getPowerStr());
}

/*	first simplifies the radicals(uses the previously defined
	function and then divides num and denom by the gcd.		*/
void Expression :: simplifyFraction()
{
	simplifyRadical();
	
	int myGCD = gcd(getNumerator(), getDenominator());
	
	setNumerator( getNumerator() / myGCD );
	setDenominator( getDenominator() / myGCD );
	
	myGCD = gcd(getDenomInSquare(), getSquareNum());
	
	setDenomInSquare( getDenomInSquare() / myGCD );
	setSquareNum( getSquareNum() / myGCD );
}

/*	Initializes termNumerator which holds the cos term for spherical, r for radial				*/
void Expression :: setTermNumerator(Term newTerm)
{
	termNumerator.setPower(newTerm.getPower());
	setNumerator(getNumerator() * newTerm.getCoeff());
	termNumerator.setCoeff(1);
	termNumerator.setBaseStr(newTerm.getBaseStr());
}

/*	Given a term that we need to put into an expression as termNumerator
	This takes that term and puts it into the expression
	properly.											*/
void Expression :: termToExpression(Term myTerm)
{
	setNumerator(getNumerator() * myTerm.getCoeff());
	myTerm.setCoeff(1);
	setTermNumerator(myTerm);
}

/*	Same idea as above but for term2 which is sin	*/
void Expression :: term2ToExpression(Term myTerm2)
{
	setNumerator(getNumerator() * myTerm2.getCoeff());
	myTerm2.setCoeff(1);
	setTermNum2(myTerm2);
}

/*	multiply two expressions, necessary to get the speherical together
	since we multiply the bracket stuff and legendre polynomials together.	*/
void Expression :: multiplyExpressions(Expression multiplier)
{
	setNumerator(multiplier.getNumerator() * getNumerator());
	setSquareNum(multiplier.getSquareNum() * getSquareNum());
	setDenominator(multiplier.getDenominator() * getDenominator());
	setDenominatorStr(multiplier.getDenominatorStr() + getDenominatorStr());
	setDenomInSquare(multiplier.getDenomInSquare() * getDenomInSquare());
	setDenomInSquareStr(multiplier.getDenomInSquareStr() + getDenomInSquareStr());
	
	setTermNumerator(multiplier.getTermNumerator() );
	setTermNum2(multiplier.getTerm2());
	setSphTerm(multiplier.getSphTerm());
		
	setDenomTerm(multiplier.getDenomTerm());
	
	setExpTerm(multiplier.getExpTerm());
}

/*	This will output(2/na)^3 as part of the radial
	wave function that is inside the square root, this 
	is where it will start.								*/
Expression cubedRadialPart(int n)
{
	Expression beforeCube;
	beforeCube.setNumerator(2);
	
	Term denominatorTerm;
	denominatorTerm.setBaseStr("a");
	denominatorTerm.setPower(1);

	beforeCube.setDenomTerm(denominatorTerm);
	beforeCube.setDenominator(n);

	beforeCube.expressionCube();

	return beforeCube;
}

/*	Calculates the inside of the square Brackets for spherical	*/
Expression insideSphericalBrackets(int l, int m)
{
	Expression sphereBrackets;
	//((2l+ 1)(l-m)!) / (4pi(l + m)!)
	int numFact = factorial(l-m);
	sphereBrackets.setNumerator((2*l+1) * numFact);
	
	int denomFact = factorial(l+m);
	sphereBrackets.setDenominator(4*denomFact);
	sphereBrackets.setDenominatorStr("pi");
	
	sphereBrackets.takeSquareRoot();
	
	return sphereBrackets;
}

/*	Actually calculates the legendrePolynomial for the specific l case
	Needed to calculate the associted legendre polynomials.	Have to take 
	multiple derivatives to get this value, after getting the (x^2-1)^l
	by using the polynomial expansion function							*/
vector<Expression> legendrePolynomial(int l)
{
	// P(cos) = (1/(2^l l!))* (d/dx)^l (x^2 - 1)^l
	vector <Expression> legendreSolution;
	
	vector<Term> derivatives = polynomialExpansion(l);
	
	for(int i = 0; i < l ; i++)
	{
		takeDerivative(derivatives);
	}
	
	legendreSolution.resize(derivatives.size());
	
	
	int denomFact = factorial(l);
	for(int i = 0; i < derivatives.size() ; i++)
	{
		legendreSolution[i].termToExpression(derivatives[i]);
		legendreSolution[i].setDenominator(denomFact * (pow(2,l)));
	}
	
	return legendreSolution;
}

/* Returns Plm FINALLY, still uses derivatives and has to take normal
	legendre polynomials first.  Also goes aheadd and sets term 2 which 
	is the sin information and takes care of the sign.					*/
vector<Expression> associatedLegendrePolynomial(int l, int m)
{
	vector<Expression> legendrePolynomials = legendrePolynomial(l);
	
	vector<Term> terms;
	terms.resize(legendrePolynomials.size());
	for(int i = 0 ; i < terms.size() ; i++)
	{
		terms[i].setTerm(legendrePolynomials[i].getTermNumerator());
	}

	for(int i = 0; i < m ; i++)
	{
		takeDerivative(terms);
	}
	
	int sign = sphericalSign(m);
	Term sinStuff;
	sinStuff.setBaseStr("sin𝜭");
	sinStuff.setPower(m);
	for(int i = 0; i < legendrePolynomials.size() ; i++)
	{
		legendrePolynomials[i].termToExpression(terms[i]);
		legendrePolynomials[i].setNumerator(legendrePolynomials[i].getNumerator() * sign);
		legendrePolynomials[i].term2ToExpression(sinStuff);
	}
	
	return legendrePolynomials;
}

/*	Puts together the rest of the spherical wave function
	mainly multiplies the bracket and associatedLegendrePolynomial
	as well as puts in the expTerm.									*/
vector<Expression> sphericalPortion(int l, int m)
{
	vector<Expression> spherical = associatedLegendrePolynomial(l,m);
	
	Expression squareBracket = insideSphericalBrackets(l,m);
	squareBracket.takeSquareRoot();
	int sign = sphericalSign(m);
	
	Term exp;
	exp.setPower(m);
	exp.setPowerStr("i𝞍");
	exp.setBaseStr("e");
	squareBracket.setNumerator(squareBracket.getNumerator() * sign);
	for(int i = 0; i < spherical.size() ; i++)
	{
		spherical[i].multiplyExpressions(squareBracket);
		spherical[i].setExpTerm(exp);
	}

	return spherical;
}

/*	Tidies up a list of terms.  Useful for after derivates
	If the numerator is 0, then the term is useless, lets 
	delete it out of this vector.							*/
void cleanUpTerms(vector<Term> &terms)
{
	int i = 0;
	while(i < terms.size())
	{
		if(terms[i].getCoeff() == 0)
		{
			terms.erase(terms.begin()+i);
		}
		else
			i++;
	}
}

/*	calls derivative function for each term in a 
	vector of terms. then calls clean up function.	*/
void takeDerivative(vector<Term> &terms)
{
	for(int i = 0 ; i < terms.size() ;i++)
	{
		terms[i] = takeDerivative(terms[i]);
	}
	
	cleanUpTerms(terms);
}

/*	
	Expands the polynomial using pascals triangle.
*/
vector<Term> polynomialExpansion(int l)
{
	//(x^2 - 1)^l
	vector <Term> expandedPolyNomial;
	
	Term currentTerm;
	currentTerm.setBaseStr("cos𝜭");
	
	int coeff = 1;

	for(int k = l ; k >= 0 ; k--)
	{
		coeff = nChooseK(l, k);
		
		coeff *= sphericalSign(l-k);
			
		currentTerm.setCoeff(coeff);
		currentTerm.setPower(2* k);
		
		expandedPolyNomial.push_back(currentTerm);
	}
	return expandedPolyNomial;
}

/*	takes the derivative of an individual term
	Uses the formula d/dx (ax^n) = (n*a)x^(n-1)	
	Since terms are inherently defined this way,
	we do not need to account for more complicated 
	types of derivatives
*/
Term takeDerivative(Term term)
{
	Term afterDerivative;
	if(term.getPower() == 0)
	{
		afterDerivative.setCoeff(0);
		return afterDerivative;
	}
	
	afterDerivative.setCoeff(term.getCoeff() * term.getPower());
	afterDerivative.setBaseStr(term.getBaseStr());
	afterDerivative.setPower(term.getPower() - 1);
	
	return afterDerivative;
}

/*	Declares a term and sets up the power and base correctly	*/
void Term :: setTerm(Term myTerm)
{
	setPowerStr(myTerm.getPowerStr());
	setPower(myTerm.getPower());
	setCoeff(myTerm.getCoeff());
	setBaseStr(myTerm.getBaseStr());
}

/*	sets the main part of the radial object to the things
	necessary for the expression.						*/
void Radial :: setPrimary(Expression exp2)
{
	primary.setNumerator(exp2.getNumerator());
	primary.setDenominator(exp2.getDenominator());
	primary.setSquareNum(exp2.getSquareNum());
	primary.setDenomInSquare(exp2.getDenomInSquare());
	primary.setNumeratorStr(exp2.getNumeratorStr());
	primary.setDenominatorStr(exp2.getDenominatorStr());
	primary.setSquareNumStr(exp2.getSquareNumStr());
	primary.setDenomInSquareStr(exp2.getDenomInSquareStr());
	primary.setTermNumerator(exp2.getTermNumerator());
	primary.setDenomTerm(exp2.getDenomTerm());
	primary.setTermNum2(exp2.getTerm2());
	primary.setExpTerm(exp2.getExpTerm());
}

/*
	Calculates Everything within the square root of the radial portion of the
	wave function. Reference genreating funcitons at the beginning for more
	information if necessary
*/
Expression getRadialSqrRootPart(int n, int l)
{
	Expression cubedPart;
	cubedPart = cubedRadialPart(n);
	cubedPart.setNumerator(cubedPart.getNumerator() * factorial(n-l-1));
	int denomMult = factorial(n+l);
	denomMult = pow(denomMult, 3);
	denomMult *= (2*n);
	
	cubedPart.setDenominator(cubedPart.getDenominator() * denomMult);
	cubedPart.takeSquareRoot();
	return cubedPart;
}

/*	Sets the laguerre polynomial part of a radial object	*/
void Radial :: setLaguerrePolynomials(vector<Expression> laguerreSolutions)
{
	laguerrePolynomials.resize(laguerreSolutions.size());
	for(int i = 0 ; i < laguerrePolynomials.size() ; i++)
	{
		setExpression1To2(laguerrePolynomials[i], laguerreSolutions[i]);
	}
}

/*	Sets expression1 to be equivalent to exp2	*/
void setExpression1To2(Expression& exp1, Expression exp2)
{
	exp1.setNumerator(exp2.getNumerator());
	exp1.setDenominator(exp2.getDenominator());
	exp1.setSquareNum(exp2.getSquareNum());
	exp1.setDenomInSquare(exp2.getDenomInSquare());
	exp1.setNumeratorStr(exp2.getNumeratorStr());
	exp1.setDenominatorStr(exp2.getDenominatorStr());
	exp1.setSquareNumStr(exp2.getSquareNumStr());
	exp1.setDenomInSquareStr(exp2.getDenomInSquareStr());
	exp1.setTermNumerator(exp2.getTermNumerator());
	exp1.setDenomTerm(exp2.getDenomTerm());
	exp1.setTermNum2(exp2.getTerm2());
	exp1.setExpTerm(exp2.getExpTerm());
	exp1.setSphTerm(exp2.getSphTerm());
}

/*	Generates the specific term during the summation with at specific
	m,j,k, and n. Uses:   (-1)^m	(j+k)!	x^m
						          -------------
							      (j-m)!(k+m)!m!						*/
Expression generateLaguerrePolynomial(int m, int j, int k, int n)
{
	Expression solution;
	if(m == 0)
	{
		solution.setNumerator(factorial(j+k));
		solution.setDenominator(factorial(j) * factorial(k));
		return solution;
	}
	
	int numFact = factorial(j+k);
	int denom = factorial(j-m);
	denom *= factorial(k+m);
	denom *= factorial(m);
	numFact *= sphericalSign(m);
	
	solution.setNumerator(numFact);
	solution.setDenominator(denom);
	
	
	
	int xNumerator = 2;
	int xDenom = n;
	int myGCD = gcd(xNumerator, xDenom);
	xNumerator /= myGCD;
	xDenom /= myGCD;
	
	Term myNumTerm;
	
	myNumTerm.setCoeff(static_cast<int>(pow(xNumerator, m)));
	myNumTerm.setBaseStr("r");
	myNumTerm.setPower(m);
	
	Term myDenomTerm;
	
	myDenomTerm.setCoeff(static_cast<int>(pow(xDenom, m)));
	myDenomTerm.setBaseStr("a");
	myDenomTerm.setPower(m);
	
	solution.setDenomTerm(myDenomTerm);
	solution.setTermNumerator(myNumTerm);
	
	return solution;
}

/*	Takes in many things that needed to be previously calculated, then compiles all
	information to create the entire radial function								*/
vector<Expression> radialTotal(Expression inSqrRoot, Term exponentialTerm, Expression powerToL, vector<Expression> laguerreSolutions, int n, int l , int m)
{
	vector <Expression> results;
	
	for(int i = 0; i < laguerreSolutions.size() ; i ++)
	{
		//outputExpression(cout, laguerreSolutions[i]);
	}

	Expression sqrexppow;					//Combination of inSqrRoot, exponentialTerm, and powertoL aspects of the radial
	inSqrRoot.multiplyExpressions(powerToL);
	setExpression1To2(sqrexppow,inSqrRoot);
	
	sqrexppow.setExpTerm(exponentialTerm);
	
	Expression result;
	Term laguerreDenom;
	Term finalDenom;
	finalDenom.setBaseStr("a");
	double aPower = 0;
	int x = 0;
	Term finalNum;
	finalNum.setBaseStr("r");
	double rPower = 0;
	
	for(int i = 0; i < laguerreSolutions.size() ; i++)
	{
		aPower = 0;
		aPower += 1.5;
		aPower += l;
		aPower += x;
		rPower = l;
		rPower += x;
		x++;
		
		aPower = 1.5 + laguerreSolutions[i].getDenomTerm().getPower() + sqrexppow.getDenomTerm().getPower();

		laguerreSolutions[i].multiplyExpressions(sqrexppow);
		
		setExpression1To2(result,laguerreSolutions[i]);

		finalDenom.setPower(aPower);
		finalNum.setPower(rPower);

		result.setDenomTerm(finalDenom);
		result.setTermNumerator(finalNum);
		
		results.push_back(result);
	}

	return results;
}

/*	Retrieves the whole laguerre solutions for a specific combination of 
	the n and l quantum numbers		*/
vector<Expression> LaguerreSolutions(int n, int l)
{
	vector<Expression> laguerreSolutions;
	int k = 2 * l + 1;
	int j = n-l-1;
	
	laguerreSolutions.resize(j+1);
	
	for(int m = 0 ; m < laguerreSolutions.size() ; m++)
	{
		setExpression1To2(laguerreSolutions[m] , generateLaguerrePolynomial(m, j, k, n));
	}
	return laguerreSolutions;
}

/*	Computes another portion of the radial wave function	*/
Expression powerToTheLRadial(int n, int l)
{
	Expression result;
	result.setNumerator(pow(2,l));
	result.setDenominator(pow(n,l));
	
	Term numTerm;
	
	numTerm.setBaseStr("r");
	numTerm.setPower(l);
	
	result.setTermNumerator(numTerm);
	
	Term denomTerm;
	
	denomTerm.setBaseStr("a");
	denomTerm.setPower(l);
	
	result.setDenomTerm(denomTerm);
	
	return result;
}

/*	Creates the exponential term found in the radial wave function	*/
Term exponentialRadialTerm(int n)
{
	Term result;
	
	result.setBaseStr("e");
	result.setPowerStr("-r/(" +to_string(n) + "a)");
	result.setPower(1);
	
	return result;
}

/*	Multiplies two terms together	*/
Term multiplyTerms(Term multiplicand, Term multiplier)
{
	string multi = multiplicand.getBaseStr();
	string plier = multiplier.getBaseStr();
	Expression result;
	Term newTerm;
	
	if(plier == multi)
	{
		double totPower = multiplicand.getPower() + multiplier.getPower();
		string strTotPower = multiplicand.getPowerStr() + " + " + multiplier.getPowerStr();
		
		newTerm.setBaseStr(multi);
		newTerm.setPower(totPower);
		newTerm.setPowerStr(strTotPower);
		
		if(plier == "r")
		{
			result.setTermNumerator(newTerm);
		}
		else if(plier == "e")
		{
			result.setExpTerm(newTerm);
		}
		else if(plier == "sin" || plier == "cos")
		{
			result.setTermNum2(newTerm);
		}
		else
			cout << endl << "Something new, need new term.  Base: " << multi << endl;
	}
	
	return newTerm;
}

/*	Significanty more complicated than above, but the same idea for expressions
	Each individual terms needs to be carefully multiplied together				*/
Expression multiplyExpressions(Expression multiplicand, Expression multiplier)
{
	int newNumerator = multiplicand.getNumerator() * multiplier.getNumerator();
	int newDenominator = multiplicand.getDenominator() * multiplier.getDenominator();
	
	int newCoeff = 1;
	
	int newNumInSquare = multiplier.getSquareNum() * multiplicand.getSquareNum();
	int newDenomSquare = multiplier.getDenomInSquare() * multiplicand.getDenomInSquare();
	
	string newNumString = multiplier.getNumeratorStr() + " * " + multiplicand.getNumeratorStr();
	string newDenomString = multiplier.getDenominatorStr() + " * " + multiplicand.getDenominatorStr();
	
	string newNumSquareStr = multiplier.getSquareNumStr() + " * " + multiplicand.getSquareNumStr();
	string newDenomSquareStr = multiplier.getDenomInSquareStr() + " * " + multiplicand.getDenomInSquareStr();
	
	Term newTermNum;
	Term newTermNum2;
	Term newExpTerm;
	Term newDenomTerm;
	Term newSphTerm;
	
	if(multiplicand.getTermNumerator().getBaseStr() == multiplier.getTermNumerator().getBaseStr())
	{
		newTermNum.setPower(multiplicand.getTermNumerator().getPower() + multiplier.getTermNumerator().getPower());
		newTermNum.setBaseStr(multiplicand.getTermNumerator().getBaseStr());
		newTermNum.setPowerStr(multiplicand.getTermNumerator().getPowerStr() + " + " + multiplier.getTermNumerator().getPowerStr());
		
		newCoeff *= multiplicand.getTermNumerator().getCoeff();
		newCoeff *= multiplier.getTermNumerator().getCoeff();
	}
	else
		cout << "Something is wrong with Numerator Term, they are not the same";
	
	if(multiplicand.getTerm2().getBaseStr() == multiplier.getTerm2().getBaseStr())
	{
		newTermNum2.setPower(multiplicand.getTerm2().getPower() + multiplier.getTerm2().getPower());
		newTermNum2.setBaseStr(multiplicand.getTerm2().getBaseStr());
		newTermNum2.setPowerStr(multiplicand.getTerm2().getPowerStr() + " + " + multiplier.getTerm2().getPowerStr());
		
		newCoeff *= multiplicand.getTerm2().getCoeff();
		newCoeff *= multiplier.getTerm2().getCoeff();
	}
	else
		cout << "Something is wrong with Numerator Term2, they are not the same";
	
	if(multiplicand.getExpTerm().getBaseStr() == multiplier.getExpTerm().getBaseStr())
	{
		newExpTerm.setPower(multiplicand.getExpTerm().getPower() + multiplier.getExpTerm().getPower());
		newExpTerm.setBaseStr(multiplicand.getExpTerm().getBaseStr());
		newExpTerm.setPowerStr(multiplicand.getExpTerm().getPowerStr() + " + " + multiplier.getExpTerm().getPowerStr());
		
		newCoeff *= multiplicand.getExpTerm().getCoeff();
		newCoeff *= multiplier.getExpTerm().getCoeff();
	}
	else
		cout << "Something is wrong with exp Term, they are not the same";
	
	if(multiplicand.getSphTerm().getBaseStr() == multiplier.getSphTerm().getBaseStr())
	{
		newSphTerm.setPower(multiplicand.getSphTerm().getPower() + multiplier.getSphTerm().getPower());
		newSphTerm.setBaseStr(multiplicand.getSphTerm().getBaseStr());
		newSphTerm.setPowerStr(multiplicand.getSphTerm().getPowerStr() + " + " + multiplier.getSphTerm().getPowerStr());
		
		newCoeff *= multiplicand.getSphTerm().getCoeff();
		newCoeff *= multiplier.getSphTerm().getCoeff();
	}
	else
		cout << "Something is wrong with sph Term, they are not the same";
	
	if(multiplicand.getDenomTerm().getBaseStr() == multiplier.getDenomTerm().getBaseStr())
	{
		newDenomTerm.setPower(multiplicand.getDenomTerm().getPower() + multiplier.getDenomTerm().getPower());
		newDenomTerm.setBaseStr(multiplicand.getDenomTerm().getBaseStr());
		newTermNum.setPowerStr(multiplicand.getDenomTerm().getPowerStr() + " + " + multiplier.getDenomTerm().getPowerStr());
		
		newCoeff *= multiplicand.getDenomTerm().getCoeff();
		newCoeff *= multiplier.getDenomTerm().getCoeff();
	}
	else
		cout << "Something is wrong with Denominator Term, they are not the same";
	
	Expression product;
	
	product.setNumerator(newNumerator);
	product.setDenominator(newDenominator);
	product.setSquareNum(newNumInSquare);
	product.setDenomInSquare(newDenomSquare);
	product.setNumeratorStr(newNumString);
	product.setDenominatorStr(newDenomString);
	product.setDenomInSquareStr(newDenomSquareStr);
	product.setSquareNumStr(newNumSquareStr);
	
	product.setTermNumerator(newTermNum);
	product.setTermNum2(newTermNum2);
	product.setExpTerm(newExpTerm);
	product.setSphTerm(newSphTerm);
	product.setDenomTerm(newDenomTerm);
	
	return product;
}

/*	Creates the spherical term e^(imΦ) necessary for the spherical
	portion of the wave function									*/
Term phiSphericalTerm(int m)
{
	Term result;
	
	result.setBaseStr("e");
	result.setPower(m);
	result.setPowerStr("iΦ");
	result.setCoeff(1);
}

/*	Same as above, but specifically within the term class	*/
void Term :: phiSphericalTerm(int m)
{
	setBaseStr("e");
	setPower(m);
	setPowerStr("iΦ");
	setCoeff(1);
}

/*	Once the individual parts of the spherical wave function gets built, they
	get sent here and are able to actually construct the final total spherical
	wave function.																*/	
vector<Expression> totalSpherical(int sign, Expression insideSphericalBrackets, vector<Expression> associatedLegendres, Term sphTerm)
{
	vector <Expression> results;
	
	sphTerm.setCoeff(sphTerm.getCoeff() * sign);
	
	insideSphericalBrackets.setSphTerm(sphTerm);
	
	Expression result;
	Term numTerm1;
	Term numTerm2;
	
	int numInSquare;
	int denomInSquare;
	
	string numInSquareStr;
	string denomInSquareStr;
	
	for(int i = 0; i < associatedLegendres.size() ; i++)
	{
		numInSquare = insideSphericalBrackets.getSquareNum();
		denomInSquare = insideSphericalBrackets.getDenomInSquare();
		numInSquareStr = insideSphericalBrackets.getSquareNumStr();
		denomInSquareStr = insideSphericalBrackets.getDenomInSquareStr();
		
		insideSphericalBrackets.multiplyExpressions(associatedLegendres[i]);

		setExpression1To2(result, associatedLegendres[i]);
		
		result.setSquareNum(numInSquare);
		result.setDenomInSquare(denomInSquare);
		result.setSquareNumStr(numInSquareStr);
		result.setDenomInSquareStr(denomInSquareStr);
		result.setSphTerm(sphTerm);
		results.push_back(result);
	}
	
	return results;
}

/*	Used to create the total radial wave function solely from the given
	quantum numbers.  Will need to create each individual part of the wave function
	before completing.  This is what is called from takeIntegral					*/
vector<Expression> createRadial(int n, int l, int m)
{
	vector <Expression> testing;
	Expression radialSqrRoot;
	Expression powerToL;
	Term expTerm;
	vector<Expression> laguerres;
	Term sphTerm;
	setExpression1To2(radialSqrRoot, getRadialSqrRootPart(n,l));
	
	radialSqrRoot.simplifyRadical();
	radialSqrRoot.simplifyFraction();		
				
	setExpression1To2(powerToL, powerToTheLRadial(n,l));
			
	powerToL.simplifyRadical();
	powerToL.simplifyFraction();
							
	expTerm.setTerm(exponentialRadialTerm(n));
				
	expTerm.setCoeff(expTerm.getCoeff() * factorial(n+l));
	laguerres = LaguerreSolutions(n,l);
				
	testing = radialTotal(radialSqrRoot, expTerm, powerToL, laguerres, n, l, m);

	return testing;
	
}