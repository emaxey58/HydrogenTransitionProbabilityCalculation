/*
	This is the code that actually takes the integral for the given two wave functions.
	We pass the 6 quantum numbers and the needed perturbation to takeIntegral.  It first calls
	on functions from generateWaveFunctions.h to actually create the neceaary wave functions which 
	we will integrate over.
	
	It operates using a Simpson's Method Approximation in which it finds the value of the wave functions
	at each small step in the variable we are currently integrating over.  Since the wave functions are 
	seperable, we integrate over r, 𝜭, and Φ seperately and then multiply the results together.
	
	If you define a new perturbation, there are two places that you need to alter this code. First input
	the name of the perturbation into the enum defined below.  Then in perterbWaveFunction, create a new
	case in which your perturbation will actually affect the wave function.  This should be relatively 
	straghtforward: Find the term in the expression which has the same base as your new perturbation and 
	change the exponent accordingly.  If a term with the same base does not exist, this becomes more difficult.
	In this case, in generateWaveFunctions.h, you will need to change the Expression class.  Add in a new terma
	and replicate the set and get functions for that new term.  Then in getValueAtPoint, you will need to add
	in the new term to actually get factored in when we find the value at a specific point for the simpson's 
	method approximation.
*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdio.h>
#define PI 3.14159265

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

using namespace std;

//Available perturbations, make sure to add in the name for yours if you create a new one.
enum Perturbation {x,y,z, none};
						//x = rsin𝜭cosΦ, y = rsin𝜭sinΦ, z = rcos𝜭

/*	Function Prototypes	*/
Expression makedV();
const long double a = .05925E-9;
bool isImaginary = false;
double getValueAtPoint( string var, double value, Expression equation);
double takeIntegral(int n1, int l1, int m1, int n2, int l2, int m2);
double takeIntegral( string var, vector<Expression> equation1, vector<Expression> equation2);
double takeIntegral( string var, vector<Expression> equation1, vector<Expression> equation2, Perturbation pert);
double takeIntegral( string var1, string var2, vector<Expression> equation1, vector<Expression> equation2, double ending1, double ending2);
double takeIntegral( string var1, vector<Expression> equation1, vector<Expression> equation2, double ending1);
Expression parseForVar2(Expression equation, string var);
vector<Expression> parseForVar2(vector<Expression> equation, string var);

/*	This takes in an equation and the desired perturbation, changing the wave function.
	The equation will be altered according to the needs of the perturbation, this is one of the 
	places that you will need to alter if you add in your perturbation.			*/
void perterbWaveFunction(Expression& equation, Perturbation pert, bool fromSph = false)
{
	bool wasCos = false;
	Term newTerm;
	switch(pert)
	{
		case x:
			if(!fromSph)
			{
				newTerm = equation.getTermNumerator();
			
				if (newTerm.getBaseStr() != "r")
					newTerm.setBaseStr("r");
				
				newTerm.setPower( newTerm.getPower() + 1);
				
				equation.setTermNumerator(newTerm);
			}
			
			
			newTerm = equation.getTerm2();
			
			
			if (newTerm.getBaseStr() != "sin𝜭")
				newTerm.setBaseStr("sin𝜭");
			
			newTerm.setPower( newTerm.getPower() + 1);
			
			
			equation.setTermNum2(newTerm);
			
			newTerm = equation.getSoidalPhiTerm();

			wasCos = true;
			if (newTerm.getBaseStr() != "cosΦ")
			{
				
				newTerm.setBaseStr("cosΦ");
				wasCos = false;
			}
			
			if(wasCos)
				newTerm.setPower( newTerm.getPower() + 1);
			
			equation.setSoidalPhiTerm(newTerm);
			
			break;
		
		case y:
			if(!fromSph)
			{
				newTerm = equation.getTermNumerator();
				
				if (newTerm.getBaseStr() != "r")
					newTerm.setBaseStr("r");
				
				newTerm.setPower( newTerm.getPower() + 1);
				
				equation.setTermNumerator(newTerm);
			}
			
			newTerm = equation.getTerm2();
			
			wasCos = true;
			if (newTerm.getBaseStr() != "sin𝜭")
			{
				wasCos = false;
				newTerm.setBaseStr("sin𝜭");
			}
			if(wasCos)
				newTerm.setPower( newTerm.getPower() + 1);
			
			equation.setTermNum2(newTerm);
			
			newTerm = equation.getSoidalPhiTerm();
			
			wasCos = true;
			if (newTerm.getBaseStr() != "sinΦ")
			{
				wasCos = false;
				newTerm.setBaseStr("sinΦ");
			}
			if(wasCos)
				newTerm.setPower( newTerm.getPower() + 1);
			
			equation.setSoidalPhiTerm(newTerm);
			break;
			
		case z:
			
			if(!fromSph)
			{
				newTerm = equation.getTermNumerator();
				
				if (newTerm.getBaseStr() != "r")
					newTerm.setBaseStr("r");
				
				newTerm.setPower( newTerm.getPower() + 1);
				
				equation.setTermNumerator(newTerm);
			}
			if(fromSph)
			{
				newTerm = equation.getTermNumerator();
			
				wasCos = true;
				if (newTerm.getBaseStr() != "cos𝜭")
				{
					wasCos = false;
					newTerm.setBaseStr("cos𝜭");
				}
				if(wasCos)
					newTerm.setPower( newTerm.getPower() + 1);
				
				equation.setTermNumerator(newTerm);
			}
			
			
			break;
	}
	
}

/* 
	Used as overflow when a wave funciton has multiple expressions associated with it
	This will simply add the result of the different terms together.
*/
double getValueAtPoint( string var, double value, vector<Expression> equation)
{
	double total = 0;
	for( int i = 0 ; i < equation.size() ; i++)
	{
		total += getValueAtPoint(var, value, equation[i]);
	}
	return total;
}

/*	Finds the value of the given wave function when a variable var = value.
	Used for the simpson's method approximation, any necessary perturbations
	should have already been applied. Will essentially ignore the other 
	variables and only use the one that we are currently integrating over.
	Runs through each coefficient in the expression and then each individual
	term that the expression has to compute.  If you need to add in a new term
	it is ESSENTIAL to change this accordingly.									*/
double getValueAtPoint( string var, double value, Expression equation)
{
	long double num = 1;
	long double denom = 1;
	
	num *= equation.getNumerator();
	 int multiplier = 1;
	if(equation.getSquareNum() < 0)
	{
		multiplier = -1;
		isImaginary = true;
	}
		num *= sqrt(multiplier * equation.getSquareNum());
	
	num *= equation.getTermNumerator().getCoeff();
	num *= equation.getTerm2().getCoeff();
	num *= equation.getExpTerm().getCoeff();
	num *= equation.getSphTerm().getCoeff();
	num *= equation.getSoidalPhiTerm().getCoeff();
	
	denom *= equation.getDenominator();
	if(equation.getDenomInSquare() < 0)
	{
		multiplier = -1;
		isImaginary = true;
	}
	denom *= sqrt(multiplier * equation.getDenomInSquare());
	
	string foundp;

	string denomStr = equation.getDenomInSquareStr();
	foundp = denomStr.find("pi");
	if(denomStr.find("pi") != std::string::npos)
		denom *= sqrt(PI);

	Term  numTerm = equation.getTermNumerator();
	int power = equation.getTermNumerator().getPower();

	string termBase = equation.getTermNumerator().getBaseStr();
	
	if( termBase == "r")
	{
		if(var == "r")
		{
			num *= pow(value, numTerm.getPower());
		}
	}
	else if(termBase == "cos𝜭")
	{
		if(var == "𝜭")
		{
			num *= pow(cos ( value ), equation.getTermNumerator().getPower());
		}
	}
	else if(termBase == "")
	{
		//its empty do nothing
	}
	else
	{
		cout << "\nTerm numerator is something other than expeccted, with base = " << termBase << endl;
	}
	
	termBase = equation.getTerm2().getBaseStr();
	
	if(termBase == "")
	{
		//its empty do nothing
	}
	else if(termBase == "sin𝜭")
	{
		if(var == "𝜭")
		{
			num *= pow(sin ( value ), equation.getTerm2().getPower());
		}
	}
	else if(termBase == "cos𝜭")
	{
		if(var == "𝜭")
		{
			num *= pow(cos ( value ), equation.getTerm2().getPower());
		}
	}
	else
	{
		cout << "\nTerm numerator2 is something other than expeccted, with base = " << termBase << endl;
	}
	
	termBase = equation.getSoidalPhiTerm().getBaseStr();
	
	if(termBase == "")
	{
		//its empty do nothing
	}
	else if(termBase == "sinΦ")
	{
		if(var == "Φ")
		{
			num *= pow(sin ( value ), equation.getSoidalPhiTerm().getPower());
		}
	}
	else if(termBase == "cosΦ")
	{
		if(var == "Φ")
		{
			num *= pow(cos ( value ), equation.getSoidalPhiTerm().getPower());
		}
	}
	else
	{
		cout << "\nSoidal Phi Term is something other than expeccted, with base = " << termBase << endl;
	}
	
	termBase = equation.getDenomTerm().getBaseStr();
	
	if(termBase == "a")
	{
		denom *= pow(a, equation.getDenomTerm().getPower());
	}
	else if(termBase == "")
	{
		//It's empty nothing
	}
	else
	{
		cout << "\nDenom Term is something other than expeccted, with base = " << termBase << endl;
	}
	
	termBase = equation.getExpTerm().getBaseStr();
	
	if( termBase == "e")
	{
		string expPower = equation.getExpTerm().getPowerStr();
		long double power = -1;
		
		int i = 0;
		string powerNum;
		string powerDenom;
		
		while(expPower[i] != '/')
		{
			powerNum += expPower[i];
			i++;
		}
		
		powerDenom = expPower.substr(i+1);
		
		string numbers = "0123456789";
		char letter;
		for(int i = 0 ; i < powerNum.size() ; i++)
		{
			letter = powerNum[i];
			if(letter == '-')
			{
				//doNothing, negative is taken care off_type
				cout <<"";
			}
			else if(letter == var[0])
				power *= value;
			else if(numbers.find(letter) != -1)//if it is a digit 0-9
			{
				power *= (letter - '0');
			}
		}
		
		long double denomPower = 1;
		
		for(int i = 0 ; i < powerDenom.size() ; i++)
		{
			letter = powerDenom[i];
			if(letter == '(' || letter == ')')
			{
				//doNothing, () is unimportant
				cout << "";
			}
			else if(numbers.find(letter) != -1)//if it is a digit 0-9
			{
				power /= (letter - '0');
				denomPower *= (letter - '0');
			}
			 else if(letter == 'a')
			 {
				denomPower *= a;
				power /= a; 
			 }
		}

		num *= exp(power);
	}
	else if(termBase == "")
	{
		//its empty do nothing
	}
	else
	{
		cout << "\nExp Term numerator is something other than expeccted, with base = " << termBase << endl;
	}
	
	termBase = equation.getSphTerm().getBaseStr();
	
	if( termBase == "e")
	{
		string expPower = equation.getExpTerm().getPowerStr();
		
		double power = 1;
		
		//cout << "expPowerStr is " << expPower << " and num, denom is " << powerNum << " and " << powerDenom<< endl;
		string numbers = "0123456789";
		char letter;
		for(int i = 0 ; i < expPower.size() ; i++)
		{
			letter = expPower[i];
			if(letter == 'i')
			{
				//doNothing, imaginary we ignore,at least for now.
				//cout <<"";
			}
			else if(letter == var[0])
				power *= value;
			else if(numbers.find(letter) != -1)//if it is a digit 0-9
			{
				power *= (letter - '0');
			}
			
		}
		
		power *= equation.getExpTerm().getPower();
		
		num *= exp(power);
	}
	else if(termBase == "")
	{
		//its empty do nothing
	}
	else
	{
		cout << "\nTerm numerator is something other than expeccted, with base = " << termBase << endl;
	}
	
	if(equation.getNumeratorStr() != "")
	{
		cout << "numerator string is " << equation.getNumeratorStr() << ", something other than nothing \n";
		exit(0);
	}
	if(equation.getDenominatorStr() != "")
	{
		cout << "denominator string is " << equation.getDenominatorStr() << ", something other than nothing \n";
		exit(0);
	}

	if(num == 0)
		return 0;
	
	return num/denom;
}

/*  Was originally used to track progress in an output format.  
	Creates a progress bar of sorts.
	Currently obsolete, but may be helpful for testing purposes	*/
void printProgress(double percentage) 
{
	return;
    int val = (int) (percentage);
    int lpad = (int) (percentage / 100 * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

/*  
	This is the funciton originally called from the main function in generateWaveFunctions.C
	It creates the waveFunctions with the specified quantum numbers, orders the perturbation 
	to be done before calling the takeIntegral function that will actually do the Simpson's
	method approximation.  It multiplies this result from the radial and spherical portions
	to get the overall integral.  Please note the slight differences in procedure for 
	spherical and radial integrations.		
*/
double takeIntegral(int n1, int l1, int m1, int n2, int l2, int m2, Perturbation pert = none)
{
	isImaginary = false;
	vector <Expression> equation1 = createRadial(n1,l1,m1);
	vector <Expression> equation2 = createRadial(n2,l2,m2);
	
	double answer = takeIntegral("r", equation1, equation2, pert);

	if(isImaginary)
		answer *= -1;
	
	vector<Expression> sphEq1 = totalSpherical(n1,l1,m1);
	vector <Expression> sphEq2 = totalSpherical(n2,l2,m2);
 
	for(int i = 0;  i < sphEq1.size() ; i++)
	{
		perterbWaveFunction(sphEq1[i], pert, true);
	}
	double spherical = takeIntegral("𝜭", "Φ", sphEq1, sphEq2, PI, 2*PI);
	
	answer *= spherical;
	
	return answer;
}

/*	This is what actually does simpson's method for the radial wave function. It
	applies the given perturbation, takes the two given wave functions and integrates
	There is a simple test to see when the values are no longer relevant and as such
	we stop the approximation.														*/
double takeIntegral( string var, vector<Expression> equation1, vector<Expression> equation2, Perturbation pert = none)
{
	vector<double> additions;
	double final = 0;
	double current = 1;
	
	double increment = .000001;//This is the step used for the approximation.  It is actually .000001a, where a is the bohr radius, so a step on the order of 10^-15
	double x = 0;
	
	/*	Obsolete right now, used for the progress bar	*/
	int iteration = 0;
	int maxIterations = 1000000000;
	int loadInterval = maxIterations / 100;
	
	/*	The following applies the r^2 from dV.
		Inherently from integrating over all space in the spherical
		coordinates, we have a dV= r^2 sin(theta) which needs to be accounted
		for which is the purpose of the following for loop.		*/
	Term newTerm;
	
	for(int i = 0;  i < equation1.size() ; i++)
	{
		newTerm = equation1[i].getTermNumerator();
		
		if (newTerm.getBaseStr() != "r")
			newTerm.setBaseStr("r");
		
		newTerm.setPower( newTerm.getPower() + 2);
		
		equation1[i].setTermNumerator(newTerm);
	}
	
	//Perturb the wave function
	for(int i = 0;  i < equation1.size() ; i++)
	{
		perterbWaveFunction(equation1[i], pert);
	}
	
	int progress = 0;
	double current2;
	int oldProg = -1;
	double totCurrent;
	
	current = getValueAtPoint(var, x*a, equation1);
	current2 = getValueAtPoint(var, x*a, equation2);

	totCurrent =1;
	
	double current3;
	double timeBelowCut = 0;

	bool relevant = true;
	double prevValue = 0;

	while(x < 100)//Arbitrary cutoff point, we assume that if we get past 100a, that the further wave funcitons will be tiny enough to be irrelevant
	{
		if(iteration % loadInterval == 0)
		{
			progress++;
			printProgress(progress);
		}
		
		current = getValueAtPoint(var, x*a, equation1);
		
		current2 = getValueAtPoint(var, x*a, equation2);

		totCurrent = current * current2; 

		if(totCurrent < 1E-5)//If the value is really freaking small and it doesnt increase, we consider it irrelevant and stop the integration.
		{
			timeBelowCut += increment;
			if(totCurrent > prevValue)
			{
				timeBelowCut = 0;
			}
			if(timeBelowCut >= 1)//1 a0
				relevant = false;
			prevValue = totCurrent;
		}
	
		additions.push_back(totCurrent);

		if(x == 0)
			totCurrent = 1;
		x += increment;
		
		iteration++;
	} 
	
	/*	
		The following is guided by Simpson's Method Approximation, which has the equation,
		increment/3(y(0) + 4y(1) + 2y(2) + ... + 2y(n-2) + 4y(n-1) + y(n))
	*/
	final += additions[0];
	final += additions[additions.size() - 1];
	
	int multiplier;

	for(int i = 0 ; i < additions.size() - 2 ; i++)
	{
		if(i % 2 == 0)
			multiplier = 2;
		else
			multiplier = 4;
		
		final += (multiplier * additions[i]);
	}
	
	final *= a;
	final *= increment;
	final /= 3;
	
	return final;
}

/*	Special case for spherical since it is one expression that is integrated over two
	seperate variables.  It takesIntegral in terms of theta normally, then uses 
	parseForVar2 to eliminate any possible redundancies that could occur with constants
	that are inherent in the expression.  Then it integrates over phi.		*/
double takeIntegral( string var1, string var2, vector<Expression> equation1, vector<Expression> equation2, double ending1, double ending2)
{
	double result;
	
	
	string firstVar = var1;
	string secVar = var2;
	vector<Expression> firstEquations = equation1;
	vector<Expression> secEquations = equation2;
	double end1 = ending1;
	double end2 = ending2;
	result = takeIntegral(firstVar, firstEquations, secEquations, end1);

	vector<Expression> new1 = parseForVar2(equation1, var2);
	vector<Expression> new2 = parseForVar2(equation2, var2);
 	double result2 = takeIntegral(secVar, new1, new2, end2);

	result *= result2;
	return result;
}

//Gets the part that only has phi
Expression parseForVar2(Expression equation, string var)
{
	Expression newExpression;
	
	newExpression.setSphTerm(equation.getSphTerm());
	
	return newExpression;
}

/*	Same as above, but with the purpose of use for when we have a vector of expressions
	and not just one simple expression	*/
vector<Expression> parseForVar2(vector<Expression> equation, string var)
{
	vector<Expression> newExpressions;
	newExpressions.resize(equation.size());
	Expression expressionI;
	for(int i = 0; i < equation.size() ; i++)
	{
		expressionI = parseForVar2(equation[i], var);
		newExpressions[i] = expressionI;
	}

	return newExpressions;
}

/*	Used to get dV with the spherical integral.	*/
Expression makedV()
{
	Expression dV;
	Term newTerm;
	newTerm.setBaseStr("r");
	newTerm.setPower(2);
	dV.setTermNumerator(newTerm);
	
	newTerm.setBaseStr("sin𝜭");
	newTerm.setPower(1);
	dV.setTermNum2(newTerm);
	
	return dV;
}

/*	Takes integral for the spherical, operates the same as the radial one, but slightly different to 
	account for changes with dV and perturbation method		*/
double takeIntegral( string var1, vector<Expression> equation1, vector<Expression> equation2, double ending1)//Spherical Integral Φ
{
	vector<double> additions;
	double final = 0;
	double current = 1;
	
	double increment = .000001;
	double x = 0;
	
	int iteration = 0;
	int maxIterations = 1000000000;
	int loadInterval = maxIterations / 100;

	//dv is r^2 sin theta dtheta dphi dr
	Expression dV = makedV();
	Term term2;
	
	if(var1 == "𝜭")
	{
		for(int i = 0 ; i < equation1.size() ; i++)
		{
			term2 = equation1[i].getTerm2();
			if(term2.getBaseStr() != "sin𝜭")
			{
				term2.setBaseStr("sin𝜭");
			}
			term2.setPower(term2.getPower() + 1);
			
			equation1[i].setTermNum2(term2);
		}
		
	}
	
	int progress = 0;
	double current2;
	double current3;
	int oldProg = -1;
	double totCurrent;
	
	current = getValueAtPoint(var1, x, equation1);
	current2 = getValueAtPoint(var1, x, equation2);

	totCurrent = current * current2;

	while(x <= ending1)//For theta, phi, we have specific endings, unlike infinity for radial
	{
		if(iteration % loadInterval == 0)
		{
			progress++;
			printProgress(progress);
		}
		
		current = getValueAtPoint(var1, x, equation1);
		current2 = getValueAtPoint(var1, x, equation2);
		
		totCurrent = current * current2; // * current3;
		
		x += increment;
		
		additions.push_back(totCurrent);
		iteration++;
	} 
	
	final += additions[0];
	final += additions[additions.size() - 1];
	
	int multiplier;
	
	//Once again, this is the implementation of simpson's method
	for(int i = 0 ; i < additions.size() - 2 ; i++)
	{
		if(i % 2 == 0)
			multiplier = 2;
		else
			multiplier = 4;
		
		final += (multiplier * additions[i]);
	}
	
	final *= increment;
	final /= 3;
	
	return final;
}