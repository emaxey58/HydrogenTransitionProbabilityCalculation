//Emily Maxey
//Calculate Transition Probabilities
//Last updated 03/23/22

/*	As part of an undergraduate Faculty Mentored Research Grant at 
	Angelo State University, the probability of transition probabilities
	between various wave functions of hydrogen after light was absorbed
	was computationally generated.  The code to undergo said process is
	below.  
	
	The general schematic works as follows:
	
	main requests the integrals with two specific combinations of quantum numbers
			|
			V
	generateWaveFunctions.h generates the wave functions
			|
			V
	takeIntegral.h perturbs the functions before using a simpson's method approximation
			|
			V
	The result of the integral is returned to main which outputs it in a formatted sense.
	
	This code is free for use!  But please credit me if you do choose to 
	use it.:) It is currently set up for the transition rates due to a light
	perturbation, but reference takeIntegral.h if you choose or need to 
	implement a new perturbation. For assistance in use or if you have any
	questions, please contact me at emilymaxey58@gmail.com
*/

#include "generateWaveFunctions.h"
#include "takeIntegral.h"
#include <thread>

//const long double a = .05925E-9;
/*	Main Function	*/

void testFunction(int n, int l)
{
	//cout << "\nn is " << n;
	//cout << "\nl is " << l;
}

bool isValidM( int m, int l)
{
	return (m >= 0 && m <= l);
}

bool isValidL(int l, int n)
{
	return (l >= 0 && l < n);
}

void outputtingWaveFunctions(int n1,int l1, int m1)
{
	
}

void testOutput(int n1,int l1, int m1, int n2, int l2, int m2)
{
	cout << endl << n1 << l1 << m1 << " to " << n2 << l2 << m2 << endl;
	vector <Expression> equation1 = createRadial(n1,l1,m1);
	vector <Expression> equation2 = createRadial(n2,l2,m2);
	
	for(int i = 0 ; i < equation1.size() ; i++)
	{
		outputExpression(cout, equation1[i]);
	}
	for(int i = 0 ; i < equation2.size() ; i++)
	{
		outputExpression(cout, equation2[i]);
	}
	
	
	vector<Expression> sphEq1 = totalSpherical(n1,l1,m1);
	vector <Expression> sphEq2 = totalSpherical(n2,l2,m2);
	
	for(int i = 0 ; i < sphEq1.size() ; i++)
	{
		outputExpression(cout, sphEq1[i]);
	}
	for(int i = 0 ; i < sphEq2.size() ; i++)
	{
		outputExpression(cout, sphEq2[i]);
	}
}

void testOutput(int n1,int l1, int m1)
{
	cout << endl << n1 << l1 << m1 << endl;
	vector <Expression> equation1 = createRadial(n1,l1,m1);
	
	vector<Expression> sphEq1 = totalSpherical(n1,l1,m1);
	
	cout << "Radial for n=" << n1 << ", l=" << l1 << ", m=" << m1 << endl;
	for(int i = 0 ; i < equation1.size() ; i++)
	{
		outputExpression(cout, equation1[i]);
	}

	cout << "Spherical for n=" << n1 << ", l=" << l1 << ", m=" << m1 << endl;
	for(int i = 0 ; i < sphEq1.size() ; i++)
	{
		outputExpression(cout, sphEq1[i]);
	}
}

int main(int argc,char** argv)
{
	vector <Expression> testing;
	Expression radialSqrRoot;
	Expression powerToL;
	Term expTerm;
	vector<Expression> laguerres;
	Term sphTerm;
	
	if(argv[1] == NULL)
		exit(0);
	
	vector<Expression> equation1;
	
	vector<Expression> equation2;
	
	equation1 = totalSpherical(2,1,0);
	equation2 = totalSpherical(2,1,0);
	
	double myIntegral;
	
	equation2 = equation1;
	bool statisticallySignificant = true;
	double minProb = 1E-20;

	int n2,l2,m2;
	string output = "";

	int n = stoi(argv[1]);
	for(int l = 0; l < n ; l++)
	{
		for(int m = 0 ; m <= l ; m++)
		{
			n2 = n;
			statisticallySignificant = true;
			
			while(statisticallySignificant)
			{
	
				n2--;
				l2 = l+1;
				if(isValidL(l2,n2))
				{
					for(int mChange = -1 ; mChange <= 1 ; mChange++)
					{
						m2 = m + mChange;
						if(isValidM(m2, l2))
						{
							string i;
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, x);
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tx" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, y);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\ty" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, z);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tz" << endl;
						}
					}
					
				}
				
				l2 = l-1;
				if(isValidL(l2, n2))
				{
					for(int mChange = -1 ; mChange <= 1 ; mChange++)
					{
						m2 = m + mChange;
						if(isValidM(m2, l2))
						{
							string i;
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, x);
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tx" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, y);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\ty" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, z);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tz" << endl;
						}
					}
				}
				
				if(n2 <= 0)
					statisticallySignificant = false;
				
			}
		}
	}
	exit(0);
	for(int l = 0; l < n ; l++)
	{
		for(int m = 0 ; m <= l ; m++)
		{
			n2 = n;
			statisticallySignificant = true;
			
			while(statisticallySignificant)
			{
	
				n2++;
				l2 = l+1;
				if(isValidL(l2,n2))
				{
					for(int mChange = -1 ; mChange <= 1 ; mChange++)
					{
						m2 = m + mChange;
						if(isValidM(m2, l2))
						{
							string i;
							//cout << n  << " , " << l << " , " << m << " , "<< n2 << " , "<< l2 << " , "<< m2<< " , " << endl;
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, x);
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tx" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, y);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\ty" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, z);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tz" << endl;
						
						}
					}
					
				}
				
				l2 = l-1;
				if(isValidL(l2, n2))
				{
					for(int mChange = -1 ; mChange <= 1 ; mChange++)
					{
						m2 = m + mChange;
						if(isValidM(m2, l2))
						{
							string i;
							
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, x);
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tx" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, y);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\ty" << endl;
							
							i = "";
							if (myIntegral < 0)
							{
								myIntegral*= -1;
								i ="i";
							}
							myIntegral = takeIntegral(n,l,m, n2, l2, m2, z);
							cout << setw(5) << n << "\t" << setw(5) <<  l << "\t" << setw(5) << m << "\t" << setw(5) << n2 << "\t" <<  setw(5) << l2 << "\t" << setw(5) << m2 << "\t" << setw(5) << myIntegral << i << "\tz" << endl;
						
						}
					}
				}
				
				if(n2 > n+7)
					statisticallySignificant = false;
				
			}
		}
	}
	
	return 0;
}