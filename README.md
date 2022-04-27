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
