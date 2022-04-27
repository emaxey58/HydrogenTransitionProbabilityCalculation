#!/bin/bash

#Kill commands:
#ps -u emaxey
#Will list all current commands: We need the PID number
#kill PID#
#ex: kill 400

#in putty: 
# dos2unix genWaveFunctions.bash
# bash genWaveFunctions.bash

valid=true
num = 1

while (( num < 7 )); do			#Will run through n = 6 at the max

	nice nice ./generateWaveFunctions $num >> integrals02-14-take2.txt &		#Make sure to replave with your own output file pathing
																				#Nice nice used to not overload system/take all processing power, may change with your system
	echo $num	#outputs num
	((num++))
done

