all: optimiseRegions.exe getTrackingInfo.exe
	
optimiseRegions.exe: optimiseRegions.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y optimiseRegions.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o optimiseRegions.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

getTrackingInfo.exe: getTrackingInfo.cpp
	g++ -g -std=c++1y getTrackingInfo.cpp -o getTrackingInfo.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
	
clean :
	rm *.exe
