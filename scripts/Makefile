RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

all: GetFOMabsSlope.exe getTrackingInfo.exe
	
GetFOMabsSlope.exe: GetFOMabsSlope.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y GetFOMabsSlope.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o GetFOMabsSlope.exe ${RAT_CONFIG}

getTrackingInfo.exe: getTrackingInfo.cpp
	g++ -g -std=c++1y getTrackingInfo.cpp -o getTrackingInfo.exe ${RAT_CONFIG}
	
clean :
	rm *.exe
