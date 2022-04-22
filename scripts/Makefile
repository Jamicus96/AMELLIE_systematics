RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

all: GetFOM.exe optimiseRegions.exe getTrackingInfo.exe
	
GetFOM.exe: GetFOM.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y GetFOM.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o GetFOM.exe ${RAT_CONFIG}

optimiseRegions.exe: optimiseRegions.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y optimiseRegions.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o optimiseRegions.exe ${RAT_CONFIG}

getTrackingInfo.exe: getTrackingInfo.cpp
	g++ -g -std=c++1y getTrackingInfo.cpp -o getTrackingInfo.exe ${RAT_CONFIG}
	
clean :
	rm *.exe
