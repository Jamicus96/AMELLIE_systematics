RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

all: GetFOMabsSlope.exe GetHists.exe
	
GetFOMabsSlope.exe: GetFOMabsSlope.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y GetFOMabsSlope.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o GetFOMabsSlope.exe ${RAT_CONFIG}

GetHists.exe: GetHists.cpp
	g++ -g -std=c++1y GetHists.cpp -o GetHists.exe ${RAT_CONFIG}
	
clean :
	rm *.exe
