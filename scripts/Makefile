RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

all: GetStats.exe GetHists.exe
	
GetStats.exe: GetStats.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++11 GetStats.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o GetStats.exe ${RAT_CONFIG}

GetHists.exe: GetHists.cpp
	g++ -g -std=c++11 GetHists.cpp -o GetHists.exe ${RAT_CONFIG}
	
clean :
	rm *.exe
