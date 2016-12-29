
INC = $(wildcard *.h)
SRC = $(wildcard *.cpp)
OPT = -Wall -g -Ofast -std=c++11

LIB += -lOpenMeshTools -lOpenMeshCore
ifeq ($(OS),Windows_NT)
	LIB += -Wl,-Bstatic -lglew32 -Wl,-Bdynamic -lglfw3 -lopengl32 -lglu32 -lgdi32
else
	LIB += -lGLEW -lglfw3 -lpthread -lGL -lXrandr -lGLU -lX11 -lXxf86vm -lXi -ldl -lXinerama -lXcursor
endif

%.o : %.cpp $(INC)
	g++ $< -c -o $@ $(OPT)
	
all: $(SRC:.cpp=.o)
	g++ -o demo $^ $(OPT) $(LIB)

clean:
	rm -f *.o
