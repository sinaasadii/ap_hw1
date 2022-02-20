CXX = g++
CXXFLAGS = -std=c++23 -Wall -I h -I /usr/local/include/gtest/ -c
LXXFLAGS = -std=c++23 -Ih -pthread
OBJECTS = ./obj/hw1.o ./obj/main.o ./obj/unit_test.o
GTEST = /usr/local/lib/libgtest.a
TARGET = main


$(TARGET): $(OBJECTS)
	$(CXX) $(LXXFLAGS) -o $(TARGET) $(OBJECTS) $(GTEST)
./obj/hw1.o: ./cpp/hw1.cpp
	$(CXX) $(CXXFLAGS) ./cpp/hw1.cpp -o ./obj/hw1.o
./obj/unit_test.o: ./cpp/unit_test.cpp
	$(CXX) $(CXXFLAGS) ./cpp/unit_test.cpp -o ./obj/unit_test.o
./obj/main.o: ./cpp/main.cpp
	$(CXX) $(CXXFLAGS) ./cpp/main.cpp -o ./obj/main.o
clean:
	rm -fv $(TARGET) $(OBJECTS)
