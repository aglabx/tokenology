CXX=g++
CXXFLAGS=-std=c++2a -pthread -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive  -O3 -rdynamic
LDLIBS=

TARGET=bpe.exe
SRCS=nlohmann/json.hpp tokens_model.hpp bpe.cpp


all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) $(LDLIBS) -o $(TARGET)

.PHONY: clean

clean:
	rm -f $(TARGET)

# sudo apt install libtbb-dev
# g++ -std=c++2a -pthread -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -O3 -rdynamic json.hpp bpe.cpp -o bpe.exe 
