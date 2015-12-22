CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic-errors -march=native

EXERCISES=\
	exercise-01 \
	exercise-02 \
	exercise-03 \
	exercise-04 \
	exercise-05 \
	exercise-06 \
	exercise-07 \
	exercise-08 \
	exercise-09 \
	exercise-10 \
	exercise-11 \
	exercise-written-examination

all: $(EXERCISES)

$(EXERCISES): %: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(EXERCISES)
