/*
Name: Kiss Sándor Ádám
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>

int randomInt(int maxValue) {
	return std::rand() % maxValue;
}

void fillArrayWithRandomNumbers(int size, int *array, int maxValue) {
	for (int i = 0; i < size; ++i) {
		array[i] = randomInt(maxValue);
	}
}

int *randomArray(int size, int maxValue) {
	int *randomNumbers = new int[size];
	
	fillArrayWithRandomNumbers(size, randomNumbers, maxValue);
	
	return randomNumbers;
}

void quickSort(int* array, int lower, int upper) {
	int pivot = array[upper];
	int i = lower;
	int j = upper;

	while (i <= j) {
		while (array[i] < pivot) {
			i++;
		}
		while (array[j] > pivot) {
			j--;
		}
		if (i <= j) {
			std::swap(array[i++], array[j--]);
		}
	}
	
	if (lower < j) {
		quickSort(array, lower, j);
	}
	if (i < upper) {
		quickSort(array, i, upper);
	}
}

void quickSort(int size,  int* array) {
	quickSort(array, 0, size - 1);
}

struct result {
	int *first;
	int firstSize;
	
	int *second;
	int secondSize;
};

struct result task8(int size, int *array) {
	struct result res;
	
	quickSort(size, array);
	
	res.first = array;
	res.firstSize = size / 2;
	
	res.second = array + (size / 2);
	res.secondSize = size / 2;
	
	return res;
}

int main () {
	std::srand(time(NULL));
	
	const int size = 10;
	const int maxValue = 100;
	
	int *array = randomArray(size, maxValue);
	
	struct result res = task8(size, array);
	std::cout << res.first << " " << res.second << std::endl;
	
	for (int i = 0; i < res.firstSize; ++i) {
		std::cout << res.first[i] << " ";
	}
	std::cout << std::endl << std::endl;
	
	for (int i = 0; i < res.secondSize; ++i) {
		std::cout << res.second[i] << " ";
	}
	std::cout << std::endl << std::endl;
	
	delete[] array;
	
	return EXIT_SUCCESS;
}
