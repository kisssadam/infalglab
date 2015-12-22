/*
developer: Kiss Sándor Ádám
e-mail: kisssandoradam@gmail.com
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stack>

class Node {
public:
	Node(int value = 0, Node *left = NULL, Node *right = NULL) {
		this->value = value;
		this->left = left;
		this->right = right;
	}

	int getValue() {
		return value;
	}

	void setValue(int value) {
		this->value = value;
	}

	Node* getLeft() {
		return left;
	}

	void setLeft(Node* left) {
		this->left = left;
	}

	Node* getRight() {
		return right;
	}

	void setRight(Node* right) {
		this->right = right;
	}

	int countChildren() {
		int childrenCount = 0;

		if (this->left != NULL) {
			++childrenCount;
		}

		if (this->right != NULL) {
			++childrenCount;
		}

		return childrenCount;
	}

private:
	int value;
	Node *left;
	Node *right;
};

std::ostream& operator<<(std::ostream& os, Node* node) {
	os << "Node: [value=" << node->getValue();
	
	os << " left=";
	if (node->getLeft() == NULL) {
		os << "NULL";
	} else {
		os << node->getLeft()->getValue();
	}

	os << " right=";
	if (node->getRight() == NULL) {
		os << "NULL";
	} else {
		os << node->getRight()->getValue();
	}

	return os << "]";
}

class BinarySearchTree{
public:
	BinarySearchTree() {
		this->root = NULL;
	}

	~BinarySearchTree() {
		postorderTraversal(root, deleteNode);
	}

	static void deleteNode(Node *node) {
		delete node;
	}

	bool insertNode(int value) {
		if (root == NULL) {
			root = new Node(value);
		} else {
			Node *p = root;
			Node *parent = NULL;
			
			while (p != NULL) {
				if (p->getValue() == value) {
					return false;
				} else {
					parent = p;

					if (value < p->getValue()) {
						p = p->getLeft();
					} else {
						p = p->getRight();
					}
				}
			}

			if (value < parent->getValue()) {
				parent->setLeft(new Node(value));
			} else {
				parent->setRight(new Node(value));
			}
		}

		return true;
	}

	bool deleteNode(int value) {
		Node *p = root;
		Node *parent = NULL;

		while ((p != NULL) and p->getValue() != value) {
			parent = p;

			if (value < p->getValue()) {
				p = p->getLeft();
			} else {
				p = p->getRight();
			}
		}

		if (p == NULL) {
			return false;
		} else {
			int childrenCount = p->countChildren();

			if (childrenCount == 0) {
				deleteLeafNode(p, parent);
			} else if (childrenCount ==1 ) {
				deleteNodeWithOneChild(p, parent);
			} else {
				Node *q = p->getRight();
				Node *t = p;

				while (q->getLeft() != NULL) {
					t = q;
					q = q->getLeft();
				}
				
				p->setValue(q->getValue());
				
				if (q->getRight() == NULL) {
					deleteLeafNode(q, t);
				} else {
					deleteNodeWithOneChild(q, t);
				}
			}

			return true;
		}
	}

	Node *findNode(int value) {
		Node *p = root;
		
		while ((p != NULL) and (p->getValue() != value)) {
			p = (value < p->getValue()) ? p->getLeft() : p->getRight();
		}

		return p;
	}

	void preorderTraversal() {
		preorderTraversal(printNode);
	}

	void inorderTraversal() {
		inorderTraversal(printNode);
	}

	void postorderTraversal() {
		postorderTraversal(printNode);
	}
	
	static void printNode(Node *node) {
		std::cout << node << std::endl;
	}

	void preorderTraversal(void (*processNodeFunction)(Node*)) {
		preorderTraversal(root, processNodeFunction);
	}

	void inorderTraversal(void (*processNodeFunction)(Node*)) {
		inorderTraversal(root, processNodeFunction);
	}

	void postorderTraversal(void (*processNodeFunction)(Node*)) {
		postorderTraversal(root, processNodeFunction);
	}

	void insertUniqueRandomValues(const int size, const int maxValue) {
		std::cout << "Inserting " << size << " unique values between [0, " << maxValue << ") into the binary search tree: " << std::endl;

		for (int i = 0; i < size; ++i) {
			int randomNumber = rand() % maxValue;

			if (this->findNode(randomNumber) == NULL) {
				std::cout << randomNumber << (i + 1 == size ? "" : " ");
				if (this->insertNode(randomNumber) == false) {
					std::cout << "Failed to insert value: " << randomNumber << ", because the tree already contains it." << std::endl;
					std::cout << "You shouldn't see this message, because findNode() method is called just before the insertion." << std::endl;
				}
			} else {
				--i;
			}
		}

		std::cout << std::endl;
	}

	void deleteUniqueRandomValues(const int uniqueRandomNumbersCount, const int maxRandomValue) {
		for (int i = 0; i < uniqueRandomNumbersCount; ++i) {
			int randomNumber = rand() % maxRandomValue;
			
			bool deleteResult = this->deleteNode(randomNumber);
			std::cout << "Deleting value " << std::setw(4) << randomNumber << ": " << (deleteResult == true ? "success" : "fail (wasn't in the tree)") << std::endl;
		}
	}

	void deleteEveryNode() {
		while(root != NULL) {
			this->deleteNode(root->getValue());
		}
	}

	int getMaxValue() {
		Node *p = root;
		
		while (p->getRight() != NULL) {
			p = p->getRight();
		}

		return p->getValue();
	}

	int getMinValue() {
		Node *p = root;
		
		while (p->getLeft() != NULL) {
			p = p->getLeft();
		}

		return p->getValue();
	}

	Node *getRoot() {
		return root;
	}

private:
	Node *root;

	void deleteLeafNode(Node *node, Node *parentNode) {
		if (parentNode == NULL) {
			root = NULL;
		} else {
			if (parentNode->getLeft() == node) {
				parentNode->setLeft(NULL);
			} else {
				parentNode->setRight(NULL);
			}
		}

		delete node;
	}

	void deleteNodeWithOneChild(Node *node, Node *parentNode) {
		Node *temp = (node->getLeft() != NULL) ? node->getLeft() : node->getRight();

		if (parentNode == NULL) {
			root = temp;
		} else {
			if (parentNode->getLeft() == node) {
				parentNode->setLeft(temp);
			} else {
				parentNode->setRight(temp);
			}
		}

		delete node;
	}

	void preorderTraversal(Node *node, void (*processNodeFunction)(Node*)) {
		if (node == NULL) {
			return;
		} else {
			processNodeFunction(node);
			postorderTraversal(node->getLeft(), processNodeFunction);
			postorderTraversal(node->getRight(), processNodeFunction);
		}
	}

	void inorderTraversal(Node *node, void (*processNodeFunction)(Node*)) {
		std::stack<Node*> *parentStack = new std::stack<Node*>();

		while (!parentStack->empty() or node != NULL) {
			if (node != NULL) {
				parentStack->push(node);
				
				node = node->getLeft();
			} else {
				node = parentStack->top();
				parentStack->pop();
				
				processNodeFunction(node);
				
				node = node->getRight();
			}
		}

		delete parentStack;
	}

	void postorderTraversal(Node *node, void (*processNodeFunction)(Node*)) {
		if (node == NULL) {
			return;
		} else {
			postorderTraversal(node->getLeft(), processNodeFunction);
			postorderTraversal(node->getRight(), processNodeFunction);
			processNodeFunction(node);
		}
	}
};

int main() {
	std::srand(time(NULL));

	BinarySearchTree *binarySearchTree = new BinarySearchTree();

	const int size = 1000;
	const int maxRandomValue = 10000;
	
	binarySearchTree->insertUniqueRandomValues(size, maxRandomValue);
	std::cout << std::endl;

	int minValue = binarySearchTree->getMinValue();
	int maxValue = binarySearchTree->getMaxValue();
	
	std::cout << "The minimum value has been inserted is: " << minValue << std::endl;
	std::cout << "The maximum value has been inserted is: " << maxValue << std::endl;
	std::cout << std::endl;

	int average = (minValue + maxValue) / 2;
	Node *node = binarySearchTree->findNode(average);

	std::cout << "The average (" << average << ") of minimum and maximum is";
	if (node == NULL) {
		std::cout << "n't";
	}
	std::cout << " in the binary search tree." << std::endl;
	std::cout << std::endl;

	binarySearchTree->deleteUniqueRandomValues(10, maxRandomValue);
	std::cout << std::endl;

	std::cout << "Deleting every node in the tree." << std::endl;
	binarySearchTree->deleteEveryNode();
	std::cout << ((binarySearchTree->getRoot() == NULL) ? "Every node has been successfully deleted." : "Failed to delete every node in the tree.") << std::endl;
	
	delete binarySearchTree;

	return EXIT_SUCCESS;
}
