#include<iostream>
#include<bits/stdc++.h>
using namespace std;

int T;
double minSupport;
double minConfidence;
vector<vector<int> > transactions;
vector<vector<int> > sortedTransactions;
vector<vector<int> > L;
vector<int> Lfreq;
vector<vector<vector<int> > > frequentSet;
vector<vector<int> > frequentSetCount;

// FP Tree Node structure
typedef struct FPnode {
	int item;
	int freq;
	struct FPnode* parent;
	struct FPnode* next;
	vector<struct FPnode*> children;
}FPnode;
FPnode* F;

// Conditional FP Tree Node structure
typedef struct condFPnode {
	int item;
	int freq;
	vector<struct condFPnode*> children;
}condFPnode; 

// FP Tree Node structure
typedef struct header {
	int item;
	int freq;
	struct FPnode* next;
}header;
vector<header> Header;

void printVector(vector<int> T) {
	for(int i = 0; i < T.size(); ++i)
		cout<<T[i] <<" ";
}

void print(vector<vector<int> > T) {
	for(int i = 0; i < T.size(); ++i) {
		for(int j = 0; j < T[i].size(); ++j)
			cout<<T[i][j] <<" ";
		cout<<"\n";
	}
	cout<<"\n";
}

void printHeader() {
	for(int i = 0; i < Header.size(); ++i) {
		FPnode* I;
		cout<<"Item: " <<Header[i].item <<" Freq: " <<Header[i].freq <<" Next-> ";
		if(Header[i].next == NULL)
			cout<<"NULL";
		else {
			I = Header[i].next;
			cout<<I->item <<": " <<I->freq <<", ";
			while(I->next) {
				I = I->next;
				cout<<I->item <<": " <<I->freq <<", ";
			}
		}
		cout<<"\n";
	}
}

// Prints Tree in Level Order
void printTree() {
	int i;
	queue<FPnode*> Q;
	FPnode* T = F;
	FPnode* tmp = new(FPnode);
	tmp->freq = -1;
	tmp->item = -10;
	tmp->next = NULL;
	tmp->parent = NULL;
	Q.push(tmp);
	while(1) {
		if(T->item == -10) {
			cout<<"\n";
			if(Q.empty())
				return;
			Q.push(tmp);
			T = Q.front();
			Q.pop();
			continue;
		}
		else {
			cout<<T->item <<"(";
			if(T->parent == NULL)
				cout<<"N";
			else
				cout<<T->parent->item;
			cout<<"):" <<T->freq <<"  ";
		}
		for(i = 0; i < T->children.size(); ++i)
			Q.push(T->children[i]);
		T = Q.front();
		Q.pop();
	}
}

// Prints Conditional Tree in Level Order
void printCondTree(condFPnode* F) {
	int i;
	queue<condFPnode*> Q;
	condFPnode* T = F;
	condFPnode* tmp = new(condFPnode);
	tmp->freq = -1;
	tmp->item = -10;
	Q.push(tmp);
	while(1) {
		if(T->item == -10) {
			cout<<"\n";
			if(Q.empty())
				return;
			Q.push(tmp);
			T = Q.front();
			Q.pop();
			continue;
		}
		else {
			cout<<T->item <<"(";
			cout<<"):" <<T->freq <<"  ";
		}
		for(i = 0; i < T->children.size(); ++i)
			Q.push(T->children[i]);
		T = Q.front();
		Q.pop();
	}
}

bool compareHeaders(header A, header B) {
	return A.freq > B.freq;
}

int findItemInHeader(int item) {
	for(int i = 0; i < Header.size(); ++i) {
		if(Header[i].item == item)
			return i;	
	}
	return -1;
}

bool compareItems(int a, int b) {
	int afreq = Header[findItemInHeader(a)].freq; 
	int bfreq = Header[findItemInHeader(b)].freq;
	return afreq > bfreq;
}

// Remove infrequent items from Header and Transactions then sort both in descending freq order
void cleanUpData() {
	int freq, i, j;
	double sup;
	// Sort Header and then transactions in descending order of item freq
	sort(Header.begin(), Header.end(), compareHeaders);
	for(i = 0; i < sortedTransactions.size(); ++i)
		sort(sortedTransactions[i].begin(), sortedTransactions[i].end(), compareItems);
	// Remove infrequent items from Header
	freq = Header[Header.size()-1].freq;
	sup = ((double) freq/T)*100;
	while(sup < minSupport) {
		Header.pop_back();
		freq = Header[Header.size()-1].freq;
		sup = ((double) freq/T)*100;
	}
	// Remove Infrequent Items from Transactions
	for(i = 0; i < sortedTransactions.size(); ++i) {
		for(j = sortedTransactions[i].size()-1; j >= 0; --j) {
			if(findItemInHeader(sortedTransactions[i][j]) == -1)
				sortedTransactions[i].pop_back();
			else
				break;
		}
	}
}

// Get Transactions from file and store in vector<vector<int> > transactions
void getTransactions(char* inputFile) {
	fstream transactionFile;
	transactionFile.open(inputFile, ios::in);
	int i, Iid, itemSize;
	header tmpHeader;
	tmpHeader.item = -1;
	tmpHeader.freq = 0;
	tmpHeader.next = NULL;
	vector<int> tmp;
	
	while(!transactionFile.eof()) {
		transactionFile>>Iid >>Iid;
		tmp.clear();
		while(Iid != -1) {
			itemSize = Header.size();
			if(Iid > itemSize) {
				for(i = itemSize; i < Iid; ++i)
					Header.push_back(tmpHeader);
			}
			Header[Iid-1].item = Iid;
			Header[Iid-1].freq++;
			tmp.push_back(Iid);
			transactionFile>>Iid;
		}
		transactions.push_back(tmp);
		sortedTransactions.push_back(tmp);
	}
	T = transactions.size();
	transactionFile.close();
	cleanUpData();
}

// Get support of itemset
double getFrequency(vector<int> itemset) {
	int freq = 0, i, j, k;
	for(i = 0; i < transactions.size(); ++i) {
		k = 0;
		for(j = 0; j < transactions[i].size(); ++j) {
			if(itemset.size() > transactions[i].size())
				break;
			if(k == itemset.size()) 
				break;
			if(itemset[k] == transactions[i][j])
				k++;
		}
		if(k == itemset.size()) 
			freq++;
	}
	return freq;
}

// Return complement of subset S
vector<int> getSubsetComplement(vector<int> itemset, vector<int> subset) {
	vector<int> comp;
	int p = 0, m;
	for(m = 0; m < itemset.size(); ++m) {
		if(p == subset.size()) {
			comp.push_back(itemset[m]);
			continue;
		}
		if(subset[p] == itemset[m]) {
			p++;
			continue;
		}
		else {
			comp.push_back(itemset[m]);
			continue;
		}
	}
	return comp;
}

// Generate subsets of length k from given itemset
vector<vector<int> > generateKlengthSubset(vector<int> itemset, int k) {
	vector<vector<int> > res;
	vector<int> s;
	int i, j, l;
	for(i = 0; i < itemset.size()-k+1; ++i) {
		l = 0;
		while(l != k-1) {
			s.push_back(itemset[i+l]);
			l++;
		}
		for(j = i+l; j < itemset.size(); ++j) {
			s.push_back(itemset[j]);
			l++;
			if(l == k) {
				res.push_back(s);
				s.pop_back();
				l--;
			}
		}
		s.clear();
		if(k == 1)
			break;
	}
	return res;
}

// Generate Association Rules
void generateAssociationRules() {
	fstream outputFile;
	outputFile.open("association.txt", ios::out);
	vector<vector<int> > s;
	vector<int> sComplement;
	double conf, supS, supL, freq;
	int i, j, k, l, m, p, ktmp, ltmp;
	
	for(i = 1; i < frequentSet.size(); ++i) {
		for(j = 0; j < frequentSet[i].size(); ++j) {
			printVector(frequentSet[i][j]);
			for(k = 1; k < frequentSet[i][j].size(); ++k) {
				cout<<"\n" <<k <<":\n";
				s = generateKlengthSubset(frequentSet[i][j], k);
				print(s);
				for(l = 0; l < s.size(); ++l) {
					sComplement = getSubsetComplement(frequentSet[i][j], s[l]);
					// Getting Confidence values
					freq = getFrequency(s[l]);
					cout<<"Freq:" <<freq;
					supS = (double) freq/T;
					supS = supS*100;
					freq = frequentSetCount[i][j];
					supL = (double) freq/T;
					supL = supL*100;
					cout<<" SupS:" <<supS <<" SupL:" <<supL <<"\n";
					if(supS == 0 || supL == 0)
						break;
					conf = (supL/supS)*100;
					// Printing to File
					if(conf >= minConfidence) {
						for(p = 0; p < s[l].size(); ++p)
							outputFile<<s[l][p] <<" ";
						outputFile<<" ---> ";
						for(p = 0; p < sComplement.size(); ++p)
							outputFile<<sComplement[p] <<" ";
						outputFile<<"\t" <<conf <<"\n";
					}
					printVector(s[l]);
					cout<<"\t";
					printVector(sComplement);
					cout<<" " <<conf <<"\n";
					sComplement.clear();
				}
				s.clear();
			}
			cout<<"\n";
		}
	}
	outputFile.close();
}

// Copy confirmed frequent sets from L to frequentSet in ascending order of size of itemset
void copyToFrequentSet(vector<vector<int> > L, vector<int> Lfreq) {
	vector<vector<int> > tmp;
	vector<int> tmpFreq;
	int i, j, k, maxSize = 1;
	for(i = 1; i <= maxSize; ++i) {
		for(j = 0; j < L.size(); ++j) {
			if(L[j].size() > maxSize)
				maxSize = L[j].size();
			if(L[j].size() == i) {
				tmp.push_back(L[j]);
				tmpFreq.push_back(Lfreq[j]);
			}	
		}
		frequentSet.push_back(tmp);
		frequentSetCount.push_back(tmpFreq);
		tmp.clear();
		tmpFreq.clear();
	}
}

// Prints frequent itemsets and their support to a file
void printFrequentItemset(char* outputFile) {
	int i, j, k;
	fstream transactionFile;
	transactionFile.open(outputFile, ios::out);
	copyToFrequentSet(L, Lfreq);
	// Print frequentSet to file
	for(i = 0; i < frequentSet.size(); ++i) {
		for(j = 0; j < frequentSet[i].size(); ++j) {
			for(k = 0; k < frequentSet[i][j].size(); ++k) 
				transactionFile<<frequentSet[i][j][k] <<" ";
			transactionFile<<" = " <<frequentSetCount[i][j] <<"\n";
		}
	}
	transactionFile.close();
}

// Check if subset is in vector of itemsets(C, L, posFreqSet)
int isSubsetInItemsetVector(vector<int>subset, vector<vector<int> > itemsets) {
	int i, j, k, flag;
	for(i = 0; i < itemsets.size(); ++i) {
		if(subset.size() > itemsets[i].size())
			continue;
		flag = 1;
		k = 0;
		for(j = 0; j < itemsets[i].size(); ++j) {
			if(itemsets[i][j] != subset[k]) {
				flag = 0;
				break;
			}
			k++;
		}
		if(flag == 1)
			return i;
	}
	return -1;
}

FPnode* initialiseTree(FPnode* F, int item = -1, int freq = -1, FPnode* next = NULL, FPnode* parent = NULL) {
	F = new(FPnode);
	F->item = item;
	F->freq = freq;
	F->next = next;
	F->parent = parent;
	return F;
}

// Add items to FPTree
void addTransactionToTree(vector<int> items, FPnode* T, int pos = 0) {
	if(pos >= items.size())
		return;
	int i, item = items[pos];
	FPnode* I;
	// If item is already a child, increment freq
	if(T->children.size() > 0) {
		for(i = 0; i < T->children.size(); ++i) {
			if(T->children[i]->item == item) {
				T->children[i]->freq++;
				addTransactionToTree(items, T->children[i], pos+1);
				return;
			}
		}
	}
	// If item is not a child, create node and add it to children
	FPnode* Tmp = new(FPnode);
	Tmp = initialiseTree(Tmp, item, 1, NULL, T);
	T->children.push_back(Tmp);
	// Add new child to header
	if(Header[findItemInHeader(item)].next == NULL)
		Header[findItemInHeader(item)].next = Tmp;
	else {
		I = Header[findItemInHeader(item)].next;
		while(I->next != NULL)
			I = I->next;
		I->next = Tmp;
	}
	// Call function recursively to add next item
	addTransactionToTree(items, T->children[T->children.size()-1], pos+1);
}

void constructFPTree() {
	int i, j, k;
	for(i = 0; i < sortedTransactions.size(); ++i) {
		addTransactionToTree(sortedTransactions[i], F);
	}
}

// Add items to conditional FPTree
void addTransactionToCondTree(vector<int> items, condFPnode* T, int freq, int pos = 0) {
	if(pos >= items.size())
		return;
	int i, item = items[pos];
	condFPnode* I;
	// If item is already a child, increment freq
	if(T->children.size() > 0) {
		for(i = 0; i < T->children.size(); ++i) {
			if(T->children[i]->item == item) {
				T->children[i]->freq += freq;
				addTransactionToCondTree(items, T->children[i], freq, pos+1);
				return;
			}
		}
	}
	// If item is not a child, create node and add it to children
	condFPnode* Tmp = new(condFPnode);
	Tmp->item = item;
	Tmp->freq = freq;
	T->children.push_back(Tmp);
	// Call function recursively to add next item
	addTransactionToCondTree(items, T->children[T->children.size()-1], freq, pos+1);
}

condFPnode* constructCondFPTree(condFPnode* F, vector<vector<int> > itemset, vector<int> freq) {
	int i, j, k;
	for(i = 0; i < itemset.size(); ++i) {
		addTransactionToCondTree(itemset[i], F, freq[i]);
	}
	return F;
}

condFPnode* conditionalFPTree(header item_head) {
	FPnode* T;
	condFPnode* C;
	vector<vector<int> > condPattern;
	vector<int> condFreq;			
	vector<int> tmp;
	vector<FPnode*> Nlist;
	int freq, i;
	T = item_head.next;
	// Get Conditional Pattern Base
	while(T != NULL) {
		Nlist.push_back(T);
		T = T->next;
	}
	for(i = 0; i < Nlist.size(); ++i) {
		freq = Nlist[i]->freq;
		T = Nlist[i]->parent;
		tmp.clear();
		while(T->item != -1) {
			tmp.push_back(T->item);
			T = T->parent;
		}
		if(!tmp.empty()) {
			reverse(tmp.begin(), tmp.end());
			condPattern.push_back(tmp);
			condFreq.push_back(freq);
		}
	}
	// Construct conditional FPTree
	C = new(condFPnode);
	C->item = -1;
	C->freq = -1;
	C = constructCondFPTree(C, condPattern, condFreq);
	return C;
}

// Go through conditional FPTree using DFS and get every path. Then get all 2 or more length subsets of each path and put in frequentSet L
void DFSMiningCondTree(header item_head, condFPnode* C, vector<int> itemset) {
	double support = ((double) C->freq/T)*100;
	int i, j;
	vector<vector<int> > subsets;
	if(support >= minSupport)
		itemset.push_back(C->item);
	// If Leaf node
	if(C->children.size() == 0) {
		for(i = 1; i <= itemset.size()-1; ++i) {
			subsets = generateKlengthSubset(itemset, i);
			for(j = 0; j < subsets.size(); ++j) {
				subsets[j].push_back(item_head.item);
				if(isSubsetInItemsetVector(subsets[j], L) == -1)
					L.push_back(subsets[j]);
			}
		}
		itemset.push_back(item_head.item);
		if(isSubsetInItemsetVector(itemset, L) == -1)
			L.push_back(itemset);
		return;
	}
	for(i = 0; i < C->children.size(); ++i) {
		DFSMiningCondTree(item_head, C->children[i], itemset);
	}
}

// Mine frequent itemsets from conditional FPTree
void mineCondTree(header item_head, condFPnode* C) {
	vector<int> itemset;
	DFSMiningCondTree(item_head, C, itemset);
}

void fpgrowth(char* inputFile, char* outputFile) {
	int i, j, k, freq;
	condFPnode* C;
	getTransactions(inputFile);
	F = initialiseTree(F);
	constructFPTree();
	for(i = Header.size()-1; i >= 0; --i) {
		L.push_back(vector<int>(1, Header[i].item));
		if(i == 0)
			break;
		// Generate conditional FPTree with conditional pattern basea as transactions
		C = conditionalFPTree(Header[i]);
		// Use conditional FPtree to get frequent sets
		mineCondTree(Header[i], C);
	}
	// Get frequency of frequent itemsets and sort each itemset
	for(i = 0; i < L.size(); ++i) {
		sort(L[i].begin(), L[i].end());
		freq = getFrequency(L[i]);
		Lfreq.push_back(freq);
	}
	printFrequentItemset(outputFile);
}

int main() {
	char* inputFile = "sampleInput.txt";
	char* outputFile = "sampleOutput.txt";
	minSupport = 22.0;
	minConfidence = 50.0;
	fpgrowth(inputFile, outputFile);
	generateAssociationRules();
	return 0;
}

