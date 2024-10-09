#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <random>
using namespace std;

struct Range
{
    uint32_t low;
    uint32_t high;
};

struct DimInfo {
    bool isCut = false;
    uint32_t numCuts = 0;
    Range range;
};

struct Rule
{
    array<Range, 5> fields;
    uint32_t priority;
    bool partiallyOverlaps(Rule rule, int dim) {
        if (fields[dim].low > rule.fields[dim].high || fields[dim].high < rule.fields[dim].low) {
            return false;
        }
        return true;
    }
    bool isInsideHyperCube(array<Range, 5> hyperCube) {
        for (int i = 0; i < 5; i++) {
            if (fields[i].low < hyperCube[i].low || fields[i].high > hyperCube[i].high) {
                return false;
            }
        }
        return true;
    }
    bool intersectsHyperCube(array<Range, 5> hyperCube) {
        for (int i = 0; i < 5; i++) {
            if (fields[i].low > hyperCube[i].high || fields[i].high < hyperCube[i].low) {
                return false;
            }
        }
        return true;
    }
};

class HyperCutsNode {
    public:
        bool isLeaf;
        vector<uint32_t>* ruleIDs;
        array<DimInfo, 5> dimInfo;
        vector <HyperCutsNode*> children;
        HyperCutsNode(vector<uint32_t>*& ruleIDs, array<DimInfo, 5> dimInfo);
        HyperCutsNode(vector<uint32_t>*& ruleIDs);
        ~HyperCutsNode();
        void printTree();
    
};
HyperCutsNode:: HyperCutsNode(vector<uint32_t>*& ruleIDs, array<DimInfo, 5> dimInfo) {
    this->ruleIDs = ruleIDs;
    this->dimInfo = dimInfo;
    //isLeaf = false;
}    
HyperCutsNode::HyperCutsNode(vector<uint32_t>*& ruleIDs) {
    this->ruleIDs = ruleIDs;
    //isLeaf = true;
}
HyperCutsNode::~HyperCutsNode() {
    if (ruleIDs != nullptr) {
        delete ruleIDs;
    }
    for (int i = 0; i < children.size(); i++) {
        delete children[i];
    }
}
void HyperCutsNode::printTree() {
    if (isLeaf) {
        cout << "Leaf node" << endl;
        cout << "Rule IDs: ";
        for (int i = 0; i < ruleIDs->size(); i++) {
            cout << ruleIDs->at(i) << " ";
        }
        cout << endl;
    }
    else {
        cout << "Non-leaf node" << endl;
        cout << "Number of children: " << children.size() << endl;
        for (int i = 0; i < children.size(); i++) {
            children[i]->printTree();
        }
    }
}
typedef array<uint32_t, 5> Packet;
typedef vector<Packet> Packets;
typedef vector<Rule> Rules;

Packets packets;
Rules rules;
uint8_t binth_val;
double spfac_val;

uint32_t determineUniqueElementsInDim(vector<uint32_t>*& ruleIDs, int dim);
void readRulesFromFile(string fileName);
void readPacketsFromFile(string fileName);
void createHyperCutsTree(HyperCutsNode*& node, vector<uint32_t>& ruleIDs, array<DimInfo, 5>& dimInfo);
void createHyperCutsTreeRecursive(HyperCutsNode*& node, vector<uint32_t>*& ruleIDs);


void generateRandomCuts(std::mt19937& gen, std::uniform_int_distribution<>& dis, int threshold, array<DimInfo,5>& dimInfo) {
    cout << "Inside generateRandomCuts: Threshold: " << threshold << endl;
    bool found = false;
    uint32_t totalCuts = 1;
    //uint32_t numCut;
    while (!found) {
        for (int i = 0; i < 5; i++) {
            if (!dimInfo[i].isCut) {
                dimInfo[i].numCuts = 0;
                continue;
            }
            //numCut = dis(gen);
            dimInfo[i].numCuts = dis(gen);
            totalCuts *= dimInfo[i].numCuts;
        }
        if (totalCuts < threshold) {
            found = true;
        }
        totalCuts = 1; //Resetting the totalCuts
    }
}


int main (int argc, char** argv) {
    if (argc != 5) {
        cout << "Usage: " << argv[0] << " <rules file> <packets file> <spfac value> <binth value>" << endl;
        return 1;
    }

    readRulesFromFile(argv[1]);
    readPacketsFromFile(argv[2]);
    cout << "Number of rules: " << rules.size() << endl;
    cout << "Number of packets: " << packets.size() << endl;
    spfac_val = atol(argv[3]);
    binth_val = atoi(argv[4]);
    HyperCutsNode* root = nullptr;

    vector<uint32_t>* ruleIDs = new vector<uint32_t>;
    for (int i = 0; i < rules.size(); i++) {
        ruleIDs->push_back(i);
    }

    createHyperCutsTreeRecursive(root, ruleIDs);


    cout << "Printing the tree" << endl;
    if (root != nullptr) {
        root->printTree();
    }
    else {
        cout << "Root is null" << endl;
    }
    return 0;
}


void createHyperCutsTree(HyperCutsNode*& node, vector<uint32_t>*& ruleIDs, array<DimInfo, 5>& dimInfo) {
    if (ruleIDs->size() <= binth_val) {
        node = new HyperCutsNode(ruleIDs, dimInfo);
        return;
    }
    //Determining the total number of cuts to be made with the ruleIDs 
}

bool areRuleIDsSame(vector<uint32_t>* ruleIDs1, vector<uint32_t>* ruleIDs2) {
    if (ruleIDs1->size() != ruleIDs2->size()) {
        return false;
    }
    for (int i = 0; i < ruleIDs1->size(); i++) {
        if (ruleIDs1->at(i) != ruleIDs2->at(i)) {
            return false;
        }
    }
    return true;
}

void createHyperCutsTreeRecursive(HyperCutsNode*& node, vector<uint32_t>*& ruleIDs) {
    if (ruleIDs->size() <= binth_val) {
        node = new HyperCutsNode(ruleIDs);
        return;
    }

    //determining the number of non overlapping ranges in each dimension
    vector<uint32_t> numUniqueVals;
    for (int i = 0; i < 5; i++) {
        uint32_t numUniqueVal = determineUniqueElementsInDim(ruleIDs, i);
        numUniqueVals.push_back(numUniqueVal);
        //cout << "Number of unique values in dimension " << i << ": " << numUniqueVal << endl;
    }

    //determining the average number of non overlapping ranges in all dimensions
    double avgNumUniqueVals = 0;
    for (int i = 0; i < 5; i++) {
        avgNumUniqueVals += numUniqueVals[i];
    }
    avgNumUniqueVals /= 5;

    //determining the fields with overlapping ranges greater than the average
    array<DimInfo, 5> dimInfo;
    uint32_t numDimsCut = 0;
    for (int i = 0; i < 5; i++) {
        if (numUniqueVals[i] >= avgNumUniqueVals) {
            dimInfo[i].isCut = true;
            numDimsCut++;
        }
        else {
            dimInfo[i].isCut = false;
        }
    }

    //Once the fields are determined, determining the number of cuts to be done in each field

    //Determining the cut threshold
    uint32_t cutThreshold = spfac_val * sqrt(ruleIDs->size()); //The threshold of the number of cuts to be made in each node.
    static std::random_device rd;
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    static std::uniform_int_distribution<> dis(0, ceil(3*cutThreshold/numDimsCut)); //Uniform distribution between 0 and size for selecting a random dimension to be cut.


    //Determining the number of cuts on each field so that the product of all cuts is less than cutThreshold
    generateRandomCuts(gen, dis, cutThreshold, dimInfo);

    //Printing dimInfo
    cout << "Printing dimInfo" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "Field " << i << ": ";
        if (dimInfo[i].isCut) {
            cout << "Cut ";
        }
        else {
            cout << "No cut ";
        }
        cout << "Number of cuts: " << dimInfo[i].numCuts << endl;
    }

    //Determine the ranges for each field
    
    for (int i = 0; i < 5; i++) {
        Range r;
        r.low = rules[0].fields[i].low;
        r.high = rules[0].fields[i].high;
        for (int j = 1; j < rules.size(); j++) {
            if (rules[j].fields[i].low < r.low) {
                r.low = rules[j].fields[i].low;
            }
            if (rules[j].fields[i].high > r.high) {
                r.high = rules[j].fields[i].high;
            }
        }
        dimInfo[i].range = r;
    }

    // //Printing the ranges for each field
    // for (int i = 0; i < 5; i++) {
    //     cout << "Range for field " << i << ": " << dimInfo[i].range.low << " - " << dimInfo[i].range.high << endl;
    // }
    // //printing the number of cuts on each field
    // for (int i = 0; i < 5; i++) {
    //     cout << "Number of cuts on field " << i << ": " << dimInfo[i].numCuts << endl;
    // }
    //Creating the non leaf node.
    node = new HyperCutsNode(ruleIDs, dimInfo);
    //Determining the children of the node based on dimInfo
    // uint32_t numChildren = 1;
    // for (int i = 0; i < 5; i++) {
    //     if (dimInfo[i].isCut) {
    //         numChildren *= (dimInfo[i].numCuts + 1);
    //     }
    // }
    //Determining the cubes for each child
    
    vector<array<Range, 5>> hyperCubes;
    for (int i = 0; i < dimInfo[0].numCuts + 1; i++){
        for (int j = 0; j < dimInfo[1].numCuts + 1; j++) {
            for (int k = 0; k < dimInfo[2].numCuts + 1; k++) {
                for (int l = 0; l < dimInfo[3].numCuts +1 ; l++) {
                    for (int m = 0; m < dimInfo[4].numCuts + 1; m++) {
                        vector<uint32_t>* childRuleIDs = new vector<uint32_t>;
                        //Define the hypercube for the child
                        array<Range, 5> hyperCube;
                        hyperCube[0].low = dimInfo[0].range.low + i*((float)(dimInfo[0].range.high - dimInfo[0].range.low)/(dimInfo[0].numCuts + 1));
                        hyperCube[0].high = hyperCube[0].low + ((float)(dimInfo[0].range.high - dimInfo[0].range.low)/(dimInfo[0].numCuts + 1));
                        hyperCube[1].low = dimInfo[1].range.low + j*((float)(dimInfo[1].range.high - dimInfo[1].range.low)/(dimInfo[1].numCuts + 1));
                        hyperCube[1].high = hyperCube[1].low + ((float)(dimInfo[1].range.high - dimInfo[1].range.low)/(dimInfo[1].numCuts + 1));
                        hyperCube[2].low = (dimInfo[2].range.low) + k*((float)(dimInfo[2].range.high - dimInfo[2].range.low)/(dimInfo[2].numCuts + 1));
                        hyperCube[2].high = hyperCube[2].low + ((float)(dimInfo[2].range.high - dimInfo[2].range.low)/(dimInfo[2].numCuts + 1));
                        hyperCube[3].low = (dimInfo[3].range.low) + l*((float)(dimInfo[3].range.high - dimInfo[3].range.low)/(dimInfo[3].numCuts + 1));
                        hyperCube[3].high = hyperCube[3].low + ((float)(dimInfo[3].range.high - dimInfo[3].range.low)/(dimInfo[3].numCuts + 1));
                        hyperCube[4].low = (dimInfo[4].range.low) + m*((float)(dimInfo[4].range.high - dimInfo[4].range.low)/(dimInfo[4].numCuts + 1));
                        hyperCube[4].high = hyperCube[4].low + ((float)(dimInfo[4].range.high - dimInfo[4].range.low)/(dimInfo[4].numCuts + 1));
                        for (int n = 0; n < ruleIDs->size(); n++) {
                            if (rules[(*ruleIDs)[n]].intersectsHyperCube(hyperCube)) {
                                childRuleIDs->push_back((*ruleIDs)[n]);
                            }
                        }
                        if (areRuleIDsSame(childRuleIDs, ruleIDs)) {
                            cout << "Child rule IDs are same as parent rule IDs" << endl;
                            //No need to create a child node.
                            delete childRuleIDs;
                            continue;
                        }
                        HyperCutsNode* childNode = nullptr;
                        node->children.push_back(childNode);
                        createHyperCutsTreeRecursive(childNode, childRuleIDs);
                    }
                }
            }
        }
    }
    return;
}
    
uint32_t determineUniqueElementsInDim(vector<uint32_t>*& ruleIDs, int dim) {
    vector<uint32_t> uniqueVals;
    for (int i = 0; i < ruleIDs->size(); i++) {
        uniqueVals.push_back(rules[(*ruleIDs)[i]].fields[dim].low);
        uniqueVals.push_back(rules[(*ruleIDs)[i]].fields[dim].high);
    }
    sort(uniqueVals.begin(), uniqueVals.end());
    uniqueVals.erase(unique(uniqueVals.begin(), uniqueVals.end()), uniqueVals.end());
    return uniqueVals.size();
}

void readPacketsFromFile(string fileName) {
	ifstream inputFile(fileName);
	//vector<Packet> packets;
	uint32_t sip, dip, sp, dp, proto;
	//int count = 0;
	while(inputFile >> sip >> dip >> sp >> dp >> proto) {
		Packet packet;
		packet[0] = sip;
		packet[1] = (dip);
		packet[2] = (sp);
		packet[3] = (dp);
		packet[4] = (proto);
		packets.push_back(packet);
	}
    inputFile.close();
}

void readRulesFromFile(string fileName) {
    static uint32_t ruleID = 0;
    ifstream inputFile(fileName);
    //vector<Rule> rules;
    uint32_t sip_low, sip_high, dip_low, dip_high, sp_low, sp_high, dp_low, dp_high, proto_low, proto_high;
    //int count = 0;
    while(inputFile >> sip_low >> sip_high >> dip_low >> dip_high >> sp_low >> sp_high >> dp_low >> dp_high >> proto_low >> proto_high) {
        Rule rule;
        rule.fields[0].low = sip_low;
        rule.fields[0].high = sip_high;
        rule.fields[1].low = dip_low;
        rule.fields[1].high = dip_high;
        rule.fields[2].low = sp_low;
        rule.fields[2].high = sp_high;
        rule.fields[3].low = dp_low;
        rule.fields[3].high = dp_high;
        rule.fields[4].low = proto_low;
        rule.fields[4].high = proto_high;
        rule.priority = ruleID++;
        rules.push_back(rule);
    }
    inputFile.close();
}