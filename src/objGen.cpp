#include <fstream> 
#include <vector> 
#include <memory.h>
#include <math.h>
#include "triangle.h"
#include "objGen.h"


using namespace std;

void SaveOBJ(std::vector <TriangleC> *v, std::string filename) {

	ofstream myfile;
	myfile.open(filename);

	myfile << "# Generated by Bedrich Benes bbenes@purdue.edu\n";
	myfile << "# vertices\n";
	for (unsigned int i = 0; i < v->size(); i++) {
		myfile << "v " << v->at(i).a.z << " " << v->at(i).a.y << " " << v->at(i).a.x << "\n";
		myfile << "v " << v->at(i).b.z << " " << v->at(i).b.y << " " << v->at(i).b.x << "\n";
		myfile << "v " << v->at(i).c.z << " " << v->at(i).c.y << " " << v->at(i).c.x << "\n";
	}
	int counter = 1;
	myfile << "# faces\n";
	for (unsigned int i = 0; i < v->size(); i++) {
		myfile << "f " << counter++;
		myfile << " " << counter++;
		myfile << " " << counter++ << " " << "\n";
	}
	myfile.close();

}


