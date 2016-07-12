#include "testClass.h"

using namespace std;

testClass::testClass(string fname, int a){
	name = fname;
	age = a;
}
void testClass::dummy(){
	cout << name << " " << age << endl;
}