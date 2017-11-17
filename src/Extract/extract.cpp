#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <iterator>
#include <string.h>
#include "/usr/local/include/pbc/pbc.h"
#include "/usr/include/json/json.h"
#include "/usr/local/include/pbc/pbc_test.h"
using namespace std;


/**
compile:  g++ -o test extract.cpp -L. -lpbc -lgmp -ljson
execute:  ./test < ../../data/param/a.param
**/

const string configPath = "../../data/config/config";
const string PPPath = "../../data/setup_data/PP";
const string MKPath = "../../data/setup_data/MK";

int N;
int M;


unsigned char *stringToChar(string str){
    char *p;
    int len = str.length();
    p=new char[len+1];
    strcpy(p,str.c_str());
   	unsigned char * sb = (unsigned char *)p;
	return sb;  
}

void getJsonValueNKey(Json::Value value,map<string,element_t> &map,int tag,pairing_t pairing){
	Json::Value::Members members;
	members = value.getMemberNames();
	int i = 0;
	for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
            string strKey = *iterMember;  
        	if(tag == 0){ 
        		element_t temp;
				element_init_G1(temp,pairing);
        		string str = value[strKey.c_str()].asString(); 
			    char *p;
			    int len = str.length();
			    p=new char[len+1];
			    strcpy(p,str.c_str());
			    unsigned char * sb = (unsigned char *)p;
				element_from_bytes_compressed(temp,sb);
				map.insert(pair<string,element_t>(strKey,temp));
        	}
           if(tag == 1){
           		element_t t;
            	element_init_Zr(t,pairing);
				string strVal = value[strKey.c_str()].asString();  
				mpz_t gmp_y;
				const char * chargmp_y = strVal.c_str();
				mpz_init_set_str(gmp_y,chargmp_y,10);
				element_set_mpz(t,gmp_y);
				map.insert(pair<string,element_t>(strKey,t));
            }
            if(tag == 2){
                int iVal = value[strKey.c_str()].asInt(); 
                if (strKey == "N"){
                	N = iVal;
                }
                else if(strKey == "M"){
                	M = iVal;
                }
            }  
     }  
}

string readFile(string path){
	stringstream ss;
	fstream readPathData(path.c_str());
	ss << readPathData.rdbuf();
	string data = ss.str();
	ss.clear();
	ss.str("");
	readPathData.close();
	return data;
}


int main(int argc, char **argv){

	string config = readFile(configPath);
	string PP = readFile(PPPath);
	string MK = readFile(MKPath);

	map<string,element_t> PPMap;
	map<string,element_t> MKMap;
	map<string,element_t> configMap;

	pairing_t pairing;
	pbc_demo_pairing_init(pairing, argc, argv);

	Json::Reader reader;
	Json::Reader reader1;
	Json::Reader reader2;

	Json::Value value;
	Json::Value value1;
	Json::Value value2;


	// ======================
	if(!reader.parse(PP,value)){
		cout << "cannot read PP file!\n";
		return 0;
	}
	 getJsonValueNKey(value,PPMap,0,pairing); 


	if(!reader1.parse(MK,value1)){
		cout << "cannot read MK file!\n";
		return 0;
	}
	getJsonValueNKey(value1,MKMap,1,pairing);
	


	if(!reader2.parse(config,value2)){
		cout << "cannot read config file!\n";
		return 0;
	}
	getJsonValueNKey(value2,configMap,2,pairing);
	// ======================

	map<string,element_t>::iterator iter;
	for(iter = PPMap.begin(); iter != PPMap.end();iter++){
		cout << iter->first;
		element_printf(":%B\n",iter->second);
	}
	for(iter = MKMap.begin(); iter != MKMap.end();iter++){
		cout << iter->first;
		element_printf(":%B\n",iter->second);
	}
	printf("N=%d,M=%d\n",N,M);

	pairing_clear(pairing);           
}