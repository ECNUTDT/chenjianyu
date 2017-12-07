#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string.h>
#include "/usr/local/include/pbc/pbc.h"
#include "/usr/include/json/json.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include <time.h>  
#include <stdlib.h>

using namespace std;


/**
compile:  g++ -o test verify.cpp -L. -lpbc -lgmp -ljson
execute:  ./test < ../../data/param/a.param
**/

const string configPath = "../../data/config/config";
const string PPPath = "../../data/setup_data/PP";
const string ExtractDataPath = "../../data/extract_data/data";
const string OmigaDataPath = "../../data/extract_data/omiga";
const string QxDataPath = "../../data/extract_data/Qx";
const string SignDataPath = "../../data/sign_data/data";

const int d = 9;
static int txIndex = 0;
int N;
//N=21 n=20
int M;
element_t g;
element_t g1;
element_t g2;
element_t v;

//
unsigned char* transfer(element_t t){

  int leng = element_length_in_bytes(t);
  unsigned char *p = new unsigned char[leng + 1];
  int writeLength = element_to_bytes_compressed(p,t);
  return p;
}

void getConfig(Json::Value value,pairing_t pairing){
	Json::Value::Members members;
	members = value.getMemberNames();
	int i = 0;
	for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
            string strKey = *iterMember;  
            int iVal = value[strKey.c_str()].asInt(); 
            if (strKey == "N"){
            	N = iVal;
            }
            else if(strKey == "M"){
            	M = iVal;
            }
    }
}

void getJsonValueNKey(Json::Value value,pairing_t pairing,element_t *firstPart,element_t *secondPart,element_t *thirdPart){
	Json::Value::Members members;
	members = value.getMemberNames();
	int i = 0;
	for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
            string strKey = *iterMember;  
    		element_t temp;
			element_init_G1(temp,pairing);
    		string str = value[strKey.c_str()].asString(); 
		    char *p;
		    int len = str.length();
		    p=new char[len+1];
		    strcpy(p,str.c_str());
		    unsigned char * sb = (unsigned char *)p;
			element_from_bytes_compressed(temp,sb);
			int pos = strKey.find("-");
			string type = strKey.substr(0,pos);
			string num = strKey.substr(pos+1,strKey.size());
			stringstream ss;
			ss << num;
			int index;
			ss >> index;
			if(type == "first"){
				element_init_G1(firstPart[index-1],pairing);
				element_from_bytes_compressed(firstPart[index-1],sb);
			}
			else if(type == "second"){
				element_init_G1(secondPart[index-1],pairing);
				element_from_bytes_compressed(secondPart[index-1],sb);
			}
			else if(type == "third"){
				element_init_G1(thirdPart[index-1],pairing);
				element_from_bytes_compressed(thirdPart[index-1],sb);
			}   	
     }  
}

void getPPValue(Json::Value value,pairing_t pairing,element_t *arrM,element_t *arrN){
	Json::Value::Members members;
	members = value.getMemberNames();
	for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
    	string strKey = *iterMember;
        element_t temp;
		element_init_G1(temp,pairing);
		string str = value[strKey.c_str()].asString(); 
	    char *p;
	    int len = str.length();
	    p=new char[len+1];
	    strcpy(p,str.c_str());
	    unsigned char * sb = (unsigned char *)p;
		element_from_bytes_compressed(temp,sb);
		int pos = strKey.find("-");
		if(pos == -1){
			if(strKey == "g"){
				element_init_G1(g,pairing);
				element_from_bytes_compressed(g,sb);
			}
			else if(strKey == "g1"){
				element_init_G1(g1,pairing);
				element_from_bytes_compressed(g1,sb);
			}
			else if(strKey == "g2"){
				element_init_G1(g2,pairing);
				element_from_bytes_compressed(g2,sb);
			}
			else if(strKey == "v"){
				element_init_G1(v,pairing);
				element_from_bytes_compressed(v,sb);
			}
		}
		else{
			string type = strKey.substr(0,pos);
			string num = strKey.substr(pos+1,strKey.size());
			stringstream ss;
			ss << num;
			int index;
			ss >> index;
			if(type == "t"){
				element_init_G1(arrN[index-1],pairing);
				element_from_bytes_compressed(arrN[index-1],sb);
			}else if(type == "v"){
				element_init_G1(arrM[index-1],pairing);
				element_from_bytes_compressed(arrM[index-1],sb);
			}
		}
    }
}

void getExtractDataValue(Json::Value value,pairing_t pairing,element_t *arrTx){
	Json::Value::Members members;
	members = value.getMemberNames();
	for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
    	string strKey = *iterMember;
    	int pos = strKey.find("-");
    	string type = strKey.substr(0,pos);
    	if(type == "T"){
    		element_t temp;
			element_init_G1(temp,pairing);
			string str = value[strKey.c_str()].asString(); 
		    char *p;
		    int len = str.length();
		    p=new char[len+1];
		    strcpy(p,str.c_str());
		    unsigned char * sb = (unsigned char *)p;
			element_from_bytes_compressed(temp,sb);
			string num = strKey.substr(pos+1,strKey.size());
			stringstream ss;
			ss << num;
			int index;
			ss >> index;
			element_init_G1(arrTx[index-1],pairing);
			element_from_bytes_compressed(arrTx[index-1],sb);
    	}
    }
}

void calculatePartThreePartTwo(int *message,pairing_t pairing,element_t *arrM,element_t fixedValue){
	element_t powTemp;
	element_init_G1(powTemp,pairing);
	element_t mulTemp;
	element_init_G1(mulTemp,pairing);
	element_t initValue;
	element_init_G1(initValue,pairing);
	element_set(initValue,v);
	for(int i = 0; i < M; ++i){
		mpz_t exp;
		mpz_init_set_ui(exp,message[i]);
		element_pow_mpz(powTemp,arrM[i],exp);
		element_mul(mulTemp,initValue,powTemp);
	}
	// fixedValue only has two result depending on the num of 1 or 0
	element_set(fixedValue,mulTemp);
}

void getOmigaValue(Json::Value value,pairing_t pairing,element_t *omiga){
	Json::Value::Members members;
	members = value.getMemberNames();
	for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
    	string strKey = *iterMember;
    	int pos = strKey.find("-");
    	string type = strKey.substr(0,pos);
    	if(type == "O"){
			string str = value[strKey.c_str()].asString(); 
			string num = strKey.substr(pos+1,strKey.size());
			stringstream ss;
			ss << num;
			int index;
			ss >> index;
			element_init_Zr(omiga[index-1],pairing);
			mpz_t mpz_element;
			mpz_init_set_str(mpz_element,str.c_str(),10);   
			element_set_mpz(omiga[index-1],mpz_element); 	
    	}
    }
}

void getQxValue(Json::Value value,pairing_t pairing,element_t *arrQx){
	Json::Value::Members members;
	members = value.getMemberNames();
	for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
    	string strKey = *iterMember;
    	int pos = strKey.find("-");
    	string type = strKey.substr(0,pos);
    	if(type == "Q"){
			string str = value[strKey.c_str()].asString(); 
			string num = strKey.substr(pos+1,strKey.size());
			stringstream ss;
			ss << num;
			int index;
			ss >> index;
			element_init_Zr(arrQx[index-1],pairing);
			mpz_t mpz_element;
			mpz_init_set_str(mpz_element,str.c_str(),10);   
			element_set_mpz(arrQx[index-1],mpz_element); 	
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

	string configData = readFile(configPath);
	string signData = readFile(SignDataPath);
	string PPData = readFile(PPPath);
	string extractData = readFile(ExtractDataPath);
	string omigaData = readFile(OmigaDataPath);
	string qxData = readFile(QxDataPath);

	pairing_t pairing;
	pbc_demo_pairing_init(pairing, argc, argv);

	Json::Reader reader;
	Json::Reader reader1;
	Json::Reader reader2;
	Json::Reader reader3;
	Json::Reader reader4;
	Json::Reader reader5;


	Json::Value value;
	Json::Value value1;
	Json::Value value2;
	Json::Value value3;
	Json::Value value4;
	Json::Value value5;


	//=============read data from file end===================
	if(!reader.parse(configData,value)){
		cout << "cannot read configData file!\n";
		return 0;
	}
	getConfig(value,pairing);


	element_t *firstPart = new element_t[N-1];
	element_t *secondPart = new element_t[N-1];
	element_t *thirdPart = new element_t[N-1];	
	if(!reader1.parse(signData,value1)){
		cout << "cannot read signData file!\n";
		return 0;
	}
	getJsonValueNKey(value1,pairing,firstPart,secondPart,thirdPart);

	element_t *arrN = new element_t[N];
	element_t *arrM = new element_t[M];
	if(!reader2.parse(PPData,value2)){
		cout << "cannot read PPData file!\n";
		return 0;
	}
	getPPValue(value2,pairing,arrM,arrN);

	element_t *arrTx = new element_t[N];
	if(!reader3.parse(extractData,value3)){
		cout << "cannot read ExtractDataPath file!\n";
		return 0;
	}
	getExtractDataValue(value3,pairing,arrTx);

	element_t *omiga = new element_t[N-1];
	if(!reader4.parse(omigaData,value4)){
		cout << "cannot read omigaData file!\n";
		return 0;
	}
	getOmigaValue(value4,pairing,omiga);


	element_t *arrQx = new element_t[N-1];
	if(!reader5.parse(qxData,value5)){
		cout << "cannot read qxData file!\n";
		return 0;
	}
	getQxValue(value5,pairing,arrQx);

	//=============read data from file end===================


	//=============generate message===================
	int *message = new int[M];
	srand((unsigned)time(NULL));
	for(int i = 0; i < M; ++i){
		message[i] = rand() % 2;
	}
	//=============generate message===================


	//=============get PartThreePartTwo===================
	element_t fixedValue;
	element_init_G1(fixedValue,pairing);
	calculatePartThreePartTwo(message,pairing,arrM,fixedValue);
	//=============get PartThreePartTwo===================


	// //=============first===================
	// element_printf("-----------------------------\n");
	// element_t outMul;
	// element_init_GT(outMul,pairing);
	// element_set1(outMul);
	// for(int i = 0; i < N-1; ++i){
	// 	element_t result;
	// 	element_init_G1(result,pairing);

	// 	element_t e1;
	// 	element_init_GT(e1,pairing);
	// 	element_t e2;
	// 	element_init_GT(e2,pairing);
	// 	element_t e3;
	// 	element_init_GT(e3,pairing);
	// 	pairing_apply(e1,firstPart[i],g,pairing);
	// 	pairing_apply(e2,secondPart[i],arrTx[i],pairing);
	// 	pairing_apply(e3,thirdPart[i],fixedValue,pairing);

	// 	element_t mulTemp1;
	// 	element_init_GT(mulTemp1,pairing);
	// 	element_mul(mulTemp1,e1,e2);
	// 	element_mul(mulTemp1,mulTemp1,e3);


	// 	//calculate delta i,s(0)
	// 	element_t innerMul;
	// 	element_init_Zr(innerMul,pairing);
	// 	element_set1(innerMul);
	// 	for(int j = 0; j < N - 1; ++j){
	// 		if(i == j){
	// 			continue;
	// 		}
	// 		element_t numerator;
	// 		element_t denominator;
	// 		element_t quotient;
	// 		element_init_Zr(numerator,pairing);
	// 		element_init_Zr(denominator,pairing);
	// 		element_init_Zr(quotient,pairing);
	// 		element_neg(numerator,omiga[j]);
	// 		element_sub(denominator,omiga[i],omiga[j]);
	// 		element_div(quotient,numerator,denominator);
	// 		element_mul(innerMul,innerMul,quotient);
	// 	}	
	// 	element_t innerPow;
	// 	element_init_GT(innerPow,pairing);
	// 	element_pow_zn(innerPow,mulTemp1,innerMul);
	// 	element_mul(outMul,outMul,innerPow);
	// }
	// element_printf("outMul = %B\n",outMul);

	// element_t A;
	// element_init_GT(A,pairing);
	// pairing_apply(A,g1,g2,pairing);

	// element_printf("A = %B\n",A);

	// element_printf("-----------------------------\n");
	// //=============first===================





	//=============second===================
	element_t outMul;
	element_init_GT(outMul,pairing);
	element_set1(outMul);
	//change N-2
	for (int i = 0; i < N - 2; ++i){
		element_t innerMul;
		element_init_Zr(innerMul,pairing);
		element_set1(innerMul);
		//change N-2
		for(int j = 0; j < N - 2; ++j){
			if(i == j){
				continue;
			}
			element_t numerator;
			element_t denominator;
			element_t quotient;
			element_init_Zr(numerator,pairing);
			element_init_Zr(denominator,pairing);
			element_init_Zr(quotient,pairing);
			element_neg(numerator,omiga[j]);
			element_sub(denominator,omiga[i],omiga[j]);
			element_div(quotient,numerator,denominator);
			element_mul(innerMul,innerMul,quotient);
		}	

		element_t GTTemp;
		element_init_GT(GTTemp,pairing);
		element_t g2qi;
		element_init_G1(g2qi,pairing);
		element_pow_zn(g2qi,g2,arrQx[i]);
		pairing_apply(GTTemp,g2qi,g,pairing);
		element_t vari;
		element_init_GT(vari,pairing);
		element_pow_zn(vari,GTTemp,innerMul);
		element_mul(outMul,outMul,vari);
	}

	element_printf("outMul = %B\n",outMul);
	element_t A;
	element_init_GT(A,pairing);
	pairing_apply(A,g1,g2,pairing);
	element_printf("A = %B\n",A);

	//=============second===================


	//useless
	delete []arrN;
	
	delete []arrQx; 
	delete []omiga;
	delete []arrM;
	delete []firstPart;
	delete []secondPart;
	delete []thirdPart;
	element_clear(g);
	element_clear(g1);
	element_clear(g2);
	element_clear(v);
	pairing_clear(pairing);           
}