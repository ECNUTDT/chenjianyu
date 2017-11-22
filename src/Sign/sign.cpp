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
compile:  g++ -o test sign.cpp -L. -lpbc -lgmp -ljson
execute:  ./test < ../../data/param/a.param
**/

const string configPath = "../../data/config/config";
const string PPPath = "../../data/setup_data/PP";
const string ExtractDataPath = "../../data/extract_data/data";

const int d = 9;

int N;
//N=21 n=20
int M;
element_t g;
element_t g1;
element_t g2;
element_t v;
element_t y;
static int qxIndex = 0;
static int txIndex = 0;

//
unsigned char* transfer(element_t t){

  int leng = element_length_in_bytes(t);
  unsigned char *p = new unsigned char[leng + 1];
  int writeLength = element_to_bytes_compressed(p,t);
  return p;
}

void getJsonValueNKey(Json::Value value,int tag,pairing_t pairing,element_t *arrN,element_t *arrM){
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
           if(tag == 1){
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
				if(type == "D"){
					element_init_G1(arrN[index-1],pairing);
					element_from_bytes_compressed(arrN[index-1],sb);
				}else if(type == "d"){
					element_init_G1(arrM[index-1],pairing);
					element_from_bytes_compressed(arrM[index-1],sb);
				}
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

void generateMLengthBitsString(int arr[]){
	srand((unsigned)time(NULL));
	for(int i = 0; i < M; ++i){
		arr[i] = rand() % 2;
		// arr[i] = 1;
	}
}

void calculateTheFirstPart(element_t v,element_t *arrM,int *message,element_t *randomSi,element_t *arrFirstPart,element_t *arrDi,pairing_t pairing){
	for(int i = 0; i < N-1; i++){
		element_t init;
		element_init_G1(init,pairing);
		element_set1(init);
		element_t mulResult;
		element_init_G1(mulResult,pairing);
		element_t powResult;
		element_init_G1(powResult,pairing);
		// element_printf("0---%B\n",init);
		for(int j = 0; j < M; j++){
			mpz_t mpz_temp;
			mpz_init_set_ui(mpz_temp,message[j]);
			element_t element_temp;
			element_init_Zr(element_temp,pairing);
			element_set_mpz(element_temp,mpz_temp);
			element_t powTemp;
			element_init_G1(powTemp,pairing);
			element_pow_zn(powTemp,arrM[j],element_temp);
			element_mul(mulResult,init,powTemp);
			element_set(init,mulResult);
		}
		element_t mulResult2;
		element_init_G1(mulResult2,pairing);
		element_mul(mulResult2,v,mulResult);
		element_pow_zn(powResult,mulResult2,randomSi[i]);
		element_t mulResult3;
		element_init_G1(mulResult3,pairing);
		element_mul(mulResult3,arrDi[i],powResult);
		element_set(arrFirstPart[i],mulResult3);
	}
}


int main(int argc, char **argv){

	//=============read data from file end===================
	string config = readFile(configPath);
	string PP = readFile(PPPath);
	string extractData = readFile(ExtractDataPath);

	pairing_t pairing;
	pbc_demo_pairing_init(pairing, argc, argv);

	Json::Reader reader;
	Json::Reader reader1;
	Json::Reader reader2;

	Json::Value value;
	Json::Value value1;
	Json::Value value2;

	if(!reader2.parse(config,value2)){
		cout << "cannot read config file!\n";
		return 0;
	}
	getJsonValueNKey(value2,2,pairing,NULL,NULL);


	element_t *arrDi = new element_t[N-1];
	element_t *arrdi = new element_t[N-1];
	if(!reader.parse(extractData,value1)){
		cout << "cannot read extractData file!\n";
		return 0;
	}
	getJsonValueNKey(value1,1,pairing,arrDi,arrdi); 


	element_t *arrN = new element_t[N];
	element_t *arrM = new element_t[M];
	if(!reader.parse(PP,value)){
		cout << "cannot read PP file!\n";
		return 0;
	}
	getJsonValueNKey(value,0,pairing,arrN,arrM); 
	//=============read data from file end===================

	//=============generate message===================
	int *message = new int[M];
	generateMLengthBitsString(message);
	// for (int i = 0; i < M; ++i){
	// 	printf("%d",message[i]);
	// }
	// printf("\n");
	//=============generate message===================

	//=============generate random===================
	element_t *randomSi = new element_t[N-1];
	for(int i = 0; i < N - 1; ++i){
		element_init_Zr(randomSi[i],pairing);
		element_random(randomSi[i]);
	}
	//=============generate random===================

	//=============calculate the first part==================
	element_t *arrFirstPart = new element_t[N-1];
	for(int i = 0; i < N - 1; ++i){
		element_init_G1(arrFirstPart[i],pairing);
	}
	calculateTheFirstPart(v,arrM,message,randomSi,arrFirstPart,arrDi,pairing);
	//arrFirstPart[]

	//=============calculate the first part===================


	//=============calculate the third part===================
	element_t *arrThirdPart = new element_t[N-1];
	for(int i = 0; i < N - 1; ++i){
		element_init_G1(arrThirdPart[i],pairing);
		element_t tempRandomSi;
		element_init_Zr(tempRandomSi,pairing);
		element_neg(tempRandomSi,randomSi[i]);
		element_pow_zn(arrThirdPart[i],g,tempRandomSi);
	}
	//arrThirdPart[]

	//=============calculate the third part===================


	//=============write to file===================
	Json::Value root;  
	Json::StyledWriter sw;  
	ofstream os;
	os.open("../../data/sign_data/data");  
	for(int i = 0; i < N-1; i++){
		string index = "";
   	 	stringstream st;
   	 	st << (i+1);
    	st >> index;
    	index ="first-" + index;
    	root[index] = Json::Value((char*)transfer(arrFirstPart[i]));  	
	}
	for(int i = 0; i < N-1; i++){
		string index = "";
   	 	stringstream st;
   	 	st << (i+1);
    	st >> index;
    	index ="second-" + index;
    	root[index] = Json::Value((char*)transfer(arrdi[i]));  	
	}
	for(int i = 0; i < N-1; i++){
		string index = "";
   	 	stringstream st;
   	 	st << (i+1);
    	st >> index;
    	index ="third-" + index;
    	root[index] = Json::Value((char*)transfer(arrThirdPart[i]));  	
	}

	os << sw.write(root);  
	os.close();  
	//=============write to file===================



	delete []arrThirdPart;
	delete []arrFirstPart;
	delete []randomSi;
	delete []message;
	delete []arrdi;
	delete []arrDi;
	delete []arrN;
	delete []arrM;
	pairing_clear(pairing);           
}