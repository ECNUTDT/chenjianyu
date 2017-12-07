#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
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
            	element_init_Zr(y,pairing);
				string strVal = value[strKey.c_str()].asString();  
				mpz_t gmp_y;
				const char * chargmp_y = strVal.c_str();
				mpz_init_set_str(gmp_y,chargmp_y,10);
				element_set_mpz(y,gmp_y);
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

//q(i)
void calculateQx(element_t i,element_t *arrD,pairing_t pairing,element_t *arrQx){
	element_t addTemp;
	element_init_Zr(addTemp,pairing);
	element_t mulTemp;
	element_init_Zr(mulTemp,pairing);
	element_t powTemp;
	element_init_Zr(powTemp,pairing);
	element_set(addTemp,y);
	for(int k = 1; k < d - 1; ++k){
		mpz_t exp;
		mpz_init_set_ui(exp,k);
		element_pow_mpz(powTemp,i,exp);
		element_mul(mulTemp,arrD[k-1],powTemp);
		element_add(addTemp,addTemp,mulTemp);
	}
	element_set(arrQx[qxIndex++],addTemp);
}

//T(x)
void calculateTx(element_t i,pairing_t pairing,element_t *arrTx,element_t *arrN){
	element_t left;
	element_init_G1(left,pairing);
	element_t ZrTemp;
	element_init_Zr(ZrTemp,pairing);
	mpz_t exp;
	mpz_init_set_ui(exp,N-1);
	element_pow_mpz(ZrTemp,i,exp);
	element_pow_zn(left,g2,ZrTemp);

	element_t right;
	element_init_Zr(right,pairing);
	element_set1(right);
	for(int k = 0; k < N; ++k){
		element_t alphai;
		element_init_Zr(alphai,pairing);
		mpz_t gmp_alphai;
		mpz_init_set_ui(gmp_alphai,k+1);
		element_set_mpz(alphai,gmp_alphai);
		element_t innerMul;
		element_init_Zr(innerMul,pairing);
		element_set1(innerMul);
		for(int j = 1; j <= N; ++j){
			if(j == k+1){
				continue;
			}
			element_t alphal;
			element_init_Zr(alphal,pairing);
			mpz_t gmp_alphal;
			mpz_init_set_ui(gmp_alphal,j);
			element_set_mpz(alphal,gmp_alphal);
			element_t numerator;
			element_t denominator;
			element_t quotient;
			element_init_Zr(numerator,pairing);
			element_init_Zr(denominator,pairing);
			element_init_Zr(quotient,pairing);
			element_sub(numerator,i,alphal);
			element_sub(denominator,alphai,alphal);
			element_div(quotient,numerator,denominator);
			element_mul(innerMul,innerMul,quotient);
		}
		element_t powTemp;
		element_init_G1(powTemp,pairing);
		element_pow_zn(powTemp,arrN[k],innerMul);
		element_mul(right,right,powTemp);
	}

	element_t tx;
	element_mul(tx,left,right);
	element_set(arrTx[txIndex++],tx);
}



int main(int argc, char **argv){

	//=============read data from file end===================
	string config = readFile(configPath);
	string PP = readFile(PPPath);
	string MK = readFile(MKPath);

	pairing_t pairing;
	pbc_demo_pairing_init(pairing, argc, argv);

	Json::Reader reader;
	Json::Reader reader1;
	Json::Reader reader2;

	Json::Value value;
	Json::Value value1;
	Json::Value value2;

	// ======================
	if(!reader2.parse(config,value2)){
		cout << "cannot read config file!\n";
		return 0;
	}
	getJsonValueNKey(value2,2,pairing,NULL,NULL);

	element_t *arrN = new element_t[N];
	element_t *arrM = new element_t[M];


	if(!reader.parse(PP,value)){
		cout << "cannot read PP file!\n";
		return 0;
	}
	 getJsonValueNKey(value,0,pairing,arrN,arrM); 


	if(!reader1.parse(MK,value1)){
		cout << "cannot read MK file!\n";
		return 0;
	}
	getJsonValueNKey(value1,1,pairing,arrN,arrM);
	// ======================

	//=============read data from file end===================
	

	//=============generate ri===================
	element_t *randomRi = new element_t[N-1];
	for(int i = 0; i < N-1; ++i){
		element_init_Zr(randomRi[i],pairing);
		element_random(randomRi[i]);
	}
	//=============generate ri===================

	Json::Value root;  
	Json::StyledWriter sw;  
	ofstream os;
	os.open("../../data/extract_data/data");  

	//=============calculate Di===================
	//d
	element_t *arrD = new element_t[d-1];
	for(int i = 0;i < d-1; ++i){
		element_init_Zr(arrD[i],pairing);
		element_random(arrD[i]);
	}
	
	element_t *omiga = new element_t[N-1];
	element_t *arrQx = new element_t[N-1];
	for(int i = 0;i < N-1; ++i){
		element_init_Zr(omiga[i],pairing);
		element_random(omiga[i]);
		element_init_Zr(arrQx[i],pairing);
		calculateQx(omiga[i],arrD,pairing,arrQx);
	}

	element_t *arrTx = new element_t[N-1];

	for(int i = 0;i < N-1; ++i){
		element_init_G1(arrTx[i],pairing);
		calculateTx(omiga[i],pairing,arrTx,arrN);
	}

	//write Tx to file
	for(int i = 0; i < N-1; ++i){
		string index = "";
	    stringstream st;
	    st << (i+1);
	    st >> index;
	    index ="T-" + index;
	    root[index] = Json::Value((char*)transfer(arrTx[i]));  
	}

	//Di
	element_t *arrDi = new element_t[N-1];
	for(int i = 0; i < N-1; ++i){
		element_init_G1(arrDi[i],pairing);
		element_t powTemp;
		element_init_G1(powTemp,pairing);
		element_pow_zn(powTemp,arrTx[i],randomRi[i]);
		element_mul(arrDi[i],arrQx[i],powTemp);
		// element_printf("D-%d:%B\n",i+1,arrDi[i]);
		string index = "";
	    stringstream st;
	    st << (i+1);
	    st >> index;
	    index ="D-" + index;
	    root[index] = Json::Value((char*)transfer(arrDi[i]));  
	}
	//=============calculate Di===================



	//=============calculate di===================
	element_t *arrdi = new element_t[N-1];
	for(int i = 0; i < N-1; ++i){
		element_init_G1(arrdi[i],pairing);
		element_t negRandomRi;
		element_init_Zr(negRandomRi,pairing);
		element_neg(negRandomRi,randomRi[i]);
		element_pow_zn(arrdi[i],g,negRandomRi);
		// element_printf("d-%d:%B\n",i+1,arrdi[i]);
		string index = "";
	    stringstream st;
	    st << (i+1);
	    st >> index;
	    index ="d-" + index;
	    root[index] = Json::Value((char*)transfer(arrdi[i])); 
	}
	//=============calculate di===================
  	

  	os << sw.write(root);  
  	os.close();  

  	 //create file omiga to keep the variables 
  	if(freopen("../../data/extract_data/omiga","w",stdout)==NULL){
    	fprintf(stderr, "error2\n");
  	}
  	printf("{\n");
  	for(int i = 0; i < N - 1; ++i){
  	  	element_printf("\"O-%d\":\"%B\"",i+1,omiga[i]);
		if(i != N-2){
  			element_printf(",\n");
  		}
  		else{
  			element_printf("\n");
  		}
  	}
  	//end writing data to file MK 
  	printf("}\n");

  	//-----------Second----------------
  	//write Qx to file
  	if(freopen("../../data/extract_data/Qx","w",stdout)==NULL){
    	fprintf(stderr, "error\n");
  	}
  	printf("{\n");
  	for(int i = 0; i < N - 1; ++i){
  	  	element_printf("\"Q-%d\":\"%B\"",i+1,arrQx[i]);
		if(i != N-2){
  			element_printf(",\n");
  		}
  		else{
  			element_printf("\n");
  		}
  	}
  	//end writing data to file MK 
  	printf("}\n");
  	fclose(stdout);

	delete []arrdi;
	delete []arrDi;
	delete []randomRi;
	delete []arrTx;
	delete []arrD;
	delete []arrQx;
	delete []omiga;
	delete []arrN;
	delete []arrM;
	pairing_clear(pairing);           
}