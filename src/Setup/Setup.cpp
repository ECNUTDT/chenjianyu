#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include "/usr/include/json/json.h"

#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


const int N = 21;
const int M = 10;


/**
compile program:  g++ -o test Setup.cpp -L. -lpbc -lgmp -ljson

execute Executable file:  ./test < ../../data/param/a.param
**/



unsigned char* transfer(element_t t){

  int leng = element_length_in_bytes(t);
  unsigned char *p = new unsigned char[leng + 1];
  int writeLength = element_to_bytes_compressed(p,t);
  return p;
}

int main(int argc, char **argv) {

  //define some variables
  pairing_t pairing;
  //random Integer
  element_t z,y;
  //elements from G
  element_t v,g,g1,g2;
  // a temp element from G
  element_t temp;

  //read a.param file when executing Executable file
  pbc_demo_pairing_init(pairing, argc, argv);

  element_t arrN[N];
  element_t arrM[M];

  //initilize variables
  element_init_G1(v,pairing);
  element_init_Zr(y,pairing);
  element_init_Zr(z,pairing);
  element_init_G1(g,pairing);
  element_init_G1(g1,pairing);
  element_init_G1(g2,pairing);
  element_init_Zr(temp,pairing);

  //generate some random necessary variables
  element_random(y);
  element_random(z);
  element_random(g);
  element_random(g2);
  //compute g1 and v
  element_pow_zn(g1, g, y);
  element_pow_zn(v, g, z);


  //根节点  
  Json::Value root;  
  Json::Value root1;  

  //根节点属性  
  Json::StyledWriter sw;  
  Json::StyledWriter sw1;  

  //输出到文件  
  ofstream os;
  ofstream os1;  

  os.open("../../data/setup_data/PP");  

  root["g"] = Json::Value((char*)transfer(g));  
  root["g1"] = Json::Value((char*)transfer(g1));  
  root["g2"] = Json::Value((char*)transfer(g2));  
  //generate elements from G
  for(int i = 0; i <= N - 1; ++i){
    element_init_G1(arrN[i],pairing);
    element_random(arrN[i]);
    string index = "";
    stringstream st;
    st << (i+1);
    st >> index;
    index ="t-" + index;
    root[index] = Json::Value((char*)transfer(arrN[i]));  
  }


  root["v"] = Json::Value((char*)transfer(v));  
  //generate elements from Zr
  for(int i = 0; i < M; ++i){
    element_init_G1(arrM[i],pairing);
    element_random(temp);
    element_pow_zn(arrM[i],g,temp);
    string index = "";
    stringstream st;
    st << (i+1);
    st >> index;
    index ="v-" + index;
    root[index] = Json::Value((char*)transfer(arrM[i]));  
  }

  os << sw.write(root);  
  os.close();  


 //create file MK to keep the variables 
  if(freopen("../../data/setup_data/MK","w",stdout)==NULL){
    fprintf(stderr, "error2\n");
  }
  printf("{\n");
  element_printf("\"y\":\"%B\"\n",y);
  //end writing data to file MK 
  printf("}\n");

  fclose(stdout);

  os1.open("../../data/config/config");  
  root1["N"] = Json::Value(N);  
  root1["M"] = Json::Value(M);  
  os1 << sw.write(root1);  
  os1.close();  


  //memory management
  element_clear(v);
  element_clear(y);
  element_clear(z);
  element_clear(temp);
  element_clear(g);
  element_clear(g1);
  element_clear(g2);
  pairing_clear(pairing);
  return 0;
}