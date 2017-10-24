#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include <malloc.h>
#include <stdio.h>

/**
compile program:  gcc -o test Setup.c -L. -lpbc -lgmp
execute Executable file:  ./test < ../../data/param/a.param
**/

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

  //assume n=21,m=10
  element_t arrN[21];
  element_t arrM[10];

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

  //create file PP to keep the variables 
  if(freopen("../../data/setup_data/PP","w",stdout)==NULL){
    fprintf(stderr, "error\n");
  }

  printf("{\n");
  element_printf("\"g1\":\"%B\",\n",g1);
  element_printf("\"g2\":\"%B\",\n",g2);

  //generate elements from G
  for(int i = 0; i <= 20; ++i){
  	element_init_G1(arrN[i],pairing);
  	element_random(arrN[i]);
    element_printf("\"t-%d\":\"%B\",\n",i+1,arrN[i]);
  }



  element_printf("\"v'\":\"%B\",\n",v);
  //generate elements from Zr
  for(int i = 0; i < 10; ++i){
    element_init_G1(arrM[i],pairing);
  	element_random(temp);
  	element_pow_zn(arrM[i],g,temp);
    element_printf("\"v-%d\":\"%B\"",i+1,arrM[i]);
    if(i != 9){
      printf(",\n");
    }
    else{
      printf("\n");
    }
  }
  //end writing data to file PP 
  printf("}\n");
  fclose(stdout);

  //create file MK to keep the variables 
  if(freopen("../../data/setup_data/MK","w",stdout)==NULL){
    fprintf(stderr, "error\n");
  }
  printf("{\n");
  element_printf("\"y\":\"%B\"\n",y);
  //end writing data to file MK 
  printf("}\n");

  fclose(stdout);

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