#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
const int N = 21;
const int M = 10;
const int d = 9;
static int qxIndex = 0;
static int txIndex = 0;
/**
compile program:  g++ -o test Setup.cpp -L. -lpbc -lgmp
execute Executable file:  ./test < ../../data/param/a.param
**/


//q(i)
void calculateQx(element_t i,element_t *coefficient,pairing_t pairing,element_t *arrQx,element_t y){
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
    element_mul(mulTemp,coefficient[k-1],powTemp);
    element_add(addTemp,addTemp,mulTemp);
  }
  element_set(arrQx[qxIndex++],addTemp);

}


//T(x)
void calculateTx(element_t i,pairing_t pairing,element_t *arrTx,element_t *arrN,element_t g2){
  element_t left;
  element_init_G1(left,pairing);
  element_t ZrTemp;
  element_init_Zr(ZrTemp,pairing);
  mpz_t exp;
  mpz_init_set_ui(exp,N-1);
  element_pow_mpz(ZrTemp,i,exp);
  element_pow_zn(left,g2,ZrTemp);

  element_t right;
  element_init_G1(right,pairing);
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

void generateMLengthBitsString(int arr[]){
  srand((unsigned)time(NULL));
  for(int i = 0; i < M; ++i){
    arr[i] = rand() % 2;
    // arr[i] = 1;
  }
}

void calculateTheFirstPart(element_t v,element_t *arrV,int *message,element_t *randomSi,element_t *arrFirstPart,element_t *arrDi,pairing_t pairing){
  for(int i = 0; i < N-1; i++){
    // element_t init;
    // element_init_G1(init,pairing);
    // element_set1(init);
    element_t mulResult;
    element_init_G1(mulResult,pairing);
    element_set1(mulResult);
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
      element_pow_zn(powTemp,arrV[j],element_temp);
      element_mul(mulResult,mulResult,powTemp);
      // element_set(init,mulResult);
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

void calculatePartThreePartTwo(element_t v,int *message,pairing_t pairing,element_t *arrV,element_t fixedValue){
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
    element_pow_mpz(powTemp,arrV[i],exp);
    element_mul(mulTemp,initValue,powTemp);
    element_set(initValue,mulTemp);
  }
  // fixedValue only has two result depending on the num of 1 or 0
  element_set(fixedValue,mulTemp);
}

int main(int argc, char **argv) {
  //定义变量
  pairing_t pairing;
  pbc_demo_pairing_init(pairing, argc, argv);

  element_t z,y;
  element_t v,g,g1,g2;
  element_t A;
  element_t arrT[N];
  element_t ZVector[M];
  element_t arrV[M];

  element_init_Zr(z,pairing);
  element_init_Zr(y,pairing);
  element_init_G1(v,pairing);
  element_init_G1(g,pairing);
  element_init_G1(g1,pairing);
  element_init_G1(g2,pairing);
  element_init_GT(A,pairing);
  element_init_G1(v,pairing);


  //===========SETUP START===================
  //生成元g
  element_random(g);
  //G中的随机元素g2
  element_random(g2);
  //主私钥y
  element_random(y);
  //g1=g的y次方
  element_pow_zn(g1, g, y);
  //z属于Zp
  element_random(z);
  //计算v'
  element_pow_zn(v, g, z);

  for(int i = 0; i <= N - 1; ++i){
    element_init_G1(arrT[i],pairing);
    element_random(arrT[i]);
  }

  for(int i = 0; i < M; ++i){
    element_init_Zr(ZVector[i],pairing);
    element_random(ZVector[i]);
  }

  for(int i = 0; i < M; ++i){
    element_init_G1(arrV[i],pairing);
    element_pow_zn(arrV[i],g,ZVector[i]);
  }
  pairing_apply(A,g1,g2,pairing);
  //setup得到g,g1,g2,n+1个t属于G1，v‘属于G1，m个v属于G1，A属于GT
  //===========SETUP END===================


  //===========EXTRACT START===================
  //确定di,其中ri为一个属于Zr的随机数
  element_t randomRi[N-1];
  for(int i = 0; i < N - 1; ++i){
    element_init_Zr(randomRi[i],pairing);
    element_random(randomRi[i]);
  }
  element_t arrdi[N-1];
  element_t tempNegRi;
  element_init_Zr(tempNegRi,pairing);
  for(int i = 0; i < N-1; ++i){
    element_init_G1(arrdi[i],pairing);
    element_neg(tempNegRi,randomRi[i]);
    element_pow_zn(arrdi[i],g,tempNegRi);
  }
  element_clear(tempNegRi);

  //确定Di
  element_t *coefficient = new element_t[d-1];
  for(int i = 0;i < d-1; ++i){
    element_init_Zr(coefficient[i],pairing);
    element_random(coefficient[i]);
  }

  element_t *omiga = new element_t[N-1];
  for(int i = 0;i < N-1; ++i){
    element_init_Zr(omiga[i],pairing);
    element_random(omiga[i]);
  }

  element_t *arrQx = new element_t[N-1];
  for(int i = 0; i < N-1; ++i){
    element_init_Zr(arrQx[i],pairing);
    calculateQx(omiga[i],coefficient,pairing,arrQx,y);
  }

  element_t *arrTx = new element_t[N-1];
  for(int i = 0;i < N-1; ++i){
    element_init_G1(arrTx[i],pairing);
    calculateTx(omiga[i],pairing,arrTx,arrT,g2);
  }

  element_t *arrDi = new element_t[N-1];
  for(int i = 0; i < N-1; ++i){
    element_init_G1(arrDi[i],pairing);
    element_t powTemp;
    element_init_G1(powTemp,pairing);
    element_pow_zn(powTemp,arrTx[i],randomRi[i]);
    element_t g2PowTemp;
    element_init_G1(g2PowTemp,pairing);
    element_pow_zn(g2PowTemp,g2,arrQx[i]);
    element_mul(arrDi[i],g2PowTemp,powTemp);
  }
  //===========EXTRACT END===================



  //===========Sign START===================
  int *message = new int[M];
  generateMLengthBitsString(message);
  element_t *randomSi = new element_t[N-1];
  for(int i = 0; i < N - 1; ++i){
    element_init_Zr(randomSi[i],pairing);
    element_random(randomSi[i]);
  }

  element_t *arrFirstPart = new element_t[N-1];
  for(int i = 0; i < N - 1; ++i){
    element_init_G1(arrFirstPart[i],pairing);
  }
  calculateTheFirstPart(v,arrV,message,randomSi,arrFirstPart,arrDi,pairing);

  //第二部分就是di

  element_t *arrThirdPart = new element_t[N-1];
  for(int i = 0; i < N - 1; ++i){
    element_init_G1(arrThirdPart[i],pairing);
    element_t tempRandomSi;
    element_init_Zr(tempRandomSi,pairing);
    element_neg(tempRandomSi,randomSi[i]);
    element_pow_zn(arrThirdPart[i],g,tempRandomSi);
  }
  //===========Sign END===================

  //===========Verify START===================

  element_t fixedValue;
  element_init_G1(fixedValue,pairing);
  calculatePartThreePartTwo(v,message,pairing,arrV,fixedValue);

  element_t outMul;
  element_init_GT(outMul,pairing);
  element_set1(outMul);
  for(int i = 0; i < N-1; ++i){
    element_t result;
    element_init_G1(result,pairing);

    element_t e1;
    element_init_GT(e1,pairing);
    element_t e2;
    element_init_GT(e2,pairing);
    element_t e3;
    element_init_GT(e3,pairing);
    pairing_apply(e1,arrFirstPart[i],g,pairing);
    pairing_apply(e2,arrdi[i],arrTx[i],pairing);
    pairing_apply(e3,arrThirdPart[i],fixedValue,pairing);

    element_t mulTemp1;
    element_init_GT(mulTemp1,pairing);
    element_mul(mulTemp1,e1,e2);
    element_mul(mulTemp1,mulTemp1,e3);


    //calculate delta i,s(0)
    element_t innerMul;
    element_init_Zr(innerMul,pairing);
    element_set1(innerMul);
    for(int j = 0; j < N - 1; ++j){
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
    element_t innerPow;
    element_init_GT(innerPow,pairing);
    element_pow_zn(innerPow,mulTemp1,innerMul);
    element_mul(outMul,outMul,innerPow);
  }

  element_printf("outMul = %B\n",outMul);
  element_printf("A = %B\n",A);
  //===========Verify END===================

  // //===========SETUP CLEAR START===================
  // delete []arrV;
  // delete []arrT;
  // delete []ZVector;
  // element_clear(v);
  // element_clear(y);
  // element_clear(z);
  // element_clear(g);
  // element_clear(g1);
  // element_clear(g2);
  // element_clear(A);
  // //===========SETUP CLEAR END===================

  // //===========EXTRACT CLEAR START===================
  // delete []randomRi;
  // delete []arrdi;
  // delete []arrDi;
  // delete []omiga;
  // delete []coefficient;
  // delete []arrQx;
  // delete []arrTx;
  // //===========EXTRACT CLEAR END===================


  // //===========SIGN CLEAR START===================
  // delete []message;
  // delete []randomSi;
  // delete []arrFirstPart;
  // delete []arrThirdPart;
  // //===========SIGN CLEAR END===================





  pairing_clear(pairing);
  return 0;
}