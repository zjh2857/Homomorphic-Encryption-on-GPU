#pragma once

#include "helper.cuh"

__global__ void fillTablePsi128(unsigned long long psi, unsigned long long q, unsigned long long psiinv, unsigned long long psiTable[], unsigned long long psiinvTable[], unsigned int nbit)
{   
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    psiTable[i] = modpow128(psi, bitReverse(i, nbit), q);
    psiinvTable[i] = modpow128(psiinv, bitReverse(i, nbit), q);
}

__global__ void fillTablePsi128Forward(unsigned long long psi, unsigned long long q, unsigned long long psiTable[], unsigned int nbit)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    psiTable[i] = modpow128(psi, bitReverse(i, nbit), q);
}

__global__ void fillTablePsi64(unsigned psi, unsigned q, unsigned psiinv, unsigned psiTable[], unsigned psiinvTable[], unsigned int nbit)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    psiTable[i] = modpow64(psi, bitReverse(i, nbit), q);
    psiinvTable[i] = modpow64(psiinv, bitReverse(i, nbit), q);
}
__global__ void fillTablePsi128Rot(unsigned long long q, unsigned long long psiTable[],unsigned long long psiRotTable[])
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    psiRotTable[i] = modpow128(psiTable[i], 3, q);
}
void getParams(unsigned long long& q, unsigned long long& psi, unsigned long long& psiinv, unsigned long long& ninv, unsigned int& q_bit, unsigned long long n)
{
    if (n == 2048)
    {
        q = 137438691329;
        psi = 22157790;
        psiinv = 88431458764;
        ninv = 137371582593;
        q_bit = 37;
    }
    else if (n == 4096)
    {
        q = 288230376135196673;
        psi = 60193018759093;
        psiinv = 236271020333049746;
        ninv = 288160007391023041;
        q_bit = 58;

        // q = 33538049;
        // psi = 2386;
        // psiinv = 26102329;
        // ninv = 33529861;
        // q_bit = 25;
    }
    else if (n == 8192)
    {
        q = 8796092858369;
        psi = 1734247217;
        psiinv = 5727406356888;
        ninv = 8795019116565;
        q_bit = 43;
    }
    else if (n == 16384)
    {
        q = 281474976546817;
        psi = 23720796222;
        psiinv = 129310633907832;
        ninv = 281457796677643;
        q_bit = 48;
    }
    else if (n == 32768)
    {
        q = 36028797017456641;
        psi = 1155186985540;
        psiinv = 31335194304461613;
        ninv = 36027697505828911;
        q_bit = 55;
    }
}

void getParams30(unsigned& q, unsigned& psi, unsigned& psiinv, unsigned& ninv, unsigned& q_bit, unsigned n)
{
    if (n == 2048)
    {
        /*q = 12931073;
        psi = 3733;
        psiinv = 10610200;
        q_bit = 24;
        ninv = 12924759;*/

        q = 536608769;
        psi = 284166;
        psiinv = 208001377;
        q_bit = 29;
        ninv = 536346753;
    }
    else if (n == 4096)
    {
        q = 33538049;
        psi = 2386;
        psiinv = 26102329;
        ninv = 33529861;
        q_bit = 25;
    }
    else if (n == 8192)
    {
        q = 8716289;
        psi = 1089;
        psiinv = 8196033;
        q_bit = 24;
        ninv = 8715225;
    }
    else if (n == 16384)
    {
        q = 13664257;
        psi = 273;
        psiinv = 8959348;
        q_bit = 24;
        ninv = 13663423;
    }
    else if (n == 32768)
    {
        q = 19070977;
        psi = 377;
        psiinv = 16642842;
        q_bit = 25;
        ninv = 19070395;
    }
    else if (n == 65536)
    {
        q = 13631489;
        psi = 13;
        psiinv = 12582913;
        q_bit = 24;
        ninv = 13631281;
    }
}
// void getParms60(unsigned& q, unsigned& psi, unsigned& psiinv, unsigned& ninv, unsigned& q_bit, unsigned n)
void getParams(unsigned long long *q, unsigned long long* psi,unsigned long long* psiinv, unsigned long long* q_bit,int polylen){
    unsigned long long q_t[8] = {1179649, 1376257, 1769473, 2424833, 2752513, 3604481, 3735553, 5308417};
    unsigned long long psi_t[8] = {1034474, 1172569, 1557013, 349058, 785782, 1977521, 627147, 4561077};
    unsigned long long psiinv_t[8] = {441827, 1160428, 233173, 928390, 2113364, 499635, 1712567, 3973130};
    unsigned long long q_bit_t[8] = {21, 21, 21, 22, 22, 22, 22, 23};
    if(polylen == 4096){

        for(int i = 0; i < 8; i++){
            q[i] = q_t[i];
            psi[i] = psi_t[i];
            psiinv[i] = psiinv_t[i];
            q_bit[i] = q_bit_t[i];
        }
        // q[0] = 288230376135196673;
        // psi[0] = 60193018759093;
        // psiinv[0] = 236271020333049746;
        // // ninv[0] = 288160007391023041;
        // q_bit[0] = 58;
    }else{
        throw "wrong polylen";
    }
}