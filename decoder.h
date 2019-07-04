#ifndef __DECODER_H__
#define __DECODER_H__

#include <map>

typedef unsigned char u8;
typedef unsigned short int u16;

enum HEADER
{
    SOI,
    APPn,
    DQT,
    SOF0,
    DHT,
    SOS,
    COM,
    EOI,
    NotDefined
};

class Scanner
{
public:
    Scanner(std::ifstream* input_obj) : _input_obj(input_obj) {};

    u16 DECODE(std::map<std::string,u16> table);
    u16 RECEIVE(size_t num);

    int EXTEND(int V, size_t T);
    int DPCM(int diff, size_t n);

private:
    std::ifstream* _input_obj;
    std::string _enc;
    std::string _dec;

    int _preDC[3] = {0};
    size_t _ptrL=0;
    size_t _ptrR=8;
};

size_t _zigzag[64] = {  0,  1,  8, 16,  9,  2,  3, 10,
                       17, 24, 32, 25, 18, 11,  4,  5,
                       12, 19, 26, 33, 40, 48, 41, 34,
                       27, 20, 13,  6,  7, 14, 21, 28,
                       35, 42, 49, 56, 57, 50, 43, 36,
                       29, 22, 15, 23, 30, 37, 44, 51,
                       58, 59, 52, 45, 38, 31, 39, 46,
                       53, 60, 61, 54, 47, 55, 62, 63};

u8* _image;

size_t _Nf = 0;
size_t _X  = 0;
size_t _Y  = 0;
size_t _plane1 = 0;
size_t _plane2 = 0;

char*   _Ci;
size_t* _H;
size_t* _V;
size_t* _Q;
size_t* _sx;
size_t* _sy;

size_t _mcuX;
size_t _mcuY;

size_t _idx=0;
size_t _idy=0;
size_t _block_x;
size_t _block_y;

size_t  _Ns;
size_t* _Hk;
size_t* _Vk;
size_t* _Qk;
size_t* _Td;
size_t* _Ta;
size_t* _stepX;
size_t* _stepY;

std::map<std::string,u16> _DC[2];
std::map<std::string,u16> _AC[2];
std::map<size_t,u8*> _YY;
std::map<size_t,u8*> _Cb;
std::map<size_t,u8*> _Cr;

std::map<u8, size_t*> _dqt;
float _IDCT[64][64] = {0};

#endif // define __DECODER_H__
