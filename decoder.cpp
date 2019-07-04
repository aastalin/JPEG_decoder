#include <stdio.h>
#include <cmath>
#include <bitset>
#include <string>
#include <fstream>
#include <iostream>
#include "decoder.h"

#define PI 3.14159265359
#define HUFBIT 16
#define AC_NUM 63
#define BL_NUM 64
#define DEBUG 0

/* Scanner */
u16 Scanner::RECEIVE(size_t num)
{
    if( num==0 ) return 0;

    std::string val="";
    if( (8-_ptrL) >= num )
    {
        val = _enc.substr(_ptrL,num);
        _ptrR = _ptrL = std::min(8,int(_ptrL+num));
    }
    else
    {
        size_t len = 0;
        while( true )
        {
            u8 byte;
            len += 8-_ptrL;
            val += _enc.substr(_ptrL,len);

            _input_obj->read( (char*)&byte, 1);

            // special mark
            if( byte==0xFF )
            {
                _input_obj->read( (char*)&byte, 1);
                if( byte!=0x00 )
                {
                    return byte | 0xFF<<8;
                }
                byte = 0xFF;
            }

            std::bitset<8> bit(byte);
            _enc = bit.to_string();

            // find
            if( (num-len)<=8 )
            {
                _ptrR = _ptrL = num-len;
                val += _enc.substr(0, num-len);
                break;
            }
            else
            {
                _ptrR = _ptrL = 0 ;
            }
        }
    }
    return std::stoi(val.c_str(),0,2);
}

u16 Scanner::DECODE(std::map<std::string,u16> table)
{
    while( true )
    {
        if( _ptrR==8 )
        {
            if( _input_obj->eof() ) break;
 
            u8 byte;
            _input_obj->read( (char*)&byte, 1);

            if(_ptrL<8)
            {
                _dec = _dec+_enc.substr(_ptrL,_ptrR-_ptrL);
            }
            else
            {
                _dec = "";
            }

            // special mark
            if( byte==0xFF )
            {
                _input_obj->read( (char*)&byte, 1);
                if( byte!=0x00 )
                {
                    return byte | 0xFF<<8;
                }
                byte = 0xFF;
            }

            _ptrR = _ptrL =0;

            std::bitset<8> bit(byte);
            _enc = bit.to_string();
        }
        // match huffman
        std::string tmp = _dec+_enc.substr(_ptrL,_ptrR-_ptrL+1);
        if( table.find(tmp)==table.end() )
        {
            _ptrR += 1;
            continue;
        }
        _dec = "";
        _ptrL = _ptrR = _ptrR + 1;
        return table[tmp];
    }
    return 0xFFFF;
}

int Scanner::EXTEND(int V, size_t T)
{
    int Vt = pow(2, T-1);
    return V<Vt ? V+((-1)<<T)+1 : V;
}

int Scanner::DPCM(int diff, size_t n)
{
    _preDC[n] += diff;
    return _preDC[n];
}


/* Utility */
void build_IDCT_table()
{
    // IDCT coefficient
    for( size_t x=0;  x<8; x++ )
    {
        for( size_t y=0; y<8; y++ )
        {
            for( size_t u=0; u<8; u++ )
            {
                for( size_t v=0; v<8; v++ )
                {
                    float CuCv = 1;
                    if( u==0 && v==0 ) CuCv = 0.5;
                    else if( u==0 || v==0 ) CuCv = 0.70710678118;
                    _IDCT[y*8+x][v*8+u] = CuCv * cos(float(2*x+1)*u*PI/16) * cos(float(2*y+1)*v*PI/16);
                }
            }
        }
    }
}

u16 parse_length( std::ifstream* input_obj )
{
    u8 bytes[2];
    u16 lp;

    input_obj->read( (char*)bytes, 2 );
    lp = bytes[1] | (bytes[0] << 8);

    if( lp==0 )
    {
        std::cerr << "data length = 0" << std::endl;
        exit(-1);
    }

    return lp;
}

u8* IDCT(int* input, size_t x, size_t y, size_t n, size_t v, size_t h)
{
    u8* output = new u8[BL_NUM];
    for( size_t ix=0; ix<8; ix++ ) {
    for( size_t iy=0; iy<8; iy++ ) {
        // calcuate IDCT
        float val = 0;
        for( size_t uv=0; uv<BL_NUM; uv++ )
        {
            val += (input[uv]) * _IDCT[iy*8+ix][uv];
        }
        output[iy*8+ix] = u8(std::min(std::max(int(std::floor(0.25*val)+128),0),255));
    }
    }

    // YCbCr data
    if( n==0 )
    {
        size_t prefix = (_Hk[0]*y+h)+ _Hk[0]*_block_x*(_Vk[0]*x+v);
        _YY[prefix] = output;
    }
    else if( n==1 )
    {
        size_t prefix = y + _block_x*x;
        _Cb[prefix] = output;
    }
    else
    {
        size_t prefix = y + _block_x*x;
        _Cr[prefix] = output;
    }
}

void output_bmp()
{
    u8 header[54] = {
        0x42, 0x4d, 0, 0, 0, 0, 0, 0, 0, 0,
        54, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0
    };
    u8 bmppad[3] = {0,0,0};

    // BMP header
    int size = 54+3*_X*_Y;
    header[ 2] = u8(size    );
    header[ 3] = u8(size>> 8);
    header[ 4] = u8(size>>16);
    header[ 5] = u8(size>>24);
    header[18] = u8(  _X    );
    header[19] = u8(  _X>> 8);
    header[20] = u8(  _X>>16);
    header[21] = u8(  _X>>24);
    header[22] = u8(  _Y    );
    header[23] = u8(  _Y>> 8);
    header[24] = u8(  _Y>>16);
    header[25] = u8(  _Y>>24);

    // initial
    for( int i=0; i<_Nf*_X*_Y; i++ )
    {
      _image[i] = 0;
    }

    for( int bx=0; bx<_block_x; bx++ ) {
    for( int by=0; by<_block_y; by++ ) {
        for( int bh=0; bh<_Vk[0]; bh++ ) {
        for( int bv=0; bv<_Hk[0]; bv++ ) {

            u8* YY = _YY[(_Hk[0]*bx+bv)+(_Hk[0]*_block_x)*(_Vk[0]*by+bh)];
            u8* Cb = _Cb[bx+_block_x*by];
            u8* Cr = _Cr[bx+_block_x*by];

            for( int lx=0; lx<8; lx++ ) {
            for( int ly=0; ly<8; ly++ ) {

                size_t ix = 8*(_Hk[0]*bx+bv) + lx;
                size_t iy = 8*(_Vk[0]*by+bh) + ly;
                if(ix>=_X || iy>=_Y) continue;

                size_t rid = lx+8*ly;
                size_t cid = ((8/_Hk[0])*bv+lx/_Hk[0])+8*((8/_Vk[0])*bh+ly/_Vk[0]);

                // YCbCr to RGB
                int b = YY[rid] +    1.772*(Cb[cid]-128);
                int g = YY[rid] - 0.344136*(Cb[cid]-128) - 0.714136*(Cr[cid]-128);
                int r = YY[rid] +    1.402*(Cr[cid]-128);
                // bgr
                _image[(ix+_X*iy)*3+0] = std::max(std::min( b, 255), 0);
                _image[(ix+_X*iy)*3+1] = std::max(std::min( g, 255), 0);
                _image[(ix+_X*iy)*3+2] = std::max(std::min( r, 255), 0);
            }
            }
        }
        }
    }
    }

    // output bmp
    FILE *f;
    f = fopen("output.bmp","wb");
    size_t pad = (4-(_X*3)%4)%4;

    fwrite(header,sizeof(u8),54,f);
    for( size_t i=0; i<_Y; i++ )
    {
        fwrite(_image+_X*(_Y-i-1)*3, 3, _X,f);
        fwrite(bmppad, 1, pad, f);
    }
    fclose(f);

}


/* MCU */
u16 parse_mcu( Scanner* scanner )
{
    for( size_t x=_idx; x<_block_y; x++ ) {
    for( size_t y=_idy; y<_block_x; y++ ) {

        if( DEBUG ) std::cout << "MCU=(" << x<< "," << y << ")" << std::endl;
        for( size_t n=0; n<_Ns; n++ )
        {
            for( size_t v=0; v<_Vk[n]; v++ ) {
            for( size_t h=0; h<_Hk[n]; h++ ) {

                if( DEBUG ) std::cout << "componenet (" << n << "," << v << "," << h << ")"<< std::endl;
                int block[BL_NUM] = {0};

                // DC first
                u16 T = scanner->DECODE(_DC[_Td[n]]);
                u16 V = scanner->RECEIVE(T);
                block[0] = scanner->DPCM( scanner->EXTEND(V, T)*_dqt[_Qk[n]][0], n);
                if( DEBUG ) std::cout << "T=" << T << " V=" << V << " diff=" << block[0] << std::endl;

                // AC
                for( size_t K=0; K<AC_NUM; )
                {
                    u8 RS = scanner->DECODE(_AC[_Ta[n]]);
                    u8 SSSS = RS & 0x0F;
                    u8 R = RS >> 4;
                    if( SSSS==0 )
                    {
                        if( R==15 ) K += 16;
                        else break;
                    }
                    else
                    {
                        K += R;
                        int ZZ = scanner->RECEIVE(SSSS);
                        block[_zigzag[K+1]] = scanner->EXTEND(ZZ,SSSS) * _dqt[_Qk[n]][_zigzag[K+1]];
                        K += 1;
                        if( DEBUG ) printf("RS=%x SSSS=%x R=%x XX=%d diff=%d\n", RS,SSSS,R,ZZ,scanner->EXTEND(ZZ,SSSS));
                    }
                }
                // inverse DCT
                IDCT(block, x, y, n, v, h);
            }
            }
        }
    }
    }
    return 0xFFFF;
}


/* Comment */
void parse_com( std::ifstream* input_obj )
{
    u16 lp = parse_length(input_obj);

    // read data
    u8 data[lp-2];
    input_obj->read( (char*)data, lp-2 );
}


/* APPn segments */
void parse_appn( std::ifstream* input_obj )
{
    u16 lp = parse_length(input_obj);

    // read data
    u8 data[lp-2];
    input_obj->read( (char*)data, lp-2 );
}


/* Quantization table */
void parse_dqt( std::ifstream* input_obj )
{
    u16 lq = parse_length(input_obj);

    // read data
    u8 data[lq-2];
    input_obj->read( (char*)data, lq-2);

    size_t times = (int(lq)-2) / 65;

    // create quantization table
    for( size_t tt=0; tt<times; tt++ )
    {
        u8 idx = data[tt*65+0] & 0x0F;
        size_t* table = new size_t[BL_NUM];

        for( int i=0; i<BL_NUM; i++ )
        {
            table[_zigzag[i]] = data[tt*65+1+i];
        }
        _dqt[idx] = table;
    }
}


/* Baseline DCT */
void parse_sof0( std::ifstream* input_obj )
{
    u16 lf = parse_length(input_obj);

    // read header
    u8 header[lf-2];
    input_obj->read( (char*)header, lf-2 );

    _Y  = header[2] | (header[1] << 8);
    _X  = header[4] | (header[3] << 8);
    _Nf = header[5];

    _Ci = new char[_Nf];
    _sx = new size_t[_Nf];
    _sy = new size_t[_Nf];
    _H  = new size_t[_Nf];
    _V  = new size_t[_Nf];
    _Q  = new size_t[_Nf];
    size_t maxH = 0;
    size_t maxV = 0;

    for( size_t i=0; i<_Nf; i++)
    {
        _Ci[i] = header[6+i*3];
        _H[i]  = header[7+i*3] >> 4;
        _V[i]  = header[7+i*3] & 0x0F;
        _Q[i]  = header[8+i*3];
	maxH = std::max(maxH, _H[i]);
        maxV = std::max(maxV, _V[i]);
    }
    for( size_t i=0; i<_Nf; i++)
    {
        _sx[i] = maxH / _H[i];
        _sy[i] = maxV / _V[i];
    }

    _idx = 0;
    _idy = 0;
    _block_x = std::ceil( float(_X)/(maxH*8) );
    _block_y = std::ceil( float(_Y)/(maxV*8) );

    _plane1 = _X * _Y;
    _plane2 = 2 * _X * _Y;
    _image = new u8[_Nf*_Y*_X];

    _mcuX = 8 * maxH;
    _mcuY = 8 * maxV;
}


/* Huffman table */
void parse_dht( std::ifstream* input_obj )
{
    u16 lh = parse_length(input_obj);

    // read data
    u8 data[lh-2];
    input_obj->read( (char*)data, lh-2);

    size_t count = 0;
    while( true )
    {
        size_t Tc = data[count+0] >> 4;
        size_t Th = data[count+0] & 0x0F;

        size_t sum[HUFBIT];
        for( size_t i=0; i<HUFBIT; i++ )
        {
            if( i==0 )
            {
                sum[i] = data[count+i+1];
            }
            else
            {
                sum[i] = sum[i-1]+data[count+i+1];
            }
        }

        u16 code = 0x00;
        size_t pre_len = 1;
        size_t cur_len = 1;
        for( size_t i=0; i<sum[HUFBIT-1]; i++ )
        {
            for( size_t j=cur_len-1; j<HUFBIT; j++ )
            {
                if( sum[j]>i )
                {
                    cur_len = j+1;
                    break;
                }
            }

            // first codeword
            if( i==0 ) 
            {
                code = 0x00;
            }
            else
            {
                code += 1;
                if( cur_len!=pre_len )
                {
                    code = code << (cur_len-pre_len);
                }
            }

            // create huffman table
            std::bitset<16> codebit( code<<(HUFBIT-cur_len) );
            std::string codeword = codebit.to_string().substr(0,cur_len);
            pre_len = cur_len;
            if( Tc==0 )
            {
                _DC[Th][codeword] = data[count+i+17];
            }
            else
            {
                _AC[Th][codeword] = data[count+i+17];
            }
        }

        // handle multiple table
        count += sum[HUFBIT-1]+17;
        if( count >= (lh-2) ) break;
    }
}


/* Scan */
u16 parse_sos( std::ifstream* input_obj )
{
    u16 ls = parse_length(input_obj);

    u8 header[ls-2];
    input_obj->read( (char*)header, ls-2);

    _Ns = header[0];
    _Hk = new size_t[_Ns];
    _Vk = new size_t[_Ns];
    _Qk = new size_t[_Ns];
    _Td = new size_t[_Ns];
    _Ta = new size_t[_Ns];
    _stepX = new size_t[_Ns];
    _stepY = new size_t[_Ns];
    for( size_t i=0; i<_Ns; i++ )
    {
        char Ck = header[1+2*i];
	u8 Td = header[2+2*i] >> 4;
        u8 Ta = header[2+2*i] & 0x0F;
        for( size_t j=0; j<_Nf;j++ )
        {
            if( Ck==_Ci[j] )
            {
		_Hk[i] = _H[j];
                _Vk[i] = _V[j];
                _Qk[i] = _Q[j];
                _Td[i] = Td;
                _Ta[i] = Ta;
                _stepX[i] = _sx[j];
                _stepY[i] = _sy[j];
                break;
            }
        }
    }

    u16 val;
    bool valid = true;
    Scanner scanner(input_obj);
    while( valid )
    {
        val = parse_mcu(&scanner);
        switch( val )
        {
            case 0xFFD0:
                std::cout << "   - Restart with modulo 8 count 0" << std::endl;
                break;
            case 0xFFD1:
                std::cout << "   - Restart with modulo 8 count 1" << std::endl;
                break;
            case 0xFFD2:
                std::cout << "   - Restart with modulo 8 count 2" << std::endl;
                break;
            case 0xFFD3:
                std::cout << "   - Restart with modulo 8 count 3" << std::endl;
                break;
            case 0xFFD4:
                std::cout << "   - Restart with modulo 8 count 4" << std::endl;
                break;
            case 0xFFD5:
                std::cout << "   - Restart with modulo 8 count 5" << std::endl;
                break;
            case 0xFFD6:
                std::cout << "   - Restart with modulo 8 count 6" << std::endl;
                break;
            case 0xFFD7:
                std::cout << "   - Restart with modulo 8 count 7" << std::endl;
                break;
            default:
                valid = false;
                break;
        }
    }
    return val;
}


/* MAIN */
int main(int argc, char *argv[]) {

    // check argument
    if( argc<2 )
    {
        std::cout << "Usage: " << argv[0] << " [input]" << std::endl;
        exit(-1);
    }
    std::cout << "JPEG decoder written by Aasta Lin" << std::endl;
    std::cout << " - file name: " << argv[1] << std::endl;

    // read file
    std::ifstream input_obj(argv[1], std::ios::binary);
    if( !input_obj.is_open() )
    {
        std::cout << "cannot open file " << argv[1] << std::endl;
        exit(-1);
    }
    build_IDCT_table();

    u16 val;
    u8 bytes[2];
    while( true )
    {
        HEADER header;

        if( val!=0xffd9 )
        {
            input_obj.read( (char*)bytes, 2 );
            val = bytes[1] | (bytes[0] << 8);
        }

        printf("%x\n", val);
        switch( val )
        {
            case 0xFFD8:
                header = HEADER::SOI;
                std::cout << "   - Start of image" << std::endl;
                break;
            case 0xFFE0:
                header = HEADER::APPn;
                std::cout << "   - Reserved for application segments" << std::endl;
                parse_appn( &input_obj );
                break;
            case 0xFFDB:
                header = HEADER::DQT;
                std::cout << "   - Define quantization table(s)" << std::endl;
                parse_dqt( &input_obj );
                break;
            case 0xFFC0:
                header = HEADER::SOF0;
                std::cout << "   - Baseline DCT" << std::endl;
                parse_sof0( &input_obj );
                break;
            case 0xFFC4:
                header = HEADER::DHT;
                std::cout << "   - Define Huffman table(s)" << std::endl;
                parse_dht( &input_obj );
                break;
            case 0xFFDA:
                header = HEADER::SOS;
                std::cout << "   - Start of scan" << std::endl;
                val = parse_sos( &input_obj );
                break;
            case 0xFFFE:
                header = HEADER::COM;
                std::cout << "   - Comment" << std::endl;
                parse_com( &input_obj );
                break;
            case 0xFFD9:
                header = HEADER::EOI;
                std::cout << "   - End of image" << std::endl;
                output_bmp();
                break;
            default:
                header = HEADER::NotDefined;
        }

        if( header==HEADER::EOI ) break;
        if( header==HEADER::NotDefined )
        {
            std::cerr << "Not defined marker" << std::endl;
            break;
        }
    }

    return 0;
}
