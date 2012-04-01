/*
 * FreeImage2Tex.h
 *
 *  Created on: Mar 6, 2012
 *      Author: per
 */

#ifndef FREEIMAGE2TEX_H_
#define FREEIMAGE2TEX_H_

#include <FreeImage.h>


#include <iostream>
#include <fstream>

class FreeImage2Texture
{
public:
    FreeImage2Texture(const std::string & _file)
    {
        //Automatocally detects the format(from over 20 formats!)
        FREE_IMAGE_FORMAT formato = FreeImage_GetFileType(_file.c_str(),0);
        FIBITMAP* imagen = FreeImage_Load(formato, _file.c_str());

        FIBITMAP* temp = imagen;
        imagen = FreeImage_ConvertTo32Bits(imagen);
        FreeImage_Unload(temp);

        w = FreeImage_GetWidth(imagen);
        h = FreeImage_GetHeight(imagen);
        //cout<<"The size of the image is: "<<textureFile<<" es "<<w<<"*"<<h<<endl; //Some debugging code

        data = new unsigned char[4*w*h];
        char* pixeles = (char*)FreeImage_GetBits(imagen);
        //FreeImage loads in BGR format, so you need to swap some bytes(Or use GL_BGR).

        //std::ofstream myfile;
        //myfile.open ("texture.txt");

        for(int j= 0; j<w*h; j++){
            data[j*4+0]= pixeles[j*4+2];
            data[j*4+1]= pixeles[j*4+1];
            data[j*4+2]= pixeles[j*4+0];
            data[j*4+3]= pixeles[j*4+3];

            //myfile << int(data[j*4+0]) << " " << int(data[j*4+1]) <<
            //        " " << int(data[j*4+2]) << " ";
        }
          //myfile.close();

        FreeImage_Unload(imagen);
    }

    ~FreeImage2Texture()
    {
        delete[] data;
    }
    unsigned char* data;
    int w,h;
};

#endif /* FREEIMAGE2TEX_H_ */
