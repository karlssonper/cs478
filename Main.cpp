/*
 * Main.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: per
 */

#include <GL/glut.h>
#include "FreeImage2Tex.h"
#include <vector>
#include <iostream>
#include <assert.h>
#include "FlashMatting.h"

std::vector<unsigned char > FlashImage;
std::vector<unsigned char > NoFlashImage;
std::vector<unsigned char > TriMapImage;
std::vector<unsigned char > ParisImage;
std::vector<FlashMatting::TrimapValue > TriMap;
std::vector<unsigned char > AlphaImage;
std::vector<unsigned char > CompImage;
int imgWidth;
int imgHeight;

FlashMatting * flashMatting;
int render = 0;

void keyPressed (unsigned char key, int x, int y)
{
    switch (key) {
        case '1' :
            render = 0;
            break;
        case '2' :
            render = 1;
            break;
        case '3' :
            render = 2;
            break;
        case '4' :
            render = 3;
            break;
        case '5' :
            render = 4 ;
            break;
        case '6' :
           render = 5 ;
           break;
        case '7' :
           render = 6 ;
           break;
    }
    glutPostRedisplay();
}
void display()
{
    glMatrixMode(GL_MODELVIEW);
    glClear(GL_COLOR_BUFFER_BIT);
    switch (render) {
        case 0:
            std::cout << "Drawing Flash..." << std::endl;
            glDrawPixels(imgWidth, imgHeight, GL_RGBA, GL_UNSIGNED_BYTE,
                    &FlashImage[0]);
            break;
        case 1:
            std::cout << "Drawing No Flash..." << std::endl;
            glDrawPixels(imgWidth, imgHeight, GL_RGBA, GL_UNSIGNED_BYTE,
                    &NoFlashImage[0]);
            break;
        case 2:
            std::cout << "Drawing Trimap..." << std::endl;
            glDrawPixels(imgWidth, imgHeight, GL_RGBA, GL_UNSIGNED_BYTE,
                    &TriMapImage[0]);
            break;
        case 3:
            AlphaImage = flashMatting->alphaMap();
            flashMatting->updateOneLayer();
            std::cout << "Drawing Alphamap..." << std::endl;
            glDrawPixels(imgWidth, imgHeight, GL_RGBA, GL_UNSIGNED_BYTE,
                    &AlphaImage[0]);
            break;
        case 4:
            AlphaImage = flashMatting->alphaMap();
            while(flashMatting->updateOneLayer());
            std::cout << "Drawing Alphamap..." << std::endl;
            glDrawPixels(imgWidth, imgHeight, GL_RGBA, GL_UNSIGNED_BYTE,
                    &AlphaImage[0]);
            break;
        case 5:
            flashMatting->median();
            AlphaImage = flashMatting->alphaMap();
            std::cout << "Drawing Alphamap..." << std::endl;
            glDrawPixels(imgWidth, imgHeight, GL_RGBA, GL_UNSIGNED_BYTE,
                    &AlphaImage[0]);
            break;
        case 6:
            for (int i = 0; i < imgWidth*imgHeight; ++i) {
                CompImage[4*i] = (1-AlphaImage[4*i]/255.0)*ParisImage[4*i] + AlphaImage[4*i]/255.0 * NoFlashImage[4*i];
                CompImage[4*i+1] = (1-AlphaImage[4*i+1]/255.0)*ParisImage[4*i+1] + AlphaImage[4*i+1]/255.0 * NoFlashImage[4*i+1];
                CompImage[4*i+2] = (1-AlphaImage[4*i+2]/255.0)*ParisImage[4*i+2] + AlphaImage[4*i+2]/255.0 * NoFlashImage[4*i+2];
                CompImage[4*i+3] = 255;
            }
            glDrawPixels(imgWidth, imgHeight, GL_RGBA, GL_UNSIGNED_BYTE,
                    &CompImage[0]);
        default:

            break;
    }
    glutSwapBuffers();
}

int main(int argc, char **argv)
{
    //FreeImage2Texture Flash("../data/cowflash.jpg");
    //FreeImage2Texture NoFlash("../data/cownoflash.jpg");
    //FreeImage2Texture trimap("../data/trimap.png");

    FreeImage2Texture Flash("../data/GT24.png");
    FreeImage2Texture NoFlash("../data/GT24noflash.png");
    FreeImage2Texture trimap("../data/GT24trimap.png");
    FreeImage2Texture paris("../data/paris.jpg");
    imgWidth = Flash.w;
    imgHeight = Flash.h;
    FlashImage = std::vector<unsigned char>(Flash.data,
            Flash.data+(imgWidth * imgHeight*4));
    NoFlashImage = std::vector<unsigned char>(NoFlash.data,
            NoFlash.data+(imgWidth * imgHeight*4));
    TriMapImage = std::vector<unsigned char>(trimap.data,
            trimap.data+(imgWidth * imgHeight*4));
    ParisImage =  std::vector<unsigned char>(paris.data,
            paris.data+(imgWidth * imgHeight*4));

    CompImage.resize(imgWidth * imgHeight*4);
    TriMap.resize(imgWidth * imgHeight);
    for (int i = 0; i < imgWidth*imgHeight; ++i) {
        //std::cout << "LOL!!!" << std::endl;
        if (TriMapImage[i*4] < 25) {
            TriMap[i] = FlashMatting::BACKGROUND;
            TriMapImage[i*4] = 0;
            TriMapImage[i*4+1] = 0;
            TriMapImage[i*4+2] = 0;

        }
        else if (TriMapImage[i*4] > 200) {
            TriMap[i] = FlashMatting::FOREGROUND;
            TriMapImage[i*4] = 255;
            TriMapImage[i*4+1] = 255;
            TriMapImage[i*4+2] = 255;
        }
        else {
            TriMap[i] = FlashMatting::UNKNOWN;
            TriMapImage[i*4] = 128;
            TriMapImage[i*4+1] = 128;
            TriMapImage[i*4+2] = 128;
        }
    }
    assert(FlashImage.size() == 4* imgWidth * imgHeight);

    flashMatting = new FlashMatting(FlashImage,NoFlashImage,imgWidth,imgHeight);
    flashMatting->loadTriMap(TriMap);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(imgWidth, imgHeight);
    glutCreateWindow("cs478 project");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyPressed);
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(0,1,0,1);
    glutMainLoop();
    return 0;
}
