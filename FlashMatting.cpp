/*
 * FlashMatting.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: per
 */

#include "FlashMatting.h"
#include <stdio.h>
#include <algorithm>

//#define PRINT_DEBUG

FlashMatting::FlashMatting( const std::vector<unsigned char> & flash,
                            const std::vector<unsigned char> & noFlash,
                            int width,
                            int height) :
                            unknowns(Compare(width))
{
    N = 5;
    Epsilon = 0.01;
    numUnknowns = 0;
    SamplesPerUnknown = 50;
    IterationsPerUnknown = 5;

    SigmaI2 = 32;
    SigmaIp2 = 32;

    nPosAlpha = 0;
    nNegAlpha = 0;
    maxAlpha = -100000.0f;
    minAlpha = 100000.0f;

    Pixel minF = Pixel(10000.f, 10000.f,10000.f);
    Pixel minB = Pixel(10000.f, 10000.f,10000.f);
    Pixel minFP = Pixel(10000.f, 10000.f,10000.f);
    Pixel maxF = Pixel(-10000.f, -10000.f,-10000.f);
    Pixel maxB = Pixel(-10000.f, -10000.f,-10000.f);
    Pixel maxFP = Pixel(-10000.f, -10000.f,-10000.f);

    W = width;
    H = height;
    //If.set(flash,W,H);
    //I.set(noFlash,W,H);

    I.reserve(W,H);
    If.reserve(W,H);
    Ip.reserve(W,H);
    for (int i = 0; i < W; ++i) {
        for (int j = 0; j < H; ++j) {
            I(i,j) = Pixel(noFlash[4*(i + width*j)    ],
                            noFlash[4*(i + width*j) + 1],
                            noFlash[4*(i + width*j) + 2]);

            Ip(i,j) = Pixel(flash[4*(i + width*j)    ],
                            flash[4*(i + width*j) + 1],
                            flash[4*(i + width*j) + 2]);

            //Ip(i,j) = If(i,j) - I(i,j);
        }
    }
}

void FlashMatting::median()
{
    std::vector<float> temp;
    temp.resize(W*H);
    int S = 1;
    for (int i = 0; i < W; ++i) {
        for (int j = 0; j < H; ++j) {
            std::vector<float> v;
            for (int x = std::max(0, i - S); x < std::min(W, i + S); ++x) {
                for (int y = std::max(0, j - S); y < std::min(H,j + S); ++y) {
                    v.push_back(Alpha(x,y));
                }
            }
            std::sort(v.begin(), v.end());
            temp[i + W * j] = v[v.size()/2];
        }
    }

    for (int i = 0; i < W; ++i) {
        for (int j = 0; j < H; ++j) {
            Alpha(i,j) = temp[i + W * j];
        }
    }
}

void FlashMatting::loadTriMap(const std::vector<FlashMatting::TrimapValue> &t)
{
    Trimap.set(t, W, H);
    numUnknowns = 0;

    //Allocate space
    F.reserve(W,H);
    B.reserve(W,H);
    Fp.reserve(W,H);
    Alpha.reserve(W,H);

    //Fill initiation data
    for (int i = 0; i < W; ++i) {
        for (int j = 0; j < H; ++j) {
            if (IsForeground(i,j)) {
                F(i,j) = I(i,j);
                Fp(i,j) = Ip(i,j);
                B(i,j) = Pixel(0.0f, 0.0f, 0.0f);
                Alpha(i,j) = 1.0f;
            } else if(IsBackground(i,j)) {
                F(i,j) = Pixel(0.0f, 0.0f, 0.0f);
                Fp(i,j) = Pixel(0.0f, 0.0f, 0.0f);
                B(i,j) = I(i,j);
                Alpha(i,j) = 0.0f;
            } else {
                F(i,j) = Pixel(0.0f, 0.0f, 0.0f);
                Fp(i,j) = Pixel(0.0f, 0.0f, 0.0f);
                B(i,j) = Pixel(0.0f, 0.0f, 0.0f);
                Alpha(i,j) = 0.5f;
                //Add to set of unknowns
                unknowns.insert(PixelPair(i,j));
                ++numUnknowns;
            }
        }
    }
}

bool FlashMatting::updateOneLayer()
{
    if (unknowns.empty()) {
        printf("Max alpha: %f\n", maxAlpha);
        printf("Min alpha: %f\n", minAlpha);

        if (nPosAlpha)
            printf("Avg positive alpha: %f\n", avgMaxAlpha/nPosAlpha);
        if (nNegAlpha)
            printf("Avg negative alpha: %f\n", avgMinAlpha/nNegAlpha);

        std::cout << std::endl << minF << std::endl;
        std::cout << std::endl << minB << std::endl;
        std::cout << std::endl << minFP << std::endl;
        std::cout << std::endl << maxF << std::endl;
        std::cout << std::endl << maxB << std::endl;
        std::cout << std::endl << maxFP << std::endl;

        return false;
    }

    //Add all unknown pixels with at least one neighbor
    PixelPairSet::const_iterator it;
    for (it = unknowns.begin(); it != unknowns.end(); ++it) {
        if (NeighborsCross(it->first, it->second) > 0) {
            queue.push(*it);
        }
    }

    while (!queue.empty()) {
        //Update Alpha and F B Fp at position pp
        PixelPair pp = queue.front();
        UpdateColors(pp.first, pp.second);
        queue.pop();
        unknowns.erase(pp);
        Trimap(pp.first,pp.second) = UPDATED_UNKNOWN;
        if (unknowns.size() % 10000 == 0)
            std::cerr << unknowns.size() << " left" << std::endl;
    }


    return true;

}

std::vector<unsigned char> FlashMatting::alphaMap()
{
    //Convert from our 32-bit single precision alpha map to 8bit one.
    std::vector<unsigned char> alphaMap;
    alphaMap.resize(W*H*4);
    for (int i = 0; i < W; ++i) {
        for (int j = 0; j < H; ++j) {
            alphaMap[4*(i + W * j) + 0] = static_cast<unsigned char>(Alpha(i,j)*255);
            alphaMap[4*(i + W * j) + 1] = static_cast<unsigned char>(Alpha(i,j)*255);
            alphaMap[4*(i + W * j) + 2] = static_cast<unsigned char>(Alpha(i,j)*255);
            alphaMap[4*(i + W * j) + 3] = 255;
        }
    }
    return alphaMap;
}


int FlashMatting::Neighbors(const int i, const int j) const
{
    int n = 0;
    for (int x = std::max(0, i - N); x < std::min(W, i + N); ++x) {
        for (int y = std::max(0, j - N); y < std::min(H,j + N); ++y) {
            if (x == i && y == j) continue;
            n += Trimap(x,y) > 0 ? 1 : 0;
        }
    }
    return n;
}

int FlashMatting::NeighborsCross(const int i, const int j) const
{
   int n = 0;

   n += Trimap(std::min(W - 1, i + 1), std::max(0    , j - 1))==FOREGROUND ||
            Trimap(std::min(W - 1, i + 1), std::max(0    , j - 1))==UPDATED_UNKNOWN
            ? 1 : 0;
    n += Trimap(std::min(W - 1, i + 1), std::min(H - 1, j + 1)) == FOREGROUND ||
            Trimap(std::min(W - 1, i + 1), std::min(H - 1, j + 1)) == UPDATED_UNKNOWN
            ? 1 : 0;
    n += Trimap(std::max(0    , i - 1), std::max(0    , j - 1)) == FOREGROUND ||
            Trimap(std::max(0    , i - 1), std::max(0    , j - 1)) == UPDATED_UNKNOWN ? 1 : 0;
    n += Trimap(std::max(0    , i - 1), std::min(H - 1, j + 1))==FOREGROUND ||
            Trimap(std::max(0    , i - 1), std::min(H - 1, j + 1))==UPDATED_UNKNOWN  ? 1 : 0;
    return n;
}

void FlashMatting::UpdateColors(const int i, const int j)
{
    //std::cerr << "This guy has: " << NeighborsCross(i,j) << std::endl;

    //InitializeAlpha(i,j);

    std::vector<PixelPair> neighborsF;
    std::vector<PixelPair> neighborsB;
    FindNeighborPoints<FOREGROUND>(neighborsF, i, j);
    FindNeighborPoints<BACKGROUND>(neighborsB, i, j);

    UpdateWeightedMeanAndCovariance<CHAN_F>(neighborsF, i, j);
    UpdateWeightedMeanAndCovariance<CHAN_B>(neighborsB, i, j);
    UpdateWeightedMeanAndCovariance<CHAN_FP>(neighborsF, i, j);

    //Start by initiating F B and FP to Weighted Mean Values
    F(i,j) = curMean[CHAN_F];
    B(i,j) = curMean[CHAN_B];
    Fp(i,j) = curMean[CHAN_FP];

    UpdateAlpha(i,j);

    float deltaAlpha;
    float prevAlpha;
    //do {
    for (int iter = 0; iter < IterationsPerUnknown; ++iter) {
        prevAlpha = Alpha(i,j);
        UpdateChannels(i,j);
        ClampPixel(F(i,j));
        ClampPixel(B(i,j));
        ClampPixel(Fp(i,j));
        UpdateAlpha(i,j);
        ClampAlpha(Alpha(i,j));
        deltaAlpha = Alpha(i,j) - prevAlpha;
        //std::cout << "Delta: " << deltaAlpha << std::endl;
    }
    //while (abs(deltaAlpha) > Epsilon);


}

void FlashMatting::InitializeAlpha(const int i, const int j)
{
    int n = 0;
    float meanAlpha = 0.f;
    for (int x = std::max(0, i - N); x < std::min(W, i + N); ++x) {
        for (int y = std::max(0, j - N); y < std::min(H,j + N); ++y) {
            if (x == i && y == j) continue;
            if (!IsUnknown(x,y)){
                meanAlpha += Alpha(x,y);
                ++n;
            }
        }
    }
    assert(n != 0);
    Alpha(i,j) = meanAlpha / n;
}

void FlashMatting::UpdateAlpha(const int i, const int j)
{
    const Pixel F_minus_B = F(i,j) - B(i,j);

    //Num
    const float a = F_minus_B.transpose() * (I(i,j)-B(i,j));
    const float b = Fp(i,j).transpose() * Ip(i,j);
    const float num = (SigmaIp2 * a + SigmaI2 * b);

    //Den
    const float c = F_minus_B.transpose() * F_minus_B;
    const float d = Fp(i,j).transpose()*Fp(i,j);
    const float den = (SigmaIp2 * c + SigmaI2 * d);

    Alpha(i,j) = num / den;

#ifdef PRINT_DEBUG
    std::cout << "F: " << std::endl << F(i,j) << std::endl << std::endl;
    std::cout << "B: " << std::endl << B(i,j) << std::endl << std::endl;
    std::cout << "FP: " << std::endl << Fp(i,j) << std::endl << std::endl;
    std::cout << "I: " << std::endl << I(i,j) << std::endl << std::endl;
    std::cout << "Ip: " << std::endl << Ip(i,j) << std::endl << std::endl;
    std::cout << "Alpha: " << Alpha(i,j) << std::endl << std::endl;
#endif
}

void FlashMatting::UpdateChannels(const int i, const int j)
{
    const Matrix3 identity = Matrix3::Identity();
    const Matrix3 K = identity * Alpha(i,j) * Alpha(i,j) / SigmaI2;
    const Matrix3 L = identity * Alpha(i,j) * (1 - Alpha(i,j)) * SigmaI2;
    Eigen::MatrixXf A(9,9);
    A.block<3,3>(0,0) = curCovInv[CHAN_F] + K;
    A.block<3,3>(0,3) = L;
    A.block<3,3>(0,6) = Matrix3::Zero();
    A.block<3,3>(3,0) = L;
    A.block<3,3>(3,3) = curCovInv[CHAN_B] + K;
    A.block<3,3>(3,6) = Matrix3::Zero();
    A.block<3,3>(6,0) = Matrix3::Zero();
    A.block<3,3>(6,3) = Matrix3::Zero();
    A.block<3,3>(6,6) = curCovInv[CHAN_FP] + K;

//#ifdef PRINT_DEBUG
    for (int x = 0; x < 9; ++x)
        for (int y = 0; y < 9; ++y){
            if (A(x,y) != A(x,y)){
                std::cout << "YEP har ocksa, for " << i << " " << j << std::endl;
                std::cout << curCovInv[CHAN_F] << std::endl;
                std::cout << curCovInv[CHAN_B] << std::endl;
                std::cout << curCovInv[CHAN_FP] << std::endl;
                std::cout << K<< std::endl;
                std::cout << L << std::endl;
                std::cout << A << std::endl;
            }
            assert(A(x,y) == A(x,y));
        }

//#endif


    const Pixel I_calc = I(i,j);
    const Pixel Ip_calc = Ip(i,j);

    Eigen::MatrixXf b(9,1);
    b.block<3,1>(0,0) = curCovInv[CHAN_F] * curMean[CHAN_F] +
                        I_calc * Alpha(i,j) / SigmaI2;
    b.block<3,1>(3,0) = curCovInv[CHAN_B] * curMean[CHAN_B] +
                        I_calc * (1.0f - Alpha(i,j)) / SigmaI2;;
    b.block<3,1>(6,0) = curCovInv[CHAN_FP] * curMean[CHAN_FP] +
                        Ip_calc * Alpha(i,j) / SigmaI2;


#ifdef PRINT_DEBUG
    std::cout << "Matrix 0 0" << std::endl << M[0][0] << std::endl << std::endl;
    std::cout << "Matrix 0 1" << std::endl << M[0][1] << std::endl << std::endl;
    std::cout << "Matrix 0 2" << std::endl << M[0][2] << std::endl << std::endl;
    std::cout << "Matrix 1 0" << std::endl << M[1][0] << std::endl << std::endl;
    std::cout << "Matrix 1 1" << std::endl << M[1][1] << std::endl << std::endl;
    std::cout << "Matrix 1 2" << std::endl << M[1][2] << std::endl << std::endl;
    std::cout << "Matrix 2 0" << std::endl << M[2][0] << std::endl << std::endl;
    std::cout << "Matrix 2 1" << std::endl << M[2][1] << std::endl << std::endl;
    std::cout << "Matrix 2 2" << std::endl << M[2][2] << std::endl << std::endl;

    std::cout << "Matrix A" << std::endl << A << std::endl << std::endl;
    std::cout << "Vector b" << std::endl << b << std::endl << std::endl;
#endif

    Eigen::MatrixXf x(9,1);
    x = A.fullPivLu().solve(b);

    F(i,j) = x.block<3,1>(0,0);
    B(i,j) = x.block<3,1>(3,0);
    Fp(i,j) = x.block<3,1>(6,0);

    if (F(i,j)(0) < minF(0)) minF(0) = F(i,j)(0);
    if (F(i,j)(1) < minF(1)) minF(1) = F(i,j)(1);
    if (F(i,j)(2) < minF(2)) minF(2) = F(i,j)(2);

    if (F(i,j)(0) > maxF(0)) maxF(0) = F(i,j)(0);
    if (F(i,j)(1) > maxF(1)) maxF(1) = F(i,j)(1);
    if (F(i,j)(2) > maxF(2)) maxF(2) = F(i,j)(2);

    if (B(i,j)(0) < minB(0)) minB(0) = B(i,j)(0);
    if (B(i,j)(1) < minB(1)) minB(1) = B(i,j)(1);
    if (B(i,j)(2) < minB(2)) minB(2) = B(i,j)(2);

    if (B(i,j)(0) > maxB(0)) maxB(0) = B(i,j)(0);
    if (B(i,j)(1) > maxB(1)) maxB(1) = B(i,j)(1);
    if (B(i,j)(2) > maxB(2)) maxB(2) = B(i,j)(2);

    if (Fp(i,j)(0) < minFP(0)) minFP(0) = Fp(i,j)(0);
    if (Fp(i,j)(1) < minFP(1)) minFP(1) = Fp(i,j)(1);
    if (Fp(i,j)(2) < minFP(2)) minFP(2) = Fp(i,j)(2);

    if (Fp(i,j)(0) > maxFP(0)) maxFP(0) = Fp(i,j)(0);
    if (Fp(i,j)(1) > maxFP(1)) maxFP(1) = Fp(i,j)(1);
    if (Fp(i,j)(2) > maxFP(2)) maxFP(2) = Fp(i,j)(2);
}



