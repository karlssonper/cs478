/*
 * FlashMatting.h
 *
 *  Created on: Mar 12, 2012
 *      Author: per
 */

#ifndef FLASHMATTING_H_
#define FLASHMATTING_H_

#include "Eigen/Core"
#include <assert.h>
#include <queue>
#include <set>
#include <vector>
#include <iostream>
#include <Eigen/Dense>


typedef Eigen::Vector3f Pixel;
typedef Eigen::Matrix3f Matrix3;

class FlashMatting
{
public:
    FlashMatting(const std::vector<unsigned char> & flash,
                 const std::vector<unsigned char> & noFlash,
                 int width,
                 int height);

    enum TrimapValue {
        UNKNOWN = 0,
        FOREGROUND = 1,
        BACKGROUND = 2,
        UPDATED_UNKNOWN =3
    };
    void generateTriMap();
    void median();
    void loadTriMap(const std::vector<TrimapValue> &trimap);
    bool updateOneLayer();

    std::vector<unsigned char> alphaMap();
private:
    enum MatrixOrder { ROW_MAJOR, COLUMN_MAJOR};

    //Class to handle 2D data with (i,j) operators
    template<class T, MatrixOrder order>
    class Array
    {
    public:
        Array() { };
        Array(const std::vector<T> & in,
                       const int  width,
                       const int height) :
                       array(in), width(width), height(height)
        {
            assert(array.size() == width*height);
        };

        const T & operator()(const int i, const int j) const
        {
            assert(i>=0 && j>=0 && i<width && j < height);
            if (order == ROW_MAJOR) {
                return array[i + width*j];
            } else {
                return array[j + height*i];
            }
        }

        T& operator()(const int i, const int j)
        {
            assert(i>=0 && j>=0 && i<width && j < height);
            if (order == ROW_MAJOR) {
                return array[i + width*j];
            } else {
                return array[j + height*i];
            }
        }

        void set(const std::vector<T> & in, int w, int h)
        {
            array = in;
            height = h;
            width = w;
        }

        void reserve(int w, int h)
        {
            width = w;
            height = h;
            array.reserve(w*h);
        }
    private:
        std::vector<T> array;
        int height;
        int width;
    };

    typedef Array<Pixel, ROW_MAJOR>        PixelArray;
    typedef Array<float, ROW_MAJOR>        AlphaArray;
    typedef Array<TrimapValue, ROW_MAJOR> TrimapArray;
    PixelArray If;
    PixelArray I;
    PixelArray Ip;

    PixelArray F;
    PixelArray B;
    PixelArray Fp;
    AlphaArray Alpha;


    float SigmaI2;
    float SigmaIp2;
    float Epsilon;

    Matrix3 curCovInv[3];
    Pixel curMean[3];

    TrimapArray Trimap;
    int numUnknowns;
    int SamplesPerUnknown;
    int IterationsPerUnknown;
    int W;
    int H;
    int N;

    float avgMinAlpha;
    float avgMaxAlpha;
    int nPosAlpha;
    int nNegAlpha;
    float minAlpha;
    float maxAlpha;
    Pixel minF;
    Pixel minB;
    Pixel minFP;
    Pixel maxF;
    Pixel maxB;
    Pixel maxFP;

    typedef std::pair<int, int> PixelPair;
    class Compare{
    public:
        Compare(int width) : w(width){};
        bool operator() (const PixelPair &lhs, const PixelPair &rhs) const {
            return lhs.first + w * lhs.second < rhs.first + w * rhs.second;
        }
    private:
        int w;
    };
    typedef std::set<PixelPair, Compare> PixelPairSet;
    PixelPairSet unknowns;
    typedef std::queue<PixelPair> Queue;
    Queue queue;

    bool IsUnknown(const int i, const int j) const {
        return Trimap(i,j) == UNKNOWN;
    };

    bool IsForeground(const int i, const int j) const {
        return Trimap(i,j) == FOREGROUND;
    };

    bool IsBackground(const int i, const int j) const {
        return Trimap(i,j) == BACKGROUND;
    };

    int Neighbors(const int i, const int j) const;
    int NeighborsCross(const int i, const int j) const;

    enum Channel { CHAN_F = 0, CHAN_B = 1, CHAN_FP = 2, NUM_CHAN = 3};
    template <Channel T_CHANNEL>
    float Weight(const int i, const int x, const int j, const int z);

    template <Channel T_CHANNEL>
    Pixel WeightedMeanValue(const std::vector<PixelPair> & neighbors,
                            const int i, const int j);

    template <Channel T_CHANNEL>
    void UpdateWeightedMeanAndCovariance(
                                       const std::vector<PixelPair> & neighbors,
                                       const int i,
                                       const int j);

    template <TrimapValue T_CHANNEL>
    void FindNeighborPoints(std::vector<PixelPair> & neighbors,
                            const int i, const int j);

    template <TrimapValue T_CHANNEL>
    void AddNeighbor(std::vector<PixelPair> & neighbors,
                                   const int i, const int j);

    void InitializeAlpha(const int i, const int j);
    void UpdateAlpha(const int i, const int j);
    void UpdateChannels(const int i, const int j);
    void UpdateColors(const int i, const int j);
    void ClampPixel(Pixel&p)
    {
        if (p(0) < 0.0) p(0) = 0.0;
        if (p(1) < 0.0) p(1) = 0.0;
        if (p(2) < 0.0) p(2) = 0.0;
        if (p(0) > 255.0) p(0) = 255.0;
        if (p(1) > 255.0) p(1) = 255.0;
        if (p(2) > 255.0) p(2) = 255.0;
    }
    void ClampAlpha(float & alpha)
    {
        //assert(alpha == alpha);
        if (alpha < 0){
            avgMinAlpha += alpha;
            ++nNegAlpha;
        }
        else{
            ++nPosAlpha;
            avgMaxAlpha += alpha;
        }

        if (alpha > maxAlpha)
            maxAlpha = alpha;
        if (alpha < minAlpha)
            minAlpha = alpha;

        if (alpha < 0.0) alpha = 0.0;
        if (alpha > 1.0) alpha = 1.0;
    }
};

template <FlashMatting::Channel T_CHANNEL>
Pixel FlashMatting::WeightedMeanValue(
        const std::vector<PixelPair> & neighbors, const int i, const int j)
{
    Pixel mean = Pixel(0.0f, 0.0f, 0.0f);
    float weight = 0.0f;
    float sumWeights = 0.0f;
    for (int idx = 0; idx < neighbors.size(); ++idx) {
        int x = neighbors[idx].first;
        int y = neighbors[idx].second;
        weight = Weight<T_CHANNEL>(i,x,j,y);
        sumWeights += weight;
        if (T_CHANNEL == CHAN_F) {
            mean += F(x,y) * weight;
        } else if (T_CHANNEL == CHAN_B) {
            mean += B(x,y) * weight;
        } else if (T_CHANNEL == CHAN_FP) {
            mean += Fp(x,y)*weight;
        }
    }
    if (sumWeights < 0.00001f) {
        return Pixel(0.0f, 0.0f, 0.0f);
    }
    return mean / sumWeights;
}

template <FlashMatting::Channel T_CHANNEL>
float FlashMatting::Weight(const int i, const int x, const int j, const int y)
{
    float g = exp(-((i-x) * (i-x) + (j-y) + (j - y))/(2.0f*8.0f*8.0f));
    float w = g;
    if (Trimap(i,j) == UPDATED_UNKNOWN) {
        if (T_CHANNEL == FlashMatting::CHAN_F ||
            T_CHANNEL == FlashMatting::CHAN_FP) {
            w*= Alpha(x,y)*Alpha(x,y);
        } else {
            w*= (1.0f - Alpha(x,y))*(1.0f - Alpha(x,y));
        }
    }
    return w;
}

template <FlashMatting::Channel T_CHANNEL>
void FlashMatting::UpdateWeightedMeanAndCovariance(
        const std::vector<PixelPair> & neighbors, const int i, const int j)
{
    Matrix3 cov = Matrix3::Identity();
    float weight = 0.0f;
    float sumWeights = 0.0f;
    curMean[T_CHANNEL] = WeightedMeanValue<T_CHANNEL>(neighbors, i,j);
    for (int x = 0; x < 3; ++x)
            assert(curMean[T_CHANNEL](x) == curMean[T_CHANNEL](x));
    for (int idx = 0; idx < neighbors.size(); ++idx) {
        int x = neighbors[idx].first;
        int y = neighbors[idx].second;
        weight = Weight<T_CHANNEL>(i,x,j,y);
        sumWeights += weight;
        Pixel temp;
        if (       T_CHANNEL == CHAN_F)  {
            temp = F(x,y) - curMean[T_CHANNEL];
        } else if (T_CHANNEL == CHAN_B)  {
            temp = B(x,y) - curMean[T_CHANNEL];
            //std::cout << std::endl << temp << std::endl;
        } else if (T_CHANNEL == CHAN_FP) {
            temp = Fp(x,y) - curMean[T_CHANNEL];
        }
        cov += weight * (temp * temp.transpose());
    }


    Matrix3 A = cov/sumWeights;

    Matrix3 B = Matrix3::Identity();
    Matrix3 X = A.fullPivLu().solve(B);
    for (int x = 0; x < 9; ++x)
            assert(X(x) == X(x));
    curCovInv[T_CHANNEL] = X;
}

template <FlashMatting::TrimapValue T_CHANNEL>
void FlashMatting::FindNeighborPoints(std::vector<PixelPair> & neighbors,
                                      const int i, const int j)
{
    int delta = 1;
    bool negx = true, negy = true, posx = true, posy = true;
    int maxnegx = 0, maxposx = 0, maxnegy = 0, maxposy = 0;
    while (neighbors.size() < SamplesPerUnknown && (negx || negy || posx || posy)) {
        //Can we travel in this direction?
        if (negx && i - delta < 0    ) negx = false;
        else if (negx) maxnegx = i - delta;
        if (posx && i + delta > W - 1) posx = false;
        else if (posx) maxposx = i + delta;
        if (negy && j - delta < 0    ) negy = false;
        else if (negy) maxnegy = j - delta;
        if (posy && j + delta > H - 1) posy = false;
        else if (posy) maxposy = j + delta;

        /*
        X = already updated
        0 = non-corner points
        C = corner points

        C O O O C
        O X X X O
        O X X X O
        0 X X X 0
        C 0 0 0 C
        */

        using std::min;
        using std::max;
        //Fill non-corner points
        if (negx) {
            for (int y = maxnegy; y < maxposy; ++y) {
                AddNeighbor<T_CHANNEL>(neighbors, i - delta, y);
            }
        }

        if (posx) {
            for (int y = maxnegy; y < maxposy; ++y) {
                AddNeighbor<T_CHANNEL>(neighbors, i + delta, y);
            }
        }

        if (negy) {
            for (int x = maxnegx; x < maxposx; ++x) {
                AddNeighbor<T_CHANNEL>(neighbors, x, j - delta);
            }
        }

        if (posy) {
            for (int x = maxnegx; x < maxposx; ++x) {
                AddNeighbor<T_CHANNEL>(neighbors, x, j + delta);
            }
        }

        //Corner points
        if (negx && negy && maxnegx == i - delta && maxnegy == j - delta) {
            AddNeighbor<T_CHANNEL>(neighbors, i - delta, j - delta);
        }

        if (negx && posy && maxnegx == i - delta && maxposy == j + delta) {
            AddNeighbor<T_CHANNEL>(neighbors, i - delta, j + delta);
        }

        if (posx && negy && maxposx == i + delta && maxnegy == j - delta) {
            AddNeighbor<T_CHANNEL>(neighbors, i + delta, j - delta);
        }

        if (posx && posy && maxposx == i + delta && maxposy == j + delta) {
            AddNeighbor<T_CHANNEL>(neighbors, i + delta, j + delta);
        }

        delta++;
    }
    neighbors.resize(SamplesPerUnknown);
}

template <FlashMatting::TrimapValue T_CHANNEL>
void FlashMatting::AddNeighbor(std::vector<PixelPair> & neighbors,
                               const int i, const int j)
{
    bool add = false;
    if (T_CHANNEL == FOREGROUND) {
        if (Trimap(i,j) == FOREGROUND || Trimap(i,j) == UPDATED_UNKNOWN) {
            add = true;
        }
    } else if (T_CHANNEL == BACKGROUND) {
        if (Trimap(i,j) == BACKGROUND || Trimap(i,j) == UPDATED_UNKNOWN) {
            add = true;
        }
    }

    if (add)
        neighbors.push_back(PixelPair(i,j));
}
#endif /* FLASHMATTING_H_ */
