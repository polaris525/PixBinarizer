
#ifndef PixAdaptiveBinarizer_h
#define PixAdaptiveBinarizer_h

#include "allheaders.h"

class PixAdaptiveBinarizer
{
public:
    bool mDebug;
    PixAdaptiveBinarizer(bool debug);
    Pix* bradleyAdaptiveThresholding(Pix*, float, int);
    void bradleyAdaptiveThresholdingInverse(Pix *, Pix *, float, Pix **, Pix **, Pix **);
    Pix* pixWindowedMeanMasked(Pix*, Pix*, int);

};

#endif /* PixAdaptiveBinarizer_h */
