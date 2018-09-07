#include "PixAdaptiveBinarizer.h"

PixAdaptiveBinarizer::PixAdaptiveBinarizer(bool debug)
{
    mDebug = debug;
}

Pix* PixAdaptiveBinarizer::bradleyAdaptiveThresholding(Pix* pixArg0, l_float32 fltArg4, int nArg8)  //Complete
{
    L_TIMER timer47 = startTimerNested();
    l_int32 height = pixGetHeight(pixArg0);
    l_int32 width = pixGetWidth(pixArg0);
    Pix* pix48 = pixCreate(width, height, 1);
    l_int32 n08 = height / 2 - 1;
    l_int32 n09 = width / 2 - 1;
    
    if (n09 > nArg8)
        n09 = nArg8;
    if (n09 >= n08)
        n09 = n08;
    L_TIMER timer10 = startTimerNested();
    Pix* pix66 = pixAddMirroredBorder(pixArg0, n09 + 1, n09 + 1, n09 + 1, n09 + 1);
    Pix* pix11 = pixWindowedMean(pix66, n09, n09, 1, 1);
    pixDestroy(&pix66);

    printTimer(timer10, "\tpixWindowMean");
    l_int32 n62, n61, n60, n59;
    pixGetDimensions(pix11, &n62, &n61, 0);
    pixGetDimensions(pix48, &n60, &n59, 0);

    l_float32 flt16 = 1.0 - fltArg4;

    l_uint32* data18 = pixGetData(pix11);
    l_uint32* data19 = pixGetData(pix48);
    l_uint32* data20 = pixGetData(pixArg0);

    l_float32* pFlt22 = new l_float32[0x100];

    for (int i = 0; i < 0x100; i++)
        pFlt22[i] = i * flt16;
    
    l_uint32* line20, line18;
    for (int i = 0; i < height; i++)
    {
        line20 = data20 + i * pixGetWpl(pix48);
        line18 = data18 + i * pixGetWpl(pix11);
        l_int32 n35 = 0;
        for (int j = 0; j < width; j++)
        {
            l_uint32 v41 = GET_DATA_BYTE(line20, j);
            l_uint32 v42 = GET_DATA_BYTE(line18, j);
            if (pFlt22[v42] != v41)
                n35 |= 1 << (j % 32);
        }

        *data19 = n35;
        data19++;
    }

    printTimer(timer47, "bradley_threshold");
    if (mDebug)
        pixWrite("pixMean.png", pix11, 3);
    pixCopyResolution(pix48, pixArg0);
    pixDestroy(pix11);

    return pix48;
}

void PixAdaptiveBinarizer::bradleyAdaptiveThresholdingInverse(Pix *pixArg0, Pix *pixArg4, float fltArg8, Pix **ppix_ArgC, Pix **ppix_Arg10, Pix **ppix_Arg14)   //Complete
{
    L_TIMER timer55 = startTimerNested();
    l_int32 height = pixGetHeight(pixArg0);
    l_int32 width = pixGetWidth(pixArg0);
    *ppix_ArgC = pixCreate(width, height, 1);
    *ppix_Arg10 = pixCreate(width, height, 1);
    *ppix_Arg14 = pixCreate(width, height, 8);
    l_int32 n12 = height / 2 - 1;
    l_int32 n11 = width / 2 - 1;
    
    if (n12 < n11)
        n11 = n12;
    if (n11 >= 255)
        n11 = 255
    L_TIMER timer13 = startTimerNested();
    Pix* pix76 = pixAddMirroredBorder(pixArg0, n11 + 1, n11 + 1, n11 + 1, n11 + 1);
    Pix* pix15 = pixWindowedMeanMasked(pix66, n09, n09, 1, 1);
    pixDestroy(&pix76);

    printTimer(timer13, "\tpixWindowedMeanMasked");

    l_uint32* data22 = pixGetData(pix15);
    l_uint32* data69 = pixGetData(*ppix_ArgC);
    l_uint32* data23 = pixGetData(pixArg0);

    l_float32 flt16 = 1.0 - fltArg4;
    l_uint32* data67 = pixGetData(*ppix_Arg10);

    l_float32* pFlt25 = new l_float32[0x100];

    for (int i = 0; i < 0x100; i++)
        pFlt25[i] = i * flt16;
    
    l_uint32* data28 = pixGetData(*ppix_Arg14);

    l_float32 flt00 = flt16 * 255;
    for (int i = 0; i < height; i++)
    {
        l_uint32* line23 = data23 + i * pixGetWpl(pixArg0);
        l_uint32* line22 = data22 + i * pixGetWpl(pix15);
        l_uint32* line69 = data69 + i * pixGetWpl(*ppix_ArgC);
        l_uint32* line67 = data67 + i * pixGetWpl(*ppix_Arg10);
        l_uint32* line28 = data28 + i * pixGetWpl(*ppix_Arg14);
        l_uint8 n33 = 0, n34 = 0;
        for (int j = 0; j < width; j++)
        {
            l_uint8 r04 = GET_DATA_BYTE(line23, j);
            l_uint8 n39 = GET_DATA_BYTE(line22, j);
            l_uint8 n40 = n39 - r04;
            if (n40 < 0) n40 = -n40;
            SET_DATA_BYTE(line28, j, n40);
            
            if (flt00 - pFlt25[n39] > r04 ^ 0xFF)
                n33 |= 1 << (j % 32);

            if (pFlt25[n39] > r04 )
                n34 |= 1 << (j % 32);
        }
        *data69 = n34;
        *data67 = n33;
        data69++;
        data67++;
    }

    printTimer(timer55, "bradley threshold inv");
    if (mDebug)
    {
        pixWrite("pixDiff.png", *ppix_Arg14, 3);
        pixWrite("pixMeanMasked.png", pix15, 3);
    }
    pixDestroy(&pix15);
}

Pix* PixAdaptiveBinarizer::pixWindowedMeanMasked(Pix* pix_Arg0, Pix* pix_Arg4, int nArg8)   //Complete
{
    if (!pix_Arg0)
    {
        if (LeptMsgSeverity > 5)
            return NULL;

        return (PIX *)returnErrorPtr("pixs not defined", "pixWindowedMean", NULL);
    }

    if (pixGetDepth(pix_Arg0) != 8)
    {
        if (LeptMsgSeverity > 5)
            return NULL;

        return returnErrorPtr("pixs not 8 or 32 bpp", "pixWindowedMean", NULL);
    }

    Pix *pix07 = pixClone(pix_Arg0);
    l_int32 w, h;
    pixGetDimensions(pix07, &w, &h, NULL);

    l_int32 wd = w - 2 * (nArg8 + 1);
    l_int32 hd = h - 2 * (nArg8 + 1);

    if (wd < 2 || hd < 2)
    {
        if (LeptMsgSeverity > 5)
            return NULL;
        
        return returnErrorPtr("w or h too small for kernel", "pixWindowedMean", NULL);
    }

    Pix *pix18 = pixCreate(wd, hd, 8);
    if (!pix18)
    {
        if (LeptMsgSeverity > 5)
            return NULL;

        return returnErrorPtr("pixd not made", "pixWindowedMean", NULL);
    }

    Pix *pix19 = pixBlockconvAccum(pix07);
    if (!pix19)
    {
        pixDestroy(&pix07);
        pixDestroy(&pix18);
        if (LeptMsgSeverity > 5)
            return NULL;

        return returnErrorPtr("pixc not made", "pixWindowedMean", NULL);
    }

    l_int32 n22 = pixGetWpl(pix19);
    l_int32 n44 = pixGetWpl(pix18);
    l_int32 n43 = pixGetWpl(pix_Arg4);
    l_uint32 *pData23 = pixGetData(pix18);
    l_uint32 *pData24 = pixGetData(pix19);

    l_int32 wincr = 2 * nArg8 + 1, hincr = 2 * nArg8 + 1;
    l_float32 flt18 = 1.0 / ((l_float32)(wincr) * hincr);
    l_uint32 *pData28 = pixGetData(pix_Arg4);
    for (int i = 0; i < hd; i++)
    {
        l_int8 n31 = 0;
        l_uint32 *linec1 = pData24 + i * n22;
        l_uint32 *linec2 = pData24 + (i + hincr) * n22;
        for (int j = 0; j < wd; j++)
        {
            l_int8 n34 = GET_DATA_BYTE(pData28 + i * n43, j);
            if (n31 != n34)
                n31 = n34;

            l_int8 r01 = linec2[j + wincr + n34] - linec2[j] - linec1[j + wincr + n34] + linec1[j];

            SET_DATA_BYTE(pData23 + i * n44, j, flt18 * r01);
        }
    }

    pixDestroy(&pix19);
    pixDestroy(&pix07);

    return pix18;
}