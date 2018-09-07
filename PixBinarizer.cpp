/*
 * PixBinarizer.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: renard
 */

#include "PixBinarizer.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include "image_processing.h"

using namespace std;

PixBinarizer::PixBinarizer(bool debug) {
    mDebug= debug;
}


int PixBinarizer::determineThresholdForTile(Pix* pixt, bool debug) {
    l_int32 start = 0, end = 255, i, error, closeSize = 0;
    l_float32 sum, moment, var, y, variance, mean, meanY, countPixels;
    
    NUMA* norm;
    NUMA* histo = pixGetGrayHistogram(pixt, 1);
    numaSetValue(histo, 255, 0); //ignore white pixels
    error = numaGetNonzeroRange(histo, 0, &start, &end);
    if (end == start || error == 1) {
        numaDestroy(&histo);
        return 0;
    }
    closeSize = end - start;
    if (closeSize % 2 == 0) {
        closeSize++;
    }
    norm = numaNormalizeHistogram(histo, 1);
    /*
     l_float32 median;
     numaGetMedian(histo,&median);
     printf("%2.1f\t", median);
     */
    
    l_float32 iMulty;
    for (sum = 0.0, moment = 0.0, var = 0.0, countPixels = 0, i = start;
         i < end; i++) {
        numaGetFValue(norm, i, &y);
        sum += y;
        iMulty = i * y;
        moment += iMulty;
        var += i * iMulty;
        numaGetFValue(histo, i, &y);
        countPixels += y;
    }
    variance = sqrt(var / sum - moment * moment / (sum * sum));
    mean = moment / sum;
    meanY = sum / (end - start);
    int result = 0;
    if (variance < 8) {
        result = 0;
    } else {
        if ((mean * 0.5) < variance) {
            result = mean;// + variance*0.4;
        } else {
            //high variance means we have probably have a very good separation between fore and background
            //in that case the threshold be closer to the mean
            result = mean - sqrt(variance) * 3;
            //result = mean - variance*0.66;
        }
        //result = mean - variance*0.66; //0.2 when median is activated, 0.66 otherwise
    }
    
    if (debug == true) {
        printf("mean = %f , variance = %f, meanY = %f \n", mean, variance,
               meanY);
        printf("%i ,%i\n", start, end);
        GPLOT *gplot;
        ostringstream name;
        name << mean;
        gplot = gplotCreate(name.str().c_str(), GPLOT_X11, name.str().c_str(),
                            "x", "y");
        ostringstream title;
        title << "mean = " << mean << ", " << "variance = " << variance
        << ", thresh = " << result;
        //gplotAddPlot(gplot, NULL, diff, GPLOT_LINES, "diff histogram");
        gplotAddPlot(gplot, NULL, norm, GPLOT_LINES, title.str().c_str());
        gplotMakeOutput(gplot);
        gplotDestroy(&gplot);
    }
    numaDestroy(&norm);
    numaDestroy(&histo);
    
    //TODO check std dev of y values. large value > foreground and background, low value > only background
    //idea: narrow histogramm means a high likelyhood for bg
    return result;
    
}

Pix* PixBinarizer::createEdgeMask(Pix* pixs) {
    L_TIMER timer = startTimerNested();
    ostringstream s;
    Pix* pixConv = pixBlockconvGray(pixs, NULL, 5, 5);
    Pix* pixConvEdges = pixSobelEdgeFilter(pixConv, L_ALL_EDGES);
    pixDestroy(&pixConv);
    pixInvert(pixConvEdges, pixConvEdges);
    s << "sobel edge detection: " << stopTimerNested(timer) << std::endl;
    timer = startTimerNested();
    if(mDebug){
        pixWrite("pixConvEdges.bmp", pixConvEdges, IFF_BMP);
    }
    
    NUMA* histo = pixGetGrayHistogram(pixConvEdges, 8);
    NUMA* norm = numaNormalizeHistogram(histo, 1.0);
    l_float32 median, mean, variance;
    numaGetHistogramStats(norm, 0, 1, &mean, &median, NULL, &variance);
    numaDestroy(&histo);
    numaDestroy(&norm);
    
    l_int32 thresh = 0;
    if (variance < 1.0) {
        thresh = 255;
    } else {
        thresh = 254;
    }
    if(mDebug){
        printf("mean = %f, median = %f, std dev = %f, thresh = %i\n", mean, median,
               variance, thresh);
    }
    
    Pix* pixForeground = pixThresholdToBinary(pixConvEdges, thresh);
    pixDestroy(&pixConvEdges);
    pixCloseSafeBrick(pixForeground,pixForeground,64,64);
    
    pixWrite("foreground.bmp", pixForeground, IFF_BMP);
    
    if(mDebug){
        s << "binarization of edge mask: " << stopTimerNested(timer) << std::endl;
    }
    timer = startTimerNested();
    
    Pix* pixacc = pixBlockconvAccum(pixForeground);
    Pix* pixRank = pixBlockrank(pixForeground, pixacc, 8, 8, 0.1);
    
    pixDestroy(&pixacc);
    pixDestroy(&pixForeground);
    pixInvert(pixRank, pixRank);
    Pix* pixResult = pixOpenBrick(NULL, pixRank, 10, 10);
    pixDestroy(&pixRank);
    if(mDebug){
        pixWrite("rank.bmp", pixResult, IFF_BMP);
        s << "mask generation: " << stopTimerNested(timer) << std::endl;
        printf("%s", s.str().c_str());
    }
    
    return pixResult;
}

/**
 * determines and applies a threshold for each tile separately
 */
Pix* PixBinarizer::binarizeTiled(Pix* pixs, const l_uint32 tileSize) {
    L_TIMER timer = startTimerNested();
    l_int32 thresh, w, h;
    Pix* pixb;
    ostringstream s;
    pixGetDimensions(pixs, &w, &h, NULL);
    l_int32 nx = L_MAX(1, w / tileSize);
    l_int32 ny = L_MAX(1, h / tileSize);
    l_int32 ox = L_MAX(1,nx/6);
    l_int32 oy = L_MAX(1,ny/6);
    PIXTILING* pt = pixTilingCreate(pixs, nx, ny, 0, 0, ox, oy);
    Pix* pixth = pixCreate(nx, ny, 8);
    Pix* pixt;
    bool debug = false;
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            pixt = pixTilingGetTile(pt, i, j);
            if (i == 0 && j == 0) {
                debug = false;
            } else {
                debug = false;
            }
            thresh = determineThresholdForTile(pixt, debug);
            pixSetPixel(pixth, j, i, thresh);
            pixDestroy(&pixt);
        }
    }
    pixTilingDestroy(&pt);
    if(mDebug){
        s << "local threshhold determination: " << stopTimerNested(timer)<< std::endl;
        timer = startTimerNested();
    }
    
    pt = pixTilingCreate(pixs, nx, ny, 0, 0, 0, 0);
    pixb = pixCreate(w, h, 1);
    l_uint32 val;
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            pixt = pixTilingGetTile(pt, i, j);
            pixGetPixel(pixth, j, i, &val);
            Pix* pixbTile = pixThresholdToBinary(pixt, val);
            pixTilingPaintTile(pixb, i, j, pixbTile, pt);
            pixDestroy(&pixt);
            pixDestroy(&pixbTile);
        }
    }
    pixTilingDestroy(&pt);
    pixDestroy(&pixth);
    if(mDebug){
        s << "local threshhold application: " << stopTimerNested(timer)<< std::endl;
        printf("%s", s.str().c_str());
    }
    
    return pixb;
}

void PixBinarizer::binarizeInternal(Pix* pixGrey, Pix* pixhm, Pix** pixb) {
    Pix* pixEdgeMask;
    l_int32 width = pixGetWidth(pixGrey);
    const l_uint32 tileSize = width/10; //size of tile during threshholding
    L_TIMER timer = startTimerNested();
    ostringstream s;
    
    pixEdgeMask = createEdgeMask(pixGrey);
    
    pixSetMasked(pixGrey, pixEdgeMask, 255);
    if (pixhm != NULL) {
        //dont allow image mask to cover text
        pixAnd(pixhm, pixhm, pixEdgeMask);
    }
    
    pixDestroy(&pixEdgeMask);
    if(mDebug){
        s << "text mask creation: " << stopTimerNested(timer) << std::endl;
        printf("%s", s.str().c_str());
        pixDisplay(pixGrey,0,0);
        pixWrite("toThresh.bmp", pixGrey, IFF_BMP);
    }
    
    /*
     timer = startTimerNested();
     Pix* pixMedian = pixRankFilterGray(pixGrey,3,3,0.25);
     printf("rank filter = %f\n",stopTimerNested(timer));
     pixWrite("rank.bmp",pixMedian,IFF_BMP);
     */
    
    *pixb = binarizeTiled(pixGrey, tileSize);
    //pixDestroy(&pixMedian);
    
    if (pixhm != NULL) {
        pixSetMasked(*pixb, pixhm, 0);
    }
}

Pix *PixBinarizer::binarize(Pix *pix, void(*previewCallBack)(Pix *)) {
    
    L_TIMER timer06 = startTimerNested();
    l_int32 width, height, depth;
    pixGetDimensions(pix, &width, &height, &depth);
    Pix* pixGrey = NULL;
    Pix* pixtContrast;
    Pix* pixBinary;
    switch(depth){
        case 1:
            return pixClone(pix);
        case 8:
            pixGrey = pixClone(pix);
            break;
        case 32:
            pixGrey = pixConvertTo8(pix, 0);
    }
    m_pix0C = pixClone(pixGrey)
    Pix* pix10 = findGrayScaleEdges(pixGrey);       //done

    calculateBinarisationThresholdParams(pix10);        //done
    Pix* pix11 = binariseGreyscaleEdges(pix10);         
    pixDestroy(pix10);  //done
    /*
    if (previewCallBack)    //NULL
    {

    }*/

    Pix *pix22 = createEdgeSizeMask(pix11);  //done
    Pix *pix38, *pix3C, *pix40;
    m_pixAdaptiveBin08.bradleyAdaptiveThresholdingInverse(pixGrey, pix22, m_flt20, &pix38, &pix3C, &pix40);

    pixDestroy(&pix22);
    pixDestroy(&pixGrey);
    pixBinary = combineBinaryImages(pix38, pix3C, pix11, pix40);

    pixDestroy(&m_pix0C);
    pixDestroy(&pix40);

    printTimer(timer06, "binarise");
    pixDestroy(&pixGrey);
    
    return pixBinary;
}

Pix *PixBinarizer::findGrayScaleEdges(Pix *pix) //Complete
{
    Pix* pixRet = smoothGreyPix(pix);
    PixEdgeDetector edgeDetector;
    Pix* pix03 = edgeDetector.makeEdges(pixRet);
    pixDestroy(pixRet);

    if (mDebug)
        pixWrite("grayEdges.png", pix03, 3);

    return pix03;
}

Pix *PixBinarizer::smoothGreyPix(Pix *pix)  //Complete
{
    L_TIMER timer_04 = startTimerNested();
    Pix* pix05 = pixBlockconvGray(pix, NULL, 2, 2);
    printTimer(timer_04, "smoothGreyPix");
    return pix05;
}

void PixBinarizer::calculateBinarisationThresholdParams(Pix *pix)   //Complete
{
    L_TIMER timer05 = startTimerNested();
    NUMA* pNuma_06 = pixGetGrayHistogram(pix, 4);
    NUMA* pNuma_38 = numaNormalizeHistogram(pNuma_06, 1.f);
    float thresh_24;
    numaSplitDistribution(pNuma_38, 0, &thresh_24, NULL, NULL, NULL, NULL, NULL);

    l_float32 mean, variance, stdd;
    numaGetHistogramStats(pNuma_06, 0, 1.f, &mean, NULL, NULL, &variance);
    numaDestroy(pNuma_06);
    numaDestroy(pNuma_38);

    stdd = sqrt(variance);
    l_float32 flt04 = (stdd + mean) * 0.5;
    double dbl02 = flt04;
    l_float32 flt_sp24 = dbl02 * 0.00416667 - 0.0333333;

    l_float32 flt_sp20 = (stdd + (float)thresh_24) * 0.0025f - 0.275f;
    l_float32 flt_r03, flt_r04;
    if (flt_sp20 < 0.3)
    {
        flt_r03 = m_flt20;
        if (flt_r03 < 0.05)
            flt_r03 = 0.05;
    }
    m_flt20 = flt_r03;
    if (flt_sp24 < 0.3)
    {
        flt_r04 = m_flt24;
        if (flt_r04 < 0.05)
            flt_r04 = 0.05;
    }
    m_flt24 = flt_r04;
    m_flt1C = 0.05;
    m_flt18 = 0.1;

    float edge_thresh = 1.32499993f, thresh = -1.08420217e-19;
    if (mDebug)
        printf("mean = %.2f, stdd = %.2f, edge thresh = %.2f, thresh = %.2f, avg component thresh = %.2f, max component thresh = %.2f, split=%i\n", 
            mean, stdd, edge_thresh, thresh, m_flt24, m_flt20, (l_int32)thresh_24);

    printTimer(timer05, "measureContrastNearEdges");
}

Pix* PixBinarizer::binariseGreyscaleEdges(Pix *pix) //Complete
{
    pixInvert(pix, pix);
    Pix* pix05 = m_pixAdaptiveBin08.bradleyAdaptiveThresholding(pix, m_flt24, 8);
    if (mDebug)
        pixWrite("binaryEdges.png", pix05, 3);

    return pix05;
}

Pix *PixBinarizer::createEdgeSizeMask(Pix *pix) //Complete
{
    L_TIMER timer04;
    if (mDebug)
        timer04 = startTimerNested();
    else
        timer04 = 0;
    
    l_int32 width, height;
    pixGetDimensions(pix, &width, &height, NULL);
    BOXA* pBoxa_18 = findConnCompSizeByRank(pix, 0.85);
    BOXA* pBoxa_12 = makeSizeMaskFromBoxa(pBoxa_18, width, height, m_edgeH_10);
    boxaDestroy(&pBoxa_18);

    if (mDebug)
    {
        printf("min edge height = %i\n", m_edgeH_10);
        printTimer(timer04, "createEdgeSizeMask");
    }

    return pBoxa_12;
}

BOXA *PixBinarizer::findConnCompSizeByRank(Pix *pix, double dblArg4)    //Complete
{
    L_TIMER timer06 = startTimerNested();
    BOXA* pBoxa_24 = pixConnCompBB(pix, 8);
    NUMA* pNuma_28, pNuma_2C:
    boxaExtractAsNuma(pBoxa_24, NULL, NULL, NULL, &pNuma_2C, &pNuma_28, NULL, 0);
    NUMA *pNuma_30 = numaMakeThresholdIndicator(pNuma_2C, 8.f, 2);
    BOXA* pBoxa_r9 = boxaSelectWithIndicator(pBoxa_24, pNuma_30, NULL);
    l_int32 nVar34, nVar38;
    boxaGetRankVals(pBoxa_r9, (l_float32)dblArg4, &nVar34, &nVar38, NULL, NULL);

    numaDestroy(&pNuma_2C);
    numaDestroy(&pNuma_2C);
    numaDestroy(&pNuma_28);
    numaDestroy(&pNuma_30);
    boxaDestroy(&pBoxa_24);

    printTimer(timer06, "finding edge height");
    if (mDebug)
        printf("edgeHeight = %i\n", m_edgeH_10);

    return pBoxa_r9;
}

NUMA *boxaMakeContainedIndicator(BOXA* pBoxa_Arg0)  //Complete
{
    l_int32 count = boxaGetCount(pBoxa_Arg0);
    NUMA *pNuma27;
    BOXA *pBoxa26 = boxaSort(pBoxa_Arg0, 10, 1, &pNuma27);
    NUMA *pNuma04 = numaMakeConstant(1.f, count);
    NUMAA *pNumaa25 = numaaCreateFull(count, 0);

    l_int32 n22, n24 = 0;

    for (int i = 0; i < count; i++)
    {
        BOX* pBox23 = boxaGetBox(pBoxa26, i, 2);
        BOX *pBox21;
        for (int j = i + 1; j < count; j++)
        {
            BOX *pBox21 = boxaGetBox(pBoxa26, j, 2);
            
            numaGetIValue(pNuma27, j, &n22);

            boxContains(pBox21, pBox23, &n24);
            if (n24 == 1)
            {
                l_int32 n02;
                numaGetIValue(pNuma27, i, &n02);
                numaaAddNumber(pNumaa25, n22, n02);
                boxDestroy(&pBox21);
                break;
            }
                
            boxDestroy(&pBox21);
        }

        boxDestroy(&pBox23);        
    }

    for (int i = 0; i < count; i++)
    {
        l_int32 n15 = numaaGetNumaCount(pNumaa25, i);
        if (n15 == 2)
        {
            numaaGetValue(pNumaa25, i, 0, NULL, NULL);
            numaSetValue(pNuma04, n24, 0);
            numaaGetValue(pNumaa25, i, 1, NULL, NULL);
            numaSetValue(pNuma04, n24, 0);
        }
        else if (n15 == 1)
        {
            numaaGetValue(pNumaa25, i, 0, NULL, NULL);
            numaSetValue(pNuma04, n24, 0);
        }
        else if (n15 >= 3)
            numaSetValue(pNuma04, n24, 0);
    }

    numaaDestroy(&pNumaa25);
    numaDestroy(&pNuma27);
    boxaDestroy(&pBoxa26);
    return pNuma04;
}

Pix *PixBinarizer::makeSizeMaskFromBoxa(BOXA *pBoxa_Arg0, int w, int h, int edge)   //Complete
{
    Pix* pix13 = pixCreate(w, h, 8);
    pixSetAllGray(pix13, edge * 5);
    NUMA *pNuma24 = boxaMakeContainedIndicator(pBoxa_Arg0);
    BOXA *pBoxa28 = boxaSelectWithIndicator(pBoxa_Arg0, pNuma_24, NULL);
    for (int i = 0; i < boxaGetCount(); i++)
    {
        BOX* pBox14 = boxaGetBox(pBoxa28, i, 2);
        boxAdjustSides(pBox14, pBox14, -1, 1, -1, 1);
        l_int32 x, y, w, h;
        boxGetGeometry(pBox14, &x, &y, &w, &h);
        if (w > edge)
            pixSetInRectArbitrary(pix13, pBox14, edge * 5);
        boxDestroy(pBox14);
    }

    if (mDebug)
    {
        renderBoxaForDebugging(m_pix0C, pBoxa28, "edgeSizeBoxa.png");
        pixWrite("sizeAreas.png", pix13, 3);
    }

    numaDestroy(&pNuma24);
    boxaDestroy(&pBoxa28);

    return pix13;
}

void PixBinarizer::renderBoxaForDebugging(Pix* pix_Arg0, BOXA* pBoxa_Arg4, const char* fname);  //Complete
{
    if (mDebug)
    {
        Pix* pix08 = pixConvertTo32(pix_Arg0);
        pixRenderBoxaArb(pix08, pBoxa_Arg4, 1, 255, 0, 0)
        pixWrite(fname, pix08, 3);
        pixDestroy(pix08);
    }
}

Pix *PixBinarizer::combineBinaryImages(Pix* pix_Arg0, Pix* pix_Arg4, Pix* pix_Arg8, Pix* pix_ArgC)  //Complete
{
    l_int32 n28, n27;

    pixGetDimensions(pix_Arg0, &n28, &n27, 0);
    Numaa* pNumaa26, *pNumaa25;
    Pixa* pixa24 = filterTextComponentsCalculateScore(pix_Arg0, pix_Arg8, pix_ArgC, &pNumaa26, "Org");
    pixDestroy(&pix_Arg0);
    Pixa* pixa18 = filterTextComponentsCalculateScore(pix_Arg4, pix_Arg8, pix_ArgC, &pNumaa25, "Inv");
    pixDestroy(&pix_Arg4);
    pixDestroy(&pix_Arg8);

    Pix* pix11 = combineTextComponents(pixa24, pixa18, n28, n27, pNumaa26, pNumma25);
    pixaDestroy(&pixa24);
    pixaDestroy(&pixa18);
    numaaDestroy(&pNumaa25);
    numaaDestroy(&pNumaa26);

    return pix11;
}

NUMA *numaSelectWithIndicator(NUMA *pNuma_Arg0, NUMA *pNuma_Arg4)   //Complete
{
    l_int32 count = numaGetCount(pNuma_Arg0);
    NUMA *pNuma05 = numaCreate(count);

    for (int i = 0; i < count; i++)
    {
        l_int32 n10;
        l_float32 f09;

        numaGetIValue(pNuma_Arg4, i, &n10);
        if (n10)
        {
            numaGetFValue(pNuma_Arg0, i, &f09);
            numaAddNumber(pNuma05, f09);
        }
    }

    return pNuma05;
}

Pix * PixBinarizer::filterTextComponentsCalculateScore(Pix* pix_Arg0, Pix* pix_Arg4, Pix* pix_Arg8, Numaa** ppNumaa_ArgC, string strArg10) //Complete
{
    L_TIMER timer09 = startTimerNested();
    L_TIMER timer10 = startTimerNested();
    l_int32 n59, n58;

    pixGetDimensions(pix_Arg0, &n59, &n58, 0);
    Pixa* pixa57;
    BOXA* pBoxa56 = pixConnCompPixa(pix_Arg0, &pixa57, 8);
    boxaDestroy(&pBoxa56);

    string str54 = "\tpixConnCompPixa" + strArg10;
    printTimer(timer09, str54);

    Pixa* pixa12 = removeNonTextComponentsFromPixa(n59, n58, pixa57);
    pixaDestroy(&pixa57);

    Numa *pNuma52, *pNuma51, *pNuma50;
    calculateComponentScore(pixa12, pix_Arg4, pix_Arg8, &pNuma52, &pNuma51, &pNuma50);
    Numa* pNuma49 = numaMakeThresholdIndicator(pNuma52, m_flt18, 4);
    Numa* pNuma48 = numaMakeThresholdIndicator(pNuma51, m_flt1C, 4);
    Numa* pNuma47 = numaMakeThresholdIndicator(pNuma50, 0.05f, 4);
    Numa* pNuma46 = numaMakeThresholdIndicator(pNuma50, 0.5f, 4);
    numaLogicalOp(pNuma49, pNuma49, pNuma47, 6);
    numaLogicalOp(pNuma49, pNuma49, pNuma48, 6);
    numaLogicalOp(pNuma49, pNuma49, pNuma46, 5);
    Pixa *pixa18 = pixaSelectWithIndicator(pixa12, pNuma49, NULL);
    *ppNumaa_ArgC = numaaCreate(0);
    NUMA *pNuma21 = numaSelectWithIndicator(pNuma51, pNuma49);
    NUMA *pNuma22 = numaSelectWithIndicator(pNuma52, pNuma49);
    NUMA *pNuma23 = numaSelectWithIndicator(pNuma50, pNuma49);
    numaaAddNuma(*ppNumaa_ArgC, pNuma21, 0);
    numaaAddNuma(*ppNumaa_ArgC, pNuma22, 0);
    numaaAddNuma(*ppNumaa_ArgC, pNuma23, 0);

    if (mDebug)
    {
        l_int32 n45, n44;
        pixGetDimensions(pix_Arg0, &n45, &n44, NULL);
        BOXA *pBoxa43 = pixaGetBoxa(pixa12, 2);
        string str26 = strArg10 + "pixBinary.png";
        pixWrite(str26, pix_Arg0, 3);
        string str37 = strArg10 + "filteredTextComponents.png";
        renderBoxaForDebugging(m_pix0C, pBoxa43, str37);
        BOXA *pBoxa36 = pixaGetBoxa(pixa18, 2);
        str37 = strArg10 + "filteredByArearFractionTextComponents.png";
        renderBoxaForDebugging(m_pix0C, pBoxa36, str37);
        boxaDestroy(&pBoxa43);
        boxaDestroy(&pBoxa36);
    }

    pixaDestroy(&pixa12);
    numaDestroy(&pNuma51);
    numaDestroy(&pNuma50);
    numaDestroy(&pNuma52);
    numaDestroy(&pNuma49);
    numaDestroy(&pNuma48);
    numaDestroy(&pNuma46);
    numaDestroy(&pNuma47);

    printTimer(timer10, "filterTextComponentsByEdgeMask");

    return pixa18;
}

BOX * boxScaleCentered(BOX *pBox_Arg0, float fltArg4)   //Complete
{
    l_int32 x, y, w, h;

    boxGetGeometry(pBox_Arg0, &x, &y, &w, &h);
    BOX* pBoxResult = boxCreate(x - (w * fltArg4 - w) / 2, y - (h * fltArg4 - h) / 2, w * fltArg4, h * fltArg4);

    return pBoxResult;
}

bool boxTouches(BOX *pBox_Arg0, BOX *pBox_Arg4) //Complete
{
    bool result = false;
    if (pBox_Arg0 && pBox_Arg4)
    {
        l_int32 n17, n15, n13, n12, n16, n14, n11, n10;
        boxGetGeometry(pBox_Arg0, &n17, &n15, &n13, &n12);
        boxGetGeometry(pBox_Arg4, &n16, &n14, &n11, &n10);

        if (n10 + n14 >= n15)
        {
            bool flag04 = false;
            if (n15 + n12 >= n14)
                flag04 = true;
            
            bool flag08 = false;
            if (n13 + n17 >= n16)
                flag08 = true;
            flag08 = flag04 & flag08;
            if (n11 + n16 >= n17)
                result = true;
            result &= flag08;
        }
    }
    else if (LeptMsgSeverity > 5)
        return true;
    else
        return (returnErrorInt("box1 and box2 not both defined", "boxTouches", 1) != 0);
    
    return result;
}

Pix * PixBinarizer::combineTextComponents(Pixa* pixa_Arg0, Pixa* pixa_Arg4, int nArg8, int nArgC, Numaa* pNumaa_Arg10, Numaa* pNumaa_Arg14) //Complete
{
    L_TIMER timer49 = startTimerNested();
    BOXA* pBoxa10 = pixaGetBoxa(pixa_Arg0, 2);
    BOXA* pBoxa11 = pixaGetBoxa(pixa_Arg4, 2);

    if (mDebug)
    {
        Pix* pix14 = pixConvertTo32(m_pix0C);
        pixRenderBoxaArb(pix14, pBoxa10, 1, 255, 0, 0);
        pixRenderBoxaArb(pix14, pBoxa11, 1, 255, 0, 0);
        pixWrite("allComponents.png", pix14, 3);
        pixDestroy(&pix14);
    }

    L_TIMER timer50 = startTimerNested();
    l_int32 n12 = boxaGetCount(pBoxa10);
    l_int32 n17 = boxaGetCount(pBoxa11);
    NUMA* pNuma72 = numaMakeConstant(1.f, n12);
    NUMA* pNuma71 = numaMakeConstant(1.f, n17);

    for (int i = 0; i < n17; i++)
    {
        BOX *pBox70 = boxaGetBox(pBoxa10, i, 2);
        BOX *pBox69 = boxScaleCentered(pBox70, 1.2f);
        l_int32 w, h;
        boxGetGeometry(pBox69, NULL, NULL, &w, &h);
        l_int32 n67;
        numaGetIValue(pNuma71, i, &n67);
        l_int32 n21 = w * h;
        if (n21 && n67)
        {
            for (int j = 0; j < n12; j++)
            {
                numaGetIValue(pNuma72, j, &n67);
                if (n67)
                {
                    BOX *pBox66 = boxaGetBox(pBoxa11, j, 2);
                    BOX *pBox65 = boxScaleCentered(pBox66, 1.2f);
                    l_int32 area2C;
                    boxOverlapArea(pBox69, pBox65, &area2C);
                    if (boxTouches(pBox69, pBox65) || area2C >= 1)
                    {
                        boxGetGeometry(pBox65, NULL, NULL, &w, &h);
                        l_int32 n23 = w * h;
                        if (n23)
                        {
                            l_float32 flt64, flt63, flt60, flt62, flt61, flt59;

                            numaaGetValue(pNumaa_Arg10, 0, i, &flt64, NULL);
                            numaaGetValue(pNumaa_Arg10, 1, i, &flt63, NULL);
                            numaaGetValue(pNumaa_Arg10, 2, i, &flt60, NULL);
                            numaaGetValue(pNumaa_Arg14, 0, j, &flt62, NULL);
                            numaaGetValue(pNumaa_Arg14, 1, j, &flt61, NULL);
                            numaaGetValue(pNumaa_Arg14, 2, j, &flt59, NULL);

                            BOX *pBox24, *pBox25;
                            if (n21 <= n23)
                            {
                                pBox24 = pBox69;
                                pBox25 = pBox65;
                            }
                            else
                            {
                                pBox24 = pBox65;
                                pBox25 = pBox69;
                            }

                            l_int32 n58;
                            boxContains(pBox25, pBox24, &n58);
                            if (flt64 * flt63 <= flt62 * flt61)

                            Numa *pNuma33; l_int32 n34;
                            if ((flt64 * flt63 <= flt62 * flt61) && (n23 < 7 | n21 > n23 | n21 > 6))
                            {
                                pNuma33 = pNuma72;
                                n34 = j;
                            }
                            else
                            {
                                pNuma33 = pNuma71;
                                n34 = i;
                            }

                            numaSetValue(pNuma33, n34, 0);
                        }
                    }

                    boxDestroy(&pBox65);
                    boxDestroy(&pBox66);
                }
            }
            boxDestroy(&pBox70);
            boxDestroy(&pBox69);
        }
        else
            boxDestroy(&pBox69);
    }

    printTimer(timer50, "handle overlaps");
    Pixa *pixa70 = pixaSelectWithIndicator(pixa_Arg0, pNuma71, NULL);
    Pixa *pixa69 = pixaSelectWithIndicator(pixa_Arg4, pNuma72, NULL);

    if (mDebug)
    {
        Pix *pix38 = pixConvertTo32(m_pix0C);
        BOXA *pBoxa68 = pixaGetBoxa(pixa70, 2);
        BOXA *pBoxa67 = pixaGetBoxa(pixa69, 2);
        
        pixRenderBoxaArb(pix38, pBoxa68, 1, 0, 255, 0);
        pixRenderBoxaArb(pix38, pBoxa67, 1, 255, 0, 0);
        pixWrite("filteredCandites.png", pix38, 3);
        pixDestroy(&pix38);
        boxaDestroy(&pBoxa68);
        boxaDestroy(&pBoxa67);
    }

    PIXAA *pixaa41 = pixaaCreate(2);
    pixaaAddPixa(pixaa41, pixa70, 2);
    pixaaAddPixa(pixaa41, pixa69, 2);
    PIXA *pixa68 = pixaaFlattenToPixa(pixaa41, NULL, 2);
    PIX *pix42 = pixaDisplay(pixa68, nArg8, nArgC);

    printTimer(timer49, "inverse combination");

    boxaDestroy(&pBoxa10);
    boxaDestroy(&pBoxa11);
    pixaDestroy(&pixa70);
    pixaDestroy(&pixa69);
    pixaaDestroy(&pixaa41);
    pixaDestroy(&pixa68);
    numaDestroy(&pNuma72);
    numaDestroy(&pNuma71);

    return pix42;
}

Pixa * PixBinarizer::removeNonTextComponentsFromPixa(int nArg0, int nArg4, Pixa* pixa_Arg8) //Complete
{
    L_TIMER timer28 = startTimerNested();
    NUMA *pNuma46, *pNuma45;
    pixaFindDimensions(pixa_Arg8, &pNuma46, &pNuma45);
    NUMA *pNuma07 = pixaFindAreaFraction(pixa_Arg8);
    NUMA *pNuma08 = pixaFindWidthHeightRatio(pixa_Arg8);

    l_int32 n06 = max(nArg0 / 2, 500);

    NUMA *pNuma42 = numaMakeThresholdIndicator(pNuma45, m_edgeH_10 / 4, 4);
    NUMA *pNuma41 = numaMakeThresholdIndicator(pNuma46, m_edgeH_10 % 4, 4);
    NUMA *pNuma39 = numaMakeThresholdIndicator(pNuma07, 0.15, 4);
    NUMA *pNuma38 = numaMakeThresholdIndicator(pNuma45, n06, 3);
    NUMA *pNuma37 = numaMakeThresholdIndicator(pNuma46, n06, 3);
    NUMA *pNuma40 = numaMakeThresholdIndicator(pNuma07, 0.5, 4);
    NUMA *pNuma36 = numaMakeThresholdIndicator(pNuma08, 0.5, 4);
    NUMA *pNuma35 = numaMakeThresholdIndicator(pNuma08, 2, 3);
    NUMA *pNuma34 = numaMakeThresholdIndicator(pNuma45, m_edgeH_10 / 4, 3);
    NUMA *pNuma33 = numaMakeThresholdIndicator(pNuma46, m_edgeH_10 % 4, 3);

    numaLogicalOp(pNuma36, pNuma36, pNuma35, 6);
    numaLogicalOp(pNuma34, pNuma34, pNuma33, 6);
    numaLogicalOp(pNuma36, pNuma36, pNuma34, 6);
    numaLogicalOp(pNuma36, pNuma36, pNuma40, 6);
    NUMA *pNuma32 = numaLogicalOp(NULL, pNuma42, pNuma41, 5);
    numaLogicalOp(pNuma32, pNuma32, pNuma39, 6);
    numaLogicalOp(pNuma32, pNuma32, pNuma38, 6);
    numaLogicalOp(pNuma32, pNuma32, pNuma37, 6);
    numaLogicalOp(pNuma32, pNuma32, pNuma36, 5);
    Pixa *pixa24 = pixaSelectWithIndicator(pixa_Arg8, pNuma32, 0);

    numaDestroy(&pNuma46);
    numaDestroy(&pNuma45);
    numaDestroy(&pNuma07);
    numaDestroy(&pNuma08);
    numaDestroy(&pNuma32);
    numaDestroy(&pNuma42);
    numaDestroy(&pNuma41);
    numaDestroy(&pNuma40);
    numaDestroy(&pNuma39);
    numaDestroy(&pNuma38);
    numaDestroy(&pNuma37);
    numaDestroy(&pNuma36);
    numaDestroy(&pNuma35);
    numaDestroy(&pNuma34);
    numaDestroy(&pNuma33);

    printTimer(timer28, "\tremoveNonTextComponentsFromPixa");

    return pixa24;
}

Pix * pixClipMaskedWithoutInvert(Pix *pix_Arg0, Pix *pix_Arg4, int nArg8, int nArgC)    //Complete
{
    Pix *pix08 = NULL;
    if (pix_Arg0)
    {
        if (pix_Arg4 && pixGetDepth(pix_Arg4) == 1)
        {
            l_int32 w, h;
            pixGetDimensions(pix_Arg4, &w, &h, 0);
            BOX *pBox12 = boxCreate(nArg8, nArgC, w, h);
            pix08 = pixClipRectangle(pix_Arg0, pBox12, NULL);
            pixRasterop(pix08, 0, 0, w, h, 8, pix_Arg4, 0, 0);
            boxDestroy(&pBox12);

            return pix08;
        }

        if (LeptMsgSeverity <= 5)
            return returnErrorPtr("pixm undefined or not 1 bpp", "pixClipMasked", NULL);
    }
    else if (LeptMsgSeverity <= 5)
        return returnErrorPtr("pixs not defined", "pixClipMasked", NULL);

    return pix08;
}

void PixBinarizer::calculateComponentScore(Pixa* pixa_Arg0, Pix* pix_Arg4, Pix* pix_Arg8, Numa** ppNuma_ArgC, Numa** ppNuma_Arg10, Numa** ppNuma_Arg14) //Complete
{
    L_TIMER timer30 = startTimerNested();
    BOXA *pBoxa45 = pixaGetBoxa(pixa_Arg0, 2);
    *ppNuma_ArgC = numaCreate(0);
    *ppNuma_Arg10 = numaCreate(0);
    *ppNuma_Arg14 = numaCreate(0);
    l_int32 count = boxaGetCount(pBoxa45);

    for (int i = 0; i < count; i++)
    {
        BOX *pBox44 = boxaGetBox(pBoxa45, i, 2);
        Pix *pix43 = pixaGetPix(pixa_Arg0, i, 2);
        l_int32 x, y, w, h;
        boxGetGeometry(pBox44, &x, &y, &w, &h);
        Pix *pix38 = pixClipMasked(pix_Arg4, pix43, x, y, 0);
        l_uint32 *pData14 = pixGetData(pix38);
        l_int32 n15 = pixGetWpl(pix38);
        l_int32 sum = 0, cnt = 0;
        l_float32 flt00 = 0;
        for (int j = 0; j < h; j++)
        {
            pData14 += j * n15;
            for (int k = 0; k < w; k++)
            {
                l_uint32 r03 = GET_DATA_BYTE(pData14, j) ^ 3;
                if (r03)
                {
                    sum += r03;
                    cnt++;
                    if (r03 > flt00)
                        flt00 = r03;
                }
            }
        }

        l_float32 flt20;
        if (cnt)
        {
            l_float32 flt02 = sum / cnt;
            flt20 = flt02 / 255.f;
        }
        else
            flt20 = 0;
        
        flt00 = flt00 / 255.f;
        numaAddNumber(*ppNuma_ArgC, flt00);
        numaAddNumber(*ppNuma_Arg10, flt20);
        Pix *pix37 = pixClipMaskedWithoutInvert(pix_Arg8, pix43, x, y);

        l_int32 n36, n35;
        pixCountPixels(pix43, &n36, NULL);
        pixCountPixels(pix37, &n35, NULL);
        numaAddNumber(*ppNuma_Arg14, n35 / n36);
        boxDestroy(&pBox44);
        pixDestroy(&pix43);
        pixDestroy(&pix37);
        pixDestroy(&pix38);
    }

    printTimer(timer30, "\tcalculate components score");
    boxaDestroy(&pBoxa45);
}

PixBinarizer::~PixBinarizer() {
}

