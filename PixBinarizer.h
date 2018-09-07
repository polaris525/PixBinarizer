/*
 * PixBinarizer.h
 *
 *  Created on: Jul 8, 2014
 *      Author: renard
 */

#ifndef PIXBINARIZER_H_
#define PIXBINARIZER_H_

#include "allheaders.h"
#include "PixAdaptiveBinarizer.h"

class PixBinarizer {
public:
	PixBinarizer(bool debug);
    Pix* binarize(Pix* pix,  void(*previewCallBack) (Pix*));

	virtual ~PixBinarizer();

private:
	void binarizeInternal(Pix* pixGrey, Pix* pixhm, Pix** pixb);
	Pix* binarizeTiled(Pix* pixs, const l_uint32 tileSize);
	Pix* createEdgeMask(Pix* pixs);
	Pix* combineBinaryImages()
	int determineThresholdForTile(Pix* pixt, bool debug);

	Pix* findGrayScaleEdges(Pix* pixs);
	Pix* smoothGreyPix(Pix*);
	void calculateBinarisationThresholdParams(Pix* pixs);
	Pix* binariseGreyscaleEdges(Pix*);
	Pix* createEdgeSizeMask(Pix* pixs);
	BOXA* findConnCompSizeByRank(Pix*, double);
	Pix* makeSizeMaskFromBoxa(BOXA*, int, int, int);
	void renderBoxaForDebugging(Pix*, BOXA*, const char*);
	Pix* combineBinaryImages(Pix*, Pix*, Pix*, Pix*);
	Pix* filterTextComponentsCalculateScore(Pix*, Pix*, Pix*, Numaa**, string);
	Pix* combineTextComponents(Pixa*, Pixa*, int, int, Numaa*, Numaa*);
	Pixa* removeNonTextComponentsFromPixa(int, int, Pixa*);
	void calculateComponentScore(Pixa*, Pix*, Pix*, Numa**, Numa**, Numa**);

	l_float32 m_flt24;
	l_float32 m_flt20;
	l_float32 m_flt1C;
	l_float32 m_flt18;
	l_int32 m_edgeH_10;
	Pix* m_pix0C;
	PixAdaptiveBinarizer m_pixAdaptiveBin08;
	bool mDebug;
};

#endif /* PIXBINARIZER_H_ */
