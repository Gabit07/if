#include "../pch.h"
#include "CInspectFunc.h"
#include "Fchain.h"
#include "../common/CommonUtil.h"
#include "../common/configure.h"
#include "../common/CommonDefine.h"
#include <math.h>
#include <omp.h>
#include <algorithm>
#include "../CMyImageBuffer.h"
#include <new>
#include <atomic>

// furie transform
PST_FURIEDATA g_ptagRowData;
PST_FURIEDATA g_ptagColData;

double *m_dCol_re;
double *m_dCol_im;

double *m_dRow_re;
double *m_dRow_im;

std::mutex mutex;

void  _fastcall CInspectFunc::CopyImage(LPBYTE psrc, LPBYTE pdst, CRoiClass* roi, int _mode)
{
	int width = roi->m_tagRect.rect.Width();
	int height = roi->m_tagRect.rect.Height();

	int left = roi->m_tagRect.rect.left;
	int top = roi->m_tagRect.rect.top;
	int right = roi->m_tagRect.rect.right;
	int bottom = roi->m_tagRect.rect.bottom;

	if (roi->m_bUsePolyRoi)
	{
		for (int y = top; y < bottom; y++)
		{
			for (int x = left; x < right; x++)
			{
				if (CInspectFunc::CnCheckArea(roi->m_polyRoi, roi->m_nPolyCnt, x, y))
					*(pdst + width * (y - top) + (x - left)) = *(psrc + IMAGE_BUFFER_SIZE_X * y + x); 
				else
					*(pdst + width * (y - top) + (x - left)) = 0;

				
				
			}

		}
	}
	else
	{
		if (roi->m_bUseDonCare)
		{
			for (int y = top; y < bottom; y++)
			{
				for (int x = left; x < right; x++)
				{
					for (int k = 0; k < roi->m_nDonCareAreaCnt; k++)
					{
						if (CInspectFunc::CnCheckArea(roi->m_pointDoncare[k], roi->m_nDoncarePtCnt[k], x, y))
						{
							if (_mode == FIND_LOW)
							{
								*(pdst + width * (y - top) + (x - left)) = 255;
							}
							else if(_mode == FIND_HIGH)
							{
								*(pdst + width * (y - top) + (x - left)) = 0;
							}
							
							break;
						}
						else
						{
							*(pdst + width * (y - top) + (x - left)) = *(psrc + IMAGE_BUFFER_SIZE_X * y + x);
						}
							
							
					}
					
				}
			}
		}
		else
		{
			

			for (int y = top; y < bottom; y++)
			{
				for (int x = left; x < right; x++)
				{
					*(pdst + width * (y - top) + (x - left)) = *(psrc + IMAGE_BUFFER_SIZE_X * y + x);
				}
			}
		}
		
	}


}

void  _fastcall CInspectFunc::CopyImage(LPBYTE psrc, LPBYTE pdst, CRect& rect, int _mode)
{
	int width = rect.Width();
	int height = rect.Height();

	int left = rect.left;
	int top = rect.top;
	int right = rect.right;
	int bottom = rect.bottom;



	for (int y = top; y < bottom; y++)
	{
		for (int x = left; x < right; x++)
		{
			*(pdst + width * (y - top) + (x - left)) = *(psrc + IMAGE_BUFFER_SIZE_X * y + x);
		}
	}

}

//void _fastcall CInspectFunc::Binary(LPBYTE pSrc, LPBYTE pDst, CRoiClass* roi, int nTh, int nMode)
//{
//	int width = roi->m_tagRect.rect.Width();
//	int height = roi->m_tagRect.rect.Height();
//
//
//	/*alignas(std::hardware_destructive_interference_size) */LPBYTE pfm = new BYTE[width * height];
//	memset(pfm, 0, sizeof(BYTE) * width * height);
//
//	{
//		std::scoped_lock<std::mutex> m(mutex);
//		CopyImage(pSrc, pfm, roi, nMode);
//	}
//
//	CMyImageBuffer Imagebuffer;
//
//	Imagebuffer.setPosition(roi->m_tagRect.rect.left, roi->m_tagRect.rect.right, roi->m_tagRect.rect.top, roi->m_tagRect.rect.bottom);
//
//	Imagebuffer.Input2DvectorImage(pfm);
//
//	Imagebuffer.binary(nTh, nMode);
//	Imagebuffer.erode();
//
//	{
//		std::scoped_lock<std::mutex> m(mutex);
//		Imagebuffer.output2DVectorImage(pDst);
//	}
//
//	Imagebuffer.deleteBuffer();
//
//	delete[] pfm;
//
//}

void _fastcall CInspectFunc::Binary(LPBYTE pSrc, LPBYTE pDst, CRect& rect, int nTh, int nMode)
{
	int width = rect.Width();//roi->m_tagRect.rect.Width();
	int height = rect.Height();//roi->m_tagRect.rect.Height();

	
	/*alignas(std::hardware_destructive_interference_size) */LPBYTE pfm = new BYTE[width * height];
	memset(pfm, 0, sizeof(BYTE) * width * height);

	{
		std::scoped_lock<std::mutex> m(mutex);
		CopyImage(pSrc, pfm, rect, nMode);
	}
	

	CMyImageBuffer Imagebuffer;
	
	Imagebuffer.setPosition(rect.left, rect.right, rect.top, rect.bottom);
	//Imagebuffer.setPosition(roi->m_tagRect.rect.left, roi->m_tagRect.rect.right, roi->m_tagRect.rect.top, roi->m_tagRect.rect.bottom);
	
	Imagebuffer.Input2DvectorImage(pfm);
	
	Imagebuffer.binary(nTh, nMode);
	Imagebuffer.erode();

	{
		std::scoped_lock<std::mutex> m(mutex);
		Imagebuffer.output2DVectorImage(pDst);
	}
	
	Imagebuffer.deleteBuffer();

	delete[] pfm;

}

void _fastcall CInspectFunc::Binary_th(LPBYTE pSrc, LPBYTE pDst, CRoiClass* roi, int nTh, int nMode)
{
	int width = roi->m_tagRect.rect.Width();;
	int height = roi->m_tagRect.rect.Height();

	LPBYTE pfm = new BYTE[width * height];
	memset(pfm, 0, sizeof(BYTE) * width * height);

	CopyImage(pSrc, pfm, roi, nMode);

	CMyImageBuffer Imagebuffer;

	Imagebuffer.setPosition(roi->m_tagRect.rect.left, roi->m_tagRect.rect.right, roi->m_tagRect.rect.top, roi->m_tagRect.rect.bottom);

	Imagebuffer.Input2DvectorImage(pfm);

	Imagebuffer.binary(nTh, nMode);
	Imagebuffer.erode();

	Imagebuffer.output2DVectorImage(pDst);

	Imagebuffer.deleteBuffer();

	delete[] pfm;

}

const int GausianMeanMask[3][3] = {{1,1,1},{1,1,1},{1,1,1}}; // x1/8
const int GausianWeightMask[3][3] = {{1,2,1},{2,4,2},{1,2,1}}; // x1/16
void _fastcall CInspectFunc::GausianMeanFilter(LPBYTE pSrc, LPBYTE pDst, CRect rcArea)
{
	int left,right,top,bottom;
	int pitch;

	left = rcArea.left;
	right = rcArea.right;
	top = rcArea.top;
	bottom = rcArea.bottom;

	int nwidth = right-left;
	int nheight = bottom-top;

	pitch = nwidth;

	for(int y=top+1; y<bottom-1; y++)
	{
		for(int x=left+1; x<right-1; x++)
		{
			int ngray=0;
			for(int masky=-1; masky <=1; masky++)
			{
				for(int maskx = -1; maskx <= 1; maskx++)
				{
					ngray += *(pSrc+pitch*(y+masky)+(x+maskx)) * GausianMeanMask[maskx+1][masky+1];
				}
			}
			ngray /= 8;
			if(ngray > 255) ngray = 255;
			*(pDst+pitch*y+x) = ngray;
		}
	}
}
/************************************************************************/
/* 중간값 필터는 노이즈에 강하나 느림, 속도 개선위해 좀 바꿔볼것    */
/************************************************************************/
void _fastcall CInspectFunc::MidFilter(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int _objectNum)
{

	CMyImageBuffer Imagebuffer;

	Imagebuffer.setPosition(rcArea.left, rcArea.right, rcArea.top, rcArea.bottom);

	Imagebuffer.Input2DvectorImage(pSrc);

	Imagebuffer.medianfilter();

	Imagebuffer.output2DVectorImage(pSrc);

	Imagebuffer.deleteBuffer();

	/*int left,right,top,bottom;
	int pitch;

	left = rcArea.left;
	right = rcArea.right;
	top = rcArea.top;
	bottom = rcArea.bottom;

	int nwidth = right-left;
	int nheight = bottom-top;

	pitch = nwidth;

	std::vector<BYTE> v;
	int ngray=0;
	int compare[9];
	ZeroMemory(compare, _countof(compare));
	int cnt = 0;
	register int y, x, maskx, masky;

	for(y=top+1; y<bottom-1; y++)
	{
		for(x=left+1; x<right-1; x++)
		{
			cnt = 0;
			for(masky=0; masky < 3; masky++)
			{
				for(maskx = 0; maskx < 3; maskx++)
				{
					ngray = *(pSrc+pitch*(y+masky)+(x+maskx));
					compare[cnt] = ngray;
					cnt++;
				}
				
			}
			std::qsort(compare, sizeof(compare)/sizeof(int), sizeof(int), CommonUtil::Compare);
			
			*(pDst+pitch*y+x) = compare[_objectNum];
		}
	}*/

}

const int SobelMask_X[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
const int SobelMask_Y[3][3] = {{-1,-2,-1},{0,0,0},{1,2,1}};
void _fastcall CInspectFunc::SobelSharpness(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int direction)
{

	CMyImageBuffer Imagebuffer;

	Imagebuffer.setPosition(rcArea.left, rcArea.right, rcArea.top, rcArea.bottom);

	Imagebuffer.Input2DvectorImage(pSrc);

	Imagebuffer.sobelfilter(DirX);

	Imagebuffer.output2DVectorImage(pSrc);

	Imagebuffer.deleteBuffer();

	/*int left,right,top,bottom;
	int pitch;

	left = rcArea.left;
	right = rcArea.right;
	top = rcArea.top;
	bottom = rcArea.bottom;

	int nwidth = right-left;
	int nheight = bottom-top;

	pitch = nwidth;

	int* pImgSobel = new int[nwidth*nheight];
	memset(pImgSobel, 0, sizeof(int)*nwidth*nheight);

	int* pImgSobelX = new int[nwidth*nheight];
	memset(pImgSobelX, 0, sizeof(int)*nwidth*nheight);

	if(direction == DirY || direction == DirXY)
	{
		for(int y=top+1; y<bottom-1; y++)
		{
			for(int x=left+1; x<right-1; x++)
			{
				int nDx=0;

				for(int masky=-1; masky <=1; masky++)
				{
					for(int maskx = -1; maskx <= 1; maskx++)
					{
						nDx += *(pSrc+pitch*(y+masky)+(x+maskx)) * SobelMask_X[maskx+1][masky+1];
					}
				}
				
				*(pImgSobelX+nwidth*(y-top)+(x-left)) = nDx;

			}
		}
	}
	
	int* pImgSobelY = new int[nwidth*nheight];
	memset(pImgSobelY, 0, sizeof(int)*nwidth*nheight);

	if(direction == DirX || direction == DirXY)
	{
		for(int y=top+1; y<bottom-1; y++)
		{
			for(int x=left+1; x<right-1; x++)
			{
				int nDy=0;

				for(int masky=-1; masky <=1; masky++)
				{
					for(int maskx = -1; maskx <= 1; maskx++)
					{
						nDy += *(pSrc+pitch*(y+masky)+(x+maskx)) * SobelMask_Y[maskx+1][masky+1];
					}
				}
				
				*(pImgSobelY+nwidth*(y-top)+(x-left)) = nDy;

			}
		}
	}
	

	int nMin = 255;
	int nMax = 0;

	for(int i=0; i<nheight; i++)
	{
		for(int j=0; j<nwidth; j++)
		{
			if(direction == DirX)
			    *(pImgSobel+nwidth*i+j) =  *(pImgSobelY+nwidth*i+j);
			else if(direction == DirY)
				*(pImgSobel+nwidth*i+j) =  *(pImgSobelX+nwidth*i+j); 
			else
				*(pImgSobel+nwidth*i+j) =  *(pImgSobelY+nwidth*i+j) + *(pImgSobelX+nwidth*i+j); 

			if(nMax < *(pImgSobel+nwidth*i+j))
				nMax = *(pImgSobel+nwidth*i+j);
			if(nMin > *(pImgSobel+nwidth*i+j))
				nMin = *(pImgSobel+nwidth*i+j);
		}
	}
	
	double dVal1 = double(255.0/double(nMax));
	double dVal2 = double(-255.0/double(nMin));


	int newVal = 0;

	for(int i=0; i<nheight; i++)
	{
		for(int j=0; j<nwidth; j++)
		{
			newVal = *(pImgSobel+nwidth*i+j);
			newVal = int(newVal*dVal1 + dVal2);
			
			if(newVal > 255) newVal = 255;
			if(newVal < 0) newVal = 0;
			*(pDst+pitch*(rcArea.top+i)+(rcArea.left+j)) = newVal;
		}
	}

	delete [] pImgSobel;
	delete [] pImgSobelX;
	delete [] pImgSobelY;*/

}

/************************************************************************/
/* negative 중앙값은 원 영상에 - 함                                                                     */
const int NegativeLaplaMask[3][3] = {{0,1,0},{1,-4,1},{0,1,0}};  // x,y 방향 엣지
const int NegativeLaplaMask2[3][3] = {{1,1,1},{1,-8,1},{1,1,1}}; // x,y,대각선 방향 엣지

/************************************************************************/
/************************************************************************/
/* 중앙값은 원 영상에 + 함                                                                     */
const int LaplaMask[3][3] = {{0,-1,0},{-1,4,-1},{0,-1,0}};  // x,y 방향 엣지
const int LaplaMask2[3][3] = {{-1,-1,-1},{-1,8,-1},{-1,-1,-1}}; // x,y,대각선 방향 엣지
/************************************************************************/
// 라플라시안은 가우시언 필터를 거쳐야 정상적인 값이 나옴
void _fastcall CInspectFunc::LaplacianEdge(LPBYTE pSrc, int* pDst, CRect rcArea)
{

	CMyImageBuffer Imagebuffer;

	Imagebuffer.setPosition(rcArea.left, rcArea.right, rcArea.top, rcArea.bottom);

	Imagebuffer.Input2DvectorImage(pSrc);

	Imagebuffer.laplacian();

	Imagebuffer.output2DVectorImage(pSrc);

	Imagebuffer.deleteBuffer();

}

int l_SumVer[IMAGE_BUFFER_SIZE_X];
int _fastcall CInspectFunc::Erode(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int nSizeX, int nSizeY)
{
	CMyImageBuffer ImageBuffer;
	ImageBuffer.setPosition(rcArea.left, rcArea.right, rcArea.top, rcArea.bottom);
	ImageBuffer.Input2DvectorImage(pSrc);
	ImageBuffer.erode();
	ImageBuffer.output2DVectorImage(pSrc);
	ImageBuffer.deleteBuffer();

	return 0;
}

int _fastcall CInspectFunc::Dilate(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int nSizeX, int nSizeY)
{

	CMyImageBuffer ImageBuffer;
	
	ImageBuffer.setPosition(rcArea.left, rcArea.right, rcArea.top, rcArea.bottom);
	
	ImageBuffer.Input2DvectorImage(pSrc);
	
	ImageBuffer.dilate();
	
	ImageBuffer.output2DVectorImage(pSrc);
	
	ImageBuffer.deleteBuffer();

	return 0;
}


void _fastcall CInspectFunc::FindBlob(LPBYTE pSrc, CRoiClass* roi, int _th, int _mode, std::vector<CRect>& v)
{
	int width = roi->m_tagRect.rect.Width();;
	int height = roi->m_tagRect.rect.Height();

	/*alignas (64) */LPBYTE pfm = new BYTE[width*height];
	memset(pfm, 0, sizeof(BYTE) * width * height);

	CopyImage(pSrc, pfm, roi, _mode);

	CMyImageBuffer ImageBuffer;

	ImageBuffer.setPosition(roi->m_tagRect.rect.left, roi->m_tagRect.rect.right, roi->m_tagRect.rect.top, roi->m_tagRect.rect.bottom);

	ImageBuffer.Input2DvectorImage(pfm);

	ImageBuffer.binary(_th, _mode);

	ImageBuffer.output2DvectorInstanceImage(pfm);

	ImageBuffer.deleteBuffer();

	int left=0,right=0,top=0,bottom=0;
	int pitch=0;

	int nSttPosX = roi->m_tagRect.rect.left;
	int nSttPosY = roi->m_tagRect.rect.top;

	int nKernel = 1;
	int nwidth = roi->m_tagRect.rect.right- roi->m_tagRect.rect.left;
	int nheight = roi->m_tagRect.rect.bottom- roi->m_tagRect.rect.top;

	CChain chain(2000, 200000);
	chain.SetChainData(1, pfm, 1, 1, 4, 20000, nwidth, nheight);
	
	int nBlobCnt = chain.FastChain(0,0,nwidth, nheight); // erode 영향으로 테두리가 희게 떠서
	int left_,right_,top_,bottom_;
	
	for(int i=0; i<nBlobCnt; i++)
	{
		left_ = chain.FindMinX(i);
		right_ = chain.FindMaxX(i);
		top_ =   chain.FindMinY(i);
		bottom_ = chain.FindMaxY(i);
		
		CRect rcLabel;
		
		rcLabel.left =   left_ + nKernel;
		rcLabel.right =  right_ - nKernel;
		rcLabel.top =    top_ + nKernel;
		rcLabel.bottom = bottom_ - nKernel;
		rcLabel.OffsetRect(CPoint(nSttPosX, nSttPosY));
		v.push_back(rcLabel);
		
	}

	delete[] pfm;
}

int _fastcall CInspectFunc::LineSlopeToEdgeX_LeftToRight(LPBYTE pSrc, CRect rcArea, int nposy, double* _data, int* nEdgePos, double dSlope)
{
	int left=0,right=0,top=0,bottom=0;
	int pitch;


	int nwidth = rcArea.right-rcArea.left;
	int nheight = rcArea.bottom-rcArea.top;
	pitch = nwidth;


	int nEdgeCnt = 0;

	for(int x=rcArea.left; x<rcArea.right; x++)
	{
		_data[nEdgeCnt] = *(pSrc+pitch*(nposy-rcArea.top)+(x-rcArea.left));	
		nEdgeCnt++;
	}

	int nMaxPos=0;
	double _dslope = 0, dMaxSlope=0;
	for(int n=2; n<nEdgeCnt-1; n++)
	{
		_dslope = _data[n+1]-_data[n-1];
		if(fabs(_dslope) >= dSlope)
		{
			if(nMaxPos < n)
			{
				nMaxPos = n;
				dMaxSlope = _dslope;
				break;

			}
		}
	}

	if(dMaxSlope == 0)
	{
		return 0;
	}


	if(nMaxPos <= 2 || nMaxPos >= nEdgeCnt-2)
	{
		*nEdgePos = nMaxPos;
	}
	else
	{
		if(fabs(_data[nMaxPos+2]-_data[nMaxPos]) < fabs(_data[nMaxPos]-_data[nMaxPos-2]))
		{
			nMaxPos-=1;
		}

		double Temp_A=(double)fabs(_data[nMaxPos+2]+_data[nMaxPos-2]-2*_data[nMaxPos]);
		double Temp_B=(double)fabs(_data[nMaxPos+3]+_data[nMaxPos-1]-2*_data[nMaxPos+1]);

		if((Temp_A+Temp_B)!=0)
			*nEdgePos=(int)(nMaxPos+Temp_A/(Temp_A+Temp_B));	
		else
			*nEdgePos=(int)nMaxPos;
	}
	return 1;
}

int _fastcall CInspectFunc::LineSlopeToEdgeX_RightToLeft(LPBYTE pSrc, CRect rcArea, int nposy, double* _data, int* nEdgePos, double dSlope)
{
	int left=0,right=0,top=0,bottom=0;
	int pitch;


	int nwidth = rcArea.right-rcArea.left;
	int nheight = rcArea.bottom-rcArea.top;
	pitch = nwidth;


	int nEdgeCnt = nwidth;

	for(int x=rcArea.right; x>rcArea.left; x--)
	{
		_data[nEdgeCnt] = *(pSrc+pitch*(nposy-rcArea.top)+(x-rcArea.left));	
		nEdgeCnt--;
	}

	int nMaxPos=0;
	double _dslope = 0, dMaxSlope=0;
	for(int n=nwidth-2; n>1; n--)
	{
		_dslope = _data[n+1]-_data[n-1];
		if(fabs(_dslope) >= dSlope)
		{
			if(nMaxPos < n)
			{
				nMaxPos = n;
				dMaxSlope = _dslope;
				break;

			}
		}
	}

	if(dMaxSlope == 0)
	{
		return 0;
	}


	if(nMaxPos <= 2 || nMaxPos >= nEdgeCnt-2)
	{
		*nEdgePos = nMaxPos;
	}
	else
	{
		if(fabs(_data[nMaxPos+2]-_data[nMaxPos]) < fabs(_data[nMaxPos]-_data[nMaxPos-2]))
		{
			nMaxPos-=1;
		}

		double Temp_A=(double)fabs(_data[nMaxPos+2]+_data[nMaxPos-2]-2*_data[nMaxPos]);
		double Temp_B=(double)fabs(_data[nMaxPos+3]+_data[nMaxPos-1]-2*_data[nMaxPos+1]);

		if((Temp_A+Temp_B)!=0)
			*nEdgePos=(int)(nMaxPos+Temp_A/(Temp_A+Temp_B));	
		else
			*nEdgePos=(int)nMaxPos;
	}
	return 1;
}

int _fastcall CInspectFunc::LineSlopeToEdgeY_TopToBottom(LPBYTE pSrc, CRect rcArea, int nposx, double* _data, int* nEdgePos, double dSlope)
{
	int left=0,right=0,top=0,bottom=0;
	int pitch;


	int nwidth = rcArea.right-rcArea.left;
	int nheight = rcArea.bottom-rcArea.top;
	pitch = nwidth;


	int nEdgeCnt = 0;

	for(int y=rcArea.top; y<rcArea.bottom; y++)
	{
		_data[nEdgeCnt] = *(pSrc+pitch*(y-rcArea.top)+(nposx-rcArea.left));	
		nEdgeCnt++;
	}

	int nMaxPos=0;
	double _dslope = 0, dMaxSlope=0;
	for(int n=2; n<nEdgeCnt-1; n++)
	{
		_dslope = _data[n+1]-_data[n-1];
		if(fabs(_dslope) >= dSlope)
		{
			if(nMaxPos < n)
			{
				nMaxPos = n;
				dMaxSlope = _dslope;
				break;
			}
		}
	}

	if(dMaxSlope == 0)
	{
		return 0;
	}


	if(nMaxPos <= 2 || nMaxPos >= nEdgeCnt-2)
	{
		*nEdgePos = nMaxPos;
	}
	else
	{
		if(fabs(_data[nMaxPos+2]-_data[nMaxPos]) < fabs(_data[nMaxPos]-_data[nMaxPos-2]))
		{
			nMaxPos-=1;
		}

		double Temp_A=(double)fabs(_data[nMaxPos+2]+_data[nMaxPos-2]-2*_data[nMaxPos]);
		double Temp_B=(double)fabs(_data[nMaxPos+3]+_data[nMaxPos-1]-2*_data[nMaxPos+1]);

		if((Temp_A+Temp_B)!=0)
			*nEdgePos=(int)(nMaxPos+Temp_A/(Temp_A+Temp_B));	
		else
			*nEdgePos=(int)nMaxPos;
	}

	return 1;
}

int _fastcall CInspectFunc::LineSlopeToEdgeY_BottomToTop(LPBYTE pSrc, CRect rcArea, int nposx, double* _data, int* nEdgePos, double dSlope)
{
	int left=0,right=0,top=0,bottom=0;
	int pitch;


	int nwidth = rcArea.right-rcArea.left;
	int nheight = rcArea.bottom-rcArea.top;
	pitch = nwidth;


	int nEdgeCnt = rcArea.bottom-rcArea.top;

	for(int y=rcArea.bottom; y>rcArea.top; y--)
	{
		_data[nEdgeCnt-1] = *(pSrc+pitch*(y-1-rcArea.top)+(nposx-rcArea.left));	
		nEdgeCnt--;
	}

	int nMaxPos=0;
	double _dslope = 0, dMaxSlope=0;
	for(int n=rcArea.bottom-rcArea.top-2; n>1; n--)
	{
		_dslope = _data[n+1]-_data[n-1];
		if(fabs(_dslope) >= dSlope)
		{
			if(nMaxPos < n)
			{
				nMaxPos = n;
				dMaxSlope = _dslope;
				break;
			}
		}
	}

	if(dMaxSlope == 0)
	{
		return 0;
	}


	if(nMaxPos <= 2 || nMaxPos >= nEdgeCnt-2)
	{
		*nEdgePos = nMaxPos;
	}
	else
	{
		if(fabs(_data[nMaxPos+2]-_data[nMaxPos]) < fabs(_data[nMaxPos]-_data[nMaxPos-2]))
		{
			nMaxPos-=1;
		}

		double Temp_A=(double)fabs(_data[nMaxPos+2]+_data[nMaxPos-2]-2*_data[nMaxPos]);
		double Temp_B=(double)fabs(_data[nMaxPos+3]+_data[nMaxPos-1]-2*_data[nMaxPos+1]);

		if((Temp_A+Temp_B)!=0)
			*nEdgePos=(int)(nMaxPos+Temp_A/(Temp_A+Temp_B));	
		else
			*nEdgePos=(int)nMaxPos;
	}

	return 1;
}

void _fastcall CInspectFunc::FindEdgeLine(LPBYTE pSrc, CRect rcArea, int N, int _direction, int _direction2, int nSlope, double* Alpha, double* Beta, double* C)
{
	int left=0,right=0,top=0,bottom=0;

	int nwidth = rcArea.right-rcArea.left;
	int nheight = rcArea.bottom-rcArea.top;
	
	double *_data = new double[nwidth*nheight];
	memset(_data,0,sizeof(double)*nwidth*nheight);

	int cnt=0;
	int nEdgePos=0;
	int nRes = -1;

	if(_direction == DirX) // x방향 탐색
	{
		double *_VerX = new double[nwidth*nheight];
		memset(_VerX,0,sizeof(double)*nwidth*nheight);

		double *_VerY = new double[nwidth*nheight];
		memset(_VerY,0,sizeof(double)*nwidth*nheight);

		if(_direction2 == eLEFT_TO_RIGHT)
		{
			for(int y=rcArea.top; y<rcArea.bottom; y++)
			{
				nRes = LineSlopeToEdgeX_LeftToRight(pSrc, rcArea, y, _data, &nEdgePos, 20);
				if(nRes)
				{
					_VerX[cnt] = nEdgePos;
					_VerY[cnt] = y;
					cnt++;
				}
			}
		}
		else if(_direction2 == eRIGHT_TO_LEFT)
		{
			for(int y=rcArea.top; y<rcArea.bottom; y++)
			{
				nRes = LineSlopeToEdgeX_RightToLeft(pSrc, rcArea, y, _data, &nEdgePos, 20);
				if(nRes)
				{
					_VerX[cnt] = nEdgePos;
					_VerY[cnt] = y;
					cnt++;
				}
			}
		}
		else
		{
			AfxMessageBox(L"잘못된 설정 입니다.");
			goto done;
		}
		

		int nHistoRange=(rcArea.right-rcArea.left)/5;
		if(nHistoRange<6) nHistoRange=6;
		else if(nHistoRange>(rcArea.right-rcArea.left)/2) nHistoRange=(rcArea.right-rcArea.left)/2;

		N = ChooseGoodData_FirstBase(cnt, _VerX, _VerY, nHistoRange);

		LineFitting(N, nwidth/2, _direction, _VerY, _VerX, Alpha, Beta, C);

		double t,a,b;
		if(_direction == DirX)
		{
			t=1;
			a=*C/(*Alpha);
			b=-(*Beta)/(*Alpha);

			*C = t;
			*Alpha = a;
			*Beta = b;
		}
done:
		delete[] _VerX;
		delete[] _VerY;
	}
	else
	{
		double *_HorX = new double[nwidth*nheight];
		memset(_HorX,0,sizeof(double)*nwidth*nheight);

		double *_HorY = new double[nwidth*nheight];
		memset(_HorY,0,sizeof(double)*nwidth*nheight);

		if(_direction2 == eTOP_TO_BOTTOM)
		{
			for(int x=rcArea.left; x<rcArea.right; x++)
			{
				nRes = LineSlopeToEdgeY_TopToBottom(pSrc, rcArea, x, _data, &nEdgePos, 25);
				if(nRes)
				{
					_HorX[cnt] = x;
					_HorY[cnt] = nEdgePos;
					cnt++;
				}
			}
		}
		else if(_direction2 == eBOTTOM_TO_TOP)
		{
			for(int x=rcArea.left; x<rcArea.right; x++)
			{
				nRes = LineSlopeToEdgeY_BottomToTop(pSrc, rcArea, x, _data, &nEdgePos, 25);
				if(nRes)
				{
					_HorX[cnt] = x;
					_HorY[cnt] = nEdgePos;
					cnt++;
				}
			}
		}
		else
		{
			AfxMessageBox(L"잘못된 설정 입니다.");
			goto done2;
		}
		

		int nHistoRange=(rcArea.bottom-rcArea.top)/5;
		if(nHistoRange<6) nHistoRange=6;
		else if(nHistoRange>(rcArea.bottom-rcArea.top)/2) nHistoRange=(rcArea.bottom-rcArea.top)/2;

		N = ChooseGoodData_FirstBase(cnt, _HorY, _HorX, nHistoRange);

		LineFitting(N, nheight/2, _direction, _HorX, _HorY, Alpha, Beta, C);

done2:
		delete[] _HorX;
		delete[] _HorY;
		
	}

	delete[] _data;
}

BOOL _fastcall CInspectFunc::LineFitting(int N, int removeCnt, int _dir, double* x, double* y, double* Alpha, double* Beta, double* C)
{
	double sigma_xx=0;
	double sigma_xy=0;
	double sigma_x=0;
	double sigma_y=0;
	double sigma_N=0;
	int MaxPos = 0;
	BOOL bSkip[IMAGE_BUFFER_SIZE_X]={FALSE,};

	for(int d=0; d<removeCnt;d++)
	{
		if(d==0)
		{
			for(int n=0; n<N; n++)
			{
				sigma_xx += x[n]*x[n];
				sigma_xy += x[n]*y[n];
				sigma_x += x[n];
				sigma_y += y[n];
				sigma_N += 1;
			}
		}
		else
		{
			sigma_xx -= *(x+MaxPos)**(x+MaxPos);
			sigma_xy -= *(x+MaxPos)**(y+MaxPos);
			sigma_x -= *(x+MaxPos);
			sigma_y -= *(y+MaxPos);
			sigma_N -= 1;
		}
		// 선형회귀 추정 방정식
		double Det = sigma_N * sigma_xx - sigma_x*sigma_x; 

		if(Det == 0) return 0; // 역행렬 존재하지 않을때

		if(fabs(Det)<0.00001)
		{
			*C=0;
			*Alpha=1;
			*Beta=-sigma_x/sigma_N;
			goto done;
		}
		else
		{
			*C=1;
			*Alpha = (sigma_N*sigma_xy -sigma_x*sigma_y) / Det;
			*Beta = ((sigma_y*sigma_xx) + (-sigma_xy*sigma_x)) / Det; 
		}

		double dd =sqrt(*Alpha**Alpha+1);
		int MaxDist=0;
		for(int j=0;j<N;j++)
		{
			if(bSkip[j]) continue;
			double dist= fabs(LineToPointDist(*Alpha, *Beta, *C, (int)x[j], (int)y[j], _dir));
			if(dist>MaxDist)
			{
				MaxDist=(int)dist;
				MaxPos=j;
			}
			if(MaxDist/dd<2) 
				goto done;
			else 
				bSkip[MaxPos] = TRUE;
			
		}
		
	}
done:

	return 1;
}

double _fastcall CInspectFunc::LineToPointDist(double Alpha, double Beta, double C, int x, int y, int direction)
{
	double dDistance;
	
	if(direction == DirX)
	{
		dDistance = x-(C*y-Beta)/Alpha;		 //Edge값 - 직선위의 Y값에 해당하는 X : +면 파인값, -면 튀어나온 값  
	}
	else
	{
		dDistance = y-(Alpha*x+Beta)/C;
	}

	return dDistance;
}

int l_histo[IMAGE_BUFFER_SIZE_X/3+10];
int _fastcall CInspectFunc::ChooseGoodData_FirstBase(int N, double *x, double *y, int space)
{
	int i;
	int nMax=0, iMax;		//nMax : Max histo Data, iMax: Max histo 발생지점
	int nTemp;				
	int nCount=0;
	int HISTO_SPACE=space/3;
	int MAX_HISTO;

	if(HISTO_SPACE<3) HISTO_SPACE=3;
	MAX_HISTO=16384/HISTO_SPACE+HISTO_SPACE;

	if(N<=0) return N;

	for(i=0;i<MAX_HISTO;i++)
		l_histo[i]=0;

	for(i=0;i<N;i++)
	{
		nTemp=(int)(x[i]/HISTO_SPACE);
		if(nTemp<MAX_HISTO)
			l_histo[nTemp]++;
	}

	for(i=0;i<MAX_HISTO;i++)
		if(l_histo[i]>nMax)
		{
			nMax=l_histo[i];
			iMax=i;
		}

		for(i=0;i<N;i++)
		{
			if((int)(x[i]/HISTO_SPACE)>=iMax-1 && (int)(x[i]/HISTO_SPACE)<=iMax+1)
			{
				x[nCount]=x[i];
				y[nCount]=y[i];
				nCount++;
			}
		}
		if(nCount>0)
			return nCount;
		else
		{
			HISTO_SPACE=10;
			MAX_HISTO=IMAGE_BUFFER_SIZE_X/HISTO_SPACE+HISTO_SPACE;

			if(N<=0) return N;

			for(i=0;i<MAX_HISTO;i++)
				l_histo[i]=0;

			for(i=0;i<N;i++){
				nTemp=(int)(x[i]/HISTO_SPACE);
				if(nTemp<MAX_HISTO)
					l_histo[nTemp]++;
			}

			for(i=0;i<MAX_HISTO;i++)
				if(l_histo[i]>nMax)
				{
					nMax=l_histo[i];
					iMax=i;
				}

				for(i=0;i<N;i++)
				{
					if((int)(x[i]/HISTO_SPACE)>=iMax-1 && (int)(x[i]/HISTO_SPACE)<=iMax+1)
					{
						x[nCount]=x[i];
						y[nCount]=y[i];
						nCount++;
					}
				}	
				return nCount;
		}
		return nCount;
}

int _fastcall CInspectFunc::FindCrossPoint(double t1, double a1, double b1, double t2, double a2, double b2, double *cx, double *cy)
{
	double local_LIMIT=1e-10;
	double dd;

	dd=-a1*t2+a2*t1;

	if(fabs(dd)<local_LIMIT) return -1;

	*cx=(t2*b1-t1*b2)/dd;
	*cy=(a2*b1-a1*b2)/dd;

	return 0;
}

void _fastcall CInspectFunc::HoughTransLineFitting(LPBYTE pSrc, CRoiClass* roi, int _direction, int _direction2, double _slope, double* A, double* B)
{
	int width = roi->m_tagRect.rect.Width();;
	int height = roi->m_tagRect.rect.Height();

	LPBYTE pfm = new BYTE[width * height];
	memset(pfm, 0, sizeof(BYTE) * width * height);

	CopyImage(pSrc, pfm, roi);

	CMyImageBuffer ImageBuffer;
	ImageBuffer.createhoughLookuptable();

	ImageBuffer.setPosition(roi->m_tagRect.rect.left, roi->m_tagRect.rect.right, roi->m_tagRect.rect.top, roi->m_tagRect.rect.bottom);

	ImageBuffer.Input2DvectorImage(pfm);

	ImageBuffer.houghlinefit(_direction,_direction2, _slope, *A, *B);

	//ImageBuffer.output2DVectorImage(pSrc);

	ImageBuffer.deleteBuffer();
	ImageBuffer.deleteLookuptable();

	delete[] pfm;

//	int left,right,top,bottom;
//	int pitch=0;
//
//	left = rcArea.left;
//	right = rcArea.right;
//	top = rcArea.top;
//	bottom = rcArea.bottom;
//
//	float fRadian = float(PI/180);
//
//	int nTh = 20;
//
//	int nwidth = right-left;
//	int nheight = bottom - top;
//	
//	double dRho = sqrt(double(nwidth*nwidth + nheight*nheight));
//	const int nRho = (int)dRho;
//
//	int** nHoughWeight = new int*[nRho+1];
//	for(int i=0; i<nRho+1; i++)
//	{
//		nHoughWeight[i] = new int[_THETA];
//		memset(nHoughWeight[i], 0, sizeof(int)*_THETA);
//	}
//
//	
//	double *dLUT_SIN = new double[_THETA];
//	double *dLUT_COS = new double[_THETA];
//
//	memset(dLUT_SIN, 0, sizeof(double)*_THETA);
//	memset(dLUT_COS, 0, sizeof(double)*_THETA);
//	
//	for(int i=0; i<_THETA; i++)
//	{
//		dLUT_SIN[i] = sin(i*fRadian);
//		dLUT_COS[i] = cos(i*fRadian);
//	}
//
//	double *_data = new double[nwidth*nheight];
//	memset(_data,0,sizeof(double)*nwidth*nheight);
//
//	int nEdgePosX=0, nEdgePosY=0, nRes=0;
//	double d=0;
//	
//	
//	if(_direction == DirX) // X방향 엣지 탐색
//	{
//		if(_direction2 == LTR)
//		{
//			for(int y=top; y<bottom; y++)
//			{
//				nRes = LineSlopeToEdgeX_LeftToRight(pSrc, rcArea, y, _data, &nEdgePosX, 20);
//				if(nRes)
//				{
//					for(int t=0; t<_THETA; t++)
//					{
//						d = (int)(nEdgePosX*dLUT_COS[t] + y * dLUT_SIN[t]);
//						if(d >= 1 && d <= nRho)
//						{
//							nHoughWeight[(int)d][t]++;
//						}
//					}	
//				}
//			}
//		}
//		else if(_direction2 == RTL)
//		{
//			for(int y=top; y<bottom; y++)
//			{
//				nRes = LineSlopeToEdgeX_RightToLeft(pSrc, rcArea, y, _data, &nEdgePosX, 20);
//				if(nRes)
//				{
//					for(int t=0; t<_THETA; t++)
//					{
//						d = (int)(nEdgePosX*dLUT_COS[t] + y * dLUT_SIN[t]);
//						if(d >= 1 && d <= nRho)
//						{
//							nHoughWeight[(int)d][t]++;
//						}
//					}	
//				}
//			}
//		}
//		else
//		{
//			AfxMessageBox(L"잘못된 설정 입니다.");
//			goto done;
//		}
//		
//	}
//	else if(_direction == DirY) // Y 방향 엣지 탐색
//	{
//		if(_direction2 == TTB)
//		{
//			for(int x=left; x<right; x++)
//			{
//				nRes = LineSlopeToEdgeY_TopToBottom(pSrc, rcArea, x, _data, &nEdgePosY, 15);
//				if(nRes)
//				{
//					for(int t=0; t<_THETA; t++)
//					{
//						d = (int)((x)*dLUT_COS[t] + (nEdgePosY) * dLUT_SIN[t]);
//						if(d >= 1 && d <= nRho)
//						{
//							nHoughWeight[(int)d][t]++;
//						}
//					}	
//				}
//			}
//		}
//		else if(_direction2 == BTT)
//		{
//			for(int x=left; x<right; x++)
//			{
//				nRes = LineSlopeToEdgeY_BottomToTop(pSrc, rcArea, x, _data, &nEdgePosY, 15);
//				if(nRes)
//				{
//					for(int t=0; t<_THETA; t++)
//					{
//						d = (int)((x)*dLUT_COS[t] + (nEdgePosY) * dLUT_SIN[t]);
//						if(d >= 1 && d <= nRho)
//						{
//							nHoughWeight[(int)d][t]++;
//						}
//					}	
//				}
//			}
//		}
//		else
//		{
//			AfxMessageBox(L"잘못된 설정 입니다.");
//			goto done;
//		}
//		
//	}
//	//else // 양방향 탐색
//	//{
//	//	for(int y=top; y<bottom; y++)
//	//	{
//	//		nRes = LineSlopeToEdgeX_LeftToRight(pSrc, rcArea, y, _data, &nEdgePosX, 15);
//	//		if(nRes)
//	//		{
//	//			for(int t=0; t<_THETA; t++)
//	//			{
//	//				d = (int)(nEdgePosX*dLUT_COS[t] + y * dLUT_SIN[t]);
//	//				if(d >= 1 && d <= nRho)
//	//				{
//	//					dHoughW[(int)d][t]++;
//	//				}
//	//			}	
//	//		}
//	//	}
//	//	for(int x=left; x<right; x++)
//	//	{
//	//		nRes = LineSlopeToEdgeY_BottomToTop(pSrc, rcArea, x, _data, &nEdgePosY, 20);
//	//		
//	//		if(nRes)
//	//		{
//	//			for(int t=0; t<_THETA; t++)
//	//			{
//	//				d = (int)((x)*dLUT_COS[t] + (nEdgePosY) * dLUT_SIN[t]);
//	//				if(d >= 1 && d <= nRho)
//	//				{
//	//					dHoughW[(int)d][t]++;
//	//				}
//	//			}	
//	//		}
//	//	}
//	//}
//	
//	int nmax = nTh;
//	int D=0;
//	int nMaxTheta=0;
//
//	for(int t=1; t<_THETA; t++)
//	{
//		for(int i=1; i<nRho; i++)
//		{
//			if(nmax < nHoughWeight[i][t])
//			{
//				nmax = nHoughWeight[i][t];
//				D = i;
//				nMaxTheta = t;
//			}
//		}
//	}
//
//	if(nMaxTheta == 0 || D==0) 
//	{
//		goto done;
//	}
//	
//	*A = -dLUT_COS[nMaxTheta]/dLUT_SIN[nMaxTheta];
//	*B = D/dLUT_SIN[nMaxTheta];
//	
//done:
//
//	delete[] _data;
//	delete[] dLUT_SIN;
//	delete[] dLUT_COS;
//	
//	for(int i=0; i<nRho+1; i++)
//		delete[] nHoughWeight[i];
//	delete[] nHoughWeight;

}


/************************************************************************/
/* Thread 분할하여 회전 하며, 예외처리 구문이 오히려 검정 선을 분할경계에 만들기 때문에
   삭제함*/
/************************************************************************/
void _fastcall CInspectFunc::RotateImage(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int centerx, int centery, double dDegree)
{
	int nwidth = rcArea.right-rcArea.left;
	int nheight = rcArea.bottom-rcArea.top;

	double dSeta = PI / (180 / dDegree);
	double dx, dy;


	//int centerx=0, centery=0;
	//centerx = (rcArea.left+rcArea.right)/2;
	//centery = (rcArea.top + rcArea.bottom)/2;

	for(int i=rcArea.top; i<rcArea.bottom; i++)
	{
		for(int j=rcArea.left; j<rcArea.right; j++)
		{
			dx = (j-centerx) * cos(dSeta) + (i-centery) * sin(dSeta) + centerx;
			dy = -(j-centerx) * sin(dSeta) + (i-centery) * cos(dSeta) + centery;

			//if(dx < 0 || dx >= rcArea.right) continue;
			//if(dy < 0 || dy >= rcArea.bottom) continue;
			
			BiLinearInterpolation(pSrc, pDst, pitch, rcArea, j, i, dx, dy); // 회전영상 보간처리
		}
	}
	
}

void _fastcall CInspectFunc::BiLinearInterpolation(LPBYTE pSrc, LPBYTE pDst, int nPitch, CRect rcArea, int x, int y, double newx, double newy)
{
	int pitch = nPitch;
	
	int x1, x2, y1, y2;
	int nGrayVal;
	double rx, ry, p, q, dVal;

	rx = newx;
	ry = newy;

	x1 = (int)rx;
	y1 = (int)ry;

	x2 = x1+1; //if(x2==nwidth) return;//x2 = nwidth-1;
	y2 = y1+1; //if(y2==nheight) return;//y2 = nheight-1;

	if(x2 <= 0) return;
	if(y2 <= 0) return;

	p = rx - x1;
	q = ry - y1;

	dVal = (1.-p)*(1.-q)**(pSrc+pitch*y1+x1)
		+ p*(1.-q)**(pSrc+pitch*y1+x2)
		+ (1.-p)*q**(pSrc+pitch*y2+x1)
		+ (p)*(q)**(pSrc+pitch*y2+x2);

	nGrayVal = int(dVal+0.5);
	
	if(nGrayVal >= 255) nGrayVal = 255;
	*(pDst+pitch*y+x) = nGrayVal;

}

void _fastcall CInspectFunc::HoughTransFindElipse(LPBYTE pSrc, CRect rcArea, int& centerx, int& centery, int& Radius)
{
	int left = rcArea.left;
	int top = rcArea.top;
	int right = rcArea.right;
	int bottom = rcArea.bottom;

	CString strTemp;

	int nWidth, nHeight, npitch = 0;
	const int nMinR = 20;
	Radius = 0;

	nWidth = right - left;
	npitch = nWidth;
	nHeight = bottom - top;

	int ntemp = (int)sqrt(pow((float)nWidth, 2) + pow((float)nHeight, 2));
	int nMaxDegree = 0;


	int fx = 0;
	int fy = 0;
	int nPosx = 0;
	int nPosy = 0;

	int* H = new int[ntemp * nWidth * nHeight];
	memset(H, 0, sizeof(int)* ntemp * nWidth * nHeight);

	//Concurrency::parallel_for(left, right, [&](int x) // 300~320ms : 100ms 정도의 속도 절약이 있음
	for (int x = left; x < right; x += 2) // 437~ 478ms
	{
		for (int y = top; y < bottom; y += 2)
		{
			if (*(pSrc + npitch * y + x) == 255)
			{
				for (int x0 = x; x0 < right; x0+=2)
				{
					for (int y0 = y; y0 < bottom; y0 += 2)
					{
						int r = (int)sqrt(pow(float(x - x0), 2) + pow(float(y - y0), 2));

						if (abs(r) > nMinR && abs(r) < ntemp / 2)
						{
							int nx = x0 - left;
							int ny = y0 - top;
							H[(r*nWidth + nx)*nHeight + ny]++;


							if (H[(r*nWidth + nx)*nHeight + ny] > 20 && H[(r*nWidth + nx)*nHeight + ny] > nMaxDegree)
							{
								nMaxDegree = H[(r*nWidth + nx)*nHeight + ny];
								Radius = r;
								centerx = nx;
								centery = ny;
							}

						}

					}
				}
			}
		}
	}//);

	delete[] H;

}
// 등록된 포인트 점, 포인트 개수, 탐색 영역 x, y 좌표
BOOL _fastcall CInspectFunc::CnCheckArea(POINT poly[], int nCount, int nX, int nY) 
{
	int i, j;
	BOOL c = FALSE;
	for (i = 0, j = nCount-1; i < nCount; j = i++) {
		if ( (((poly[i].y<=nY) && (nY<poly[j].y)) ||
			((poly[j].y<=nY) && (nY<poly[i].y))) &&
			(nX < (poly[j].x - poly[i].x) * (nY - poly[i].y) /
			(poly[j].y - poly[i].y) + poly[i].x))
			c = !c;
	}
	return c;
}

void _fastcall CInspectFunc::InspectHerrisCorner(LPBYTE pSrc, CRect rcArea, std::vector<CPoint>& v)
{
	int left			= rcArea.left;
	int top				= rcArea.top;
	int right			= rcArea.right;
	int bottom			= rcArea.bottom;

	int npitch = right-left;
	int w = right-left;
	int h = bottom-top;


	float* dx2 = new float[w*h];
	float* dy2 = new float[w*h];
	float* dxy = new float[w*h];

	memset(dx2, 0, sizeof(float)*w*h);
	memset(dy2, 0, sizeof(float)*w*h);
	memset(dxy, 0, sizeof(float)*w*h);

	for(int i=left; i<right; i++)
	{
		for(int j=top; j<bottom; j++)
		{
			int val = *(pSrc + npitch * j + i);

			if(val < 100)
				*(pSrc + npitch * j + i) = 255;
			else
				*(pSrc + npitch * j + i) = 0;


		}
	}

	float tx, ty=0;

	for(int j= top+1; j<bottom-1; j++)
	{
		for(int i=left + 1; i<right-1; i++)
		{
			tx = (float)((*(pSrc+npitch*(j-1)+(i+1)) + *(pSrc+npitch*(j)+(i+1)) + *(pSrc+npitch*(j+1)+(i+1)) - *(pSrc+npitch*(j-1)+(i-1)) - *(pSrc+npitch*(j)+(i-1)) - *(pSrc+npitch*(j+1)+(i+1)))/6.0);
			ty = (float)((*(pSrc+npitch*(j+1)+(i-1)) + *(pSrc+npitch*(j+1)+(i)) + *(pSrc+npitch*(j+1)+(i+1)) - *(pSrc+npitch*(j-1)+(i-1)) - *(pSrc+npitch*(j-1)+(i)) - *(pSrc+npitch*(j-1)+(i+1)))/6.0);

			dx2[w*(j-top)+(i-left)] = tx * tx;
			dy2[w*(j-top)+(i-left)] = ty * ty;
			dxy[w*(j-top)+(i-left)] = tx * ty;
		}
	}

	float g[5][5] = { { 1, 4, 6, 4, 1 },{ 4, 16, 24, 16, 4 },
	{ 6, 24, 36, 24, 6 },{ 4, 16, 24, 16, 4 },{ 1, 4, 6, 4, 1 } };

	for (int y = 0; y < 5; y++)
	{
		for (int x = 0; x < 5; x++)
		{
			g[y][x] /= 256.f;
		}
	}

	float* gdx2 = new float[w*h];
	float* gdy2 = new float[w*h];
	float* gdxy = new float[w*h];

	memset(gdx2, 0, sizeof(float)*w*h);
	memset(gdy2, 0, sizeof(float)*w*h);
	memset(gdxy, 0, sizeof(float)*w*h);


	float tx2, ty2, txy;
	for (int j = top + 2; j < bottom - 2; j++)
	{
		for (int i = left + 2; i < right - 2; i++)
		{
			tx2 = ty2 = txy = 0;
			for (int y = 0; y < 5; y++)
			{
				for (int x = 0; x < 5; x++)
				{
					tx2 += (dx2[w * (j-top + y - 2) + (i-left + x - 2)] * g[y][x]);
					ty2 += (dy2[w * (j-top + y - 2) + (i-left + x - 2)] * g[y][x]);
					txy += (dxy[w * (j-top + y - 2) + (i-left + x - 2)] * g[y][x]);
				}

				gdx2[w*(j-top)+(i-left)] = tx2;
				gdy2[w*(j-top)+(i-left)] = ty2;
				gdxy[w*(j-top)+(i-left)] = txy;
			}
		}
	}


	float* crf = new float[w*h];
	memset(crf, 0, sizeof(float)*w*h);

	float k = 0.04f;
	for (int j = top + 2; j < bottom - 2; j++)
	{
		for (int i = left+2; i < right - 2; i++)
		{
			crf[w * (j-top) + (i-left)] = (gdx2[w * (j-top) + (i-left)] * gdy2[w * (j-top) + (i-left)] - gdxy[w * (j-top) + (i-left)] * gdxy[w * (j-top) + (i-left)])
				- k * (gdx2[w * (j-top) + (i-left)] + gdy2[w * (j-top) + (i-left)]) * (gdx2[w * (j-top) + (i-left)] + gdy2[w * (j-top) + (i-left)]);
		}
	}

	float cvf_value;
	float maxval = 0;

	v.clear();

	for (int j = top + 2; j < bottom - 2; j++)
	{
		for (int i = left+2; i < right - 2; i++)
		{
			cvf_value = crf[w * (j-top) + (i-left)];
			if(cvf_value > maxval)
				maxval = cvf_value;
			if (cvf_value > 20000)
			{
				if (cvf_value > crf[w * (j-top - 1) + (i-left)] && cvf_value > crf[w * (j-top - 1) + (i-left + 1)] &&
					cvf_value > crf[w * (j-top) + (i-left + 1)] && cvf_value > crf[w * (j-top + 1) + (i-left + 1)] &&
					cvf_value > crf[w * (j-top + 1) + (i-left)] && cvf_value > crf[w * (j-top + 1) + (i-left - 1)] &&
					cvf_value > crf[w * (j-top) + (i-left - 1)] && cvf_value > crf[w * (j-top - 1) + (i-left - 1)])
				{
					CPoint point(i,j);
					v.push_back(point);
					
				}
			}
		}
	}

	delete[] dx2;
	delete[] dy2;
	delete[] dxy;
	delete[] gdx2;
	delete[] gdy2;
	delete[] gdxy;
	delete[] crf;

	return;
}

void _fastcall CInspectFunc::FindRoundLine(LPBYTE pSrc, CRect rcArea, int _dir, double &A0, double &A1, double &A2)
{
	int nwidth = rcArea.right-rcArea.left;
	int nheight = rcArea.bottom-rcArea.top;

	int pitch = nwidth;
	
	int centerx=0, centery=0;
	
	int nRes = 0, nCnt=0;
	int nEdgePos;
	centerx = (rcArea.left+rcArea.right)/2;
	centery = (rcArea.top + rcArea.bottom)/2;

	double *_data = new double[nwidth*nheight];
	memset(_data,0,sizeof(double)*nwidth*nheight);
 
	int nC=0; int nMaxX=0;int nY=0;
	std::vector<double> v_x;
	std::vector<double> v_y;
	if(_dir == DirX)
	{
		for(int y=rcArea.top; y<rcArea.bottom; y++)
		{
			nRes = LineSlopeToEdgeX_LeftToRight(pSrc, rcArea, y, _data, &nEdgePos, 20);
			if(nRes)
			{
				v_x.push_back(nEdgePos);
				v_y.push_back(y);
			}
		}
		get2ndOrderRegression(&v_y, &v_x, &A0, &A1, &A2);
	}
	else
	{
		for(int x=rcArea.left; x<rcArea.right; x++)
		{
			nRes = LineSlopeToEdgeY_TopToBottom(pSrc, rcArea, x, _data, &nEdgePos, 20);
			if(nRes)
			{
				v_x.push_back(x);
				v_y.push_back(nEdgePos);
			}
		}
		get2ndOrderRegression(&v_x, &v_y, &A0, &A1, &A2);
	}
	
	//CalcThreePoint(CPoint(_VerX[0],_VerY[0]), CPoint(nMaxX,nY), CPoint(0, nC), &A, &B, &C);
	
	delete[] _data;
}

void _fastcall CInspectFunc::CalcThreePoint(CPoint p1, CPoint p2, CPoint p3, double *Alpha, double *Beta, double *C)
{
	int a, b, c, d;

	a = p1.x*p1.x; b=p2.x*p2.x;
	c = p1.x;      d= p2.x;

	double det = a*d-b*c;
	*Alpha = (d*(p1.y-p3.y) - b*(p2.y-p3.y))/det;
	*Beta = (-c*(p1.y-p3.y) + a*(p2.y-p3.y))/det;
    *C = p3.y;
	
}


double _fastcall CInspectFunc::square(double init, double x)
{
	return init + x*x;
}

double _fastcall CInspectFunc::cubic(double init, double x)
{
	return init + x*x*x;
}

double _fastcall CInspectFunc::forth_power(double init, double x)
{
	return init+x*x*x*x;
}

// Yi = a0 + a1*Xi + a2*Xi^2 + ei 으로, sum(e_i)를 최소화시키게끔

// regression 된 변수를 찾도록 정규방정식을 편미분해서 0 이 되는

// 각각의 bi를 구한다
bool _fastcall CInspectFunc::get2ndOrderRegression(std::vector<double>* srcX, std::vector<double>* srcY, double *a0, double *a1, double* a2)
{
	double Y = std::accumulate(srcY->begin(), srcY->end(), 0.0);

	double X = std::accumulate(srcX->begin(), srcX->end(), 0.0);

	double X2 = std::accumulate(srcX->begin(), srcX->end(), 0.0, square);

	double X3 = std::accumulate(srcX->begin(), srcX->end(), 0.0, cubic);

	double X4 = std::accumulate(srcX->begin(), srcX->end(), 0.0, forth_power);

	double K = 0.0;

	double L = 0.0;

	int i = 0;

	int n = (int)srcX->size();

	for(i = 0; i<n; i++)
	{

		K += ((*srcY)[i]*(*srcX)[i]*(*srcX)[i]);

		L += ((*srcY)[i]*(*srcX)[i]);

	}

	double denominator = -n*X4*X2 + X4*X*X + X2*X2*X2 + X3*X3*n - 2*X3*X*X2;

	double a0p = -(Y*X4*X2 - Y*X3*X3 - X*L*X4 + X*X3*K - X2*X2*K + X2*X3*L);

	double a1p = X*Y*X4 - X*K*X2 - L*n*X4 + X3*n*K - Y*X2*X3 + X2*X2*L;

	double a2p = -(K*n*X2 - K*X*X - X2*X2*Y - X3*n*L + X3*X*Y + X*X2*L);

	*a0 = a0p/denominator;

	*a1 = a1p/denominator;

	*a2 = a2p/denominator;

	return true;
}

void _fastcall CInspectFunc::GetAveBrightValue(LPBYTE pSrc, CRect rcArea, int nJump, double &dAveVal)
{
	int left = rcArea.left; int right = rcArea.right;
	int top = rcArea.top;   int bottom = rcArea.bottom;
	int nPitch = right-left;

	double val = 0;
	int nCount= 0;//(right-left)*(bottom-top);

	for(int j=top; j<bottom; j+=nJump)
	{
		for(int i=left; i<right; i+=nJump)
		{
			if(j >= bottom) break;
			if(i >= right) break;
			val += *(pSrc + nPitch*j+i);
			nCount++;
		}
	}
	dAveVal = val / nCount;
	
}

void _fastcall CInspectFunc::CheckBlackAndWhite(LPBYTE pSrc, LPBYTE pDst, CRect rcArea)
{
	int left,right,top,bottom;
	int pitch;

	left = rcArea.left;
	right = rcArea.right;
	top = rcArea.top;
	bottom = rcArea.bottom;

	int nwidth = right-left;
	int nheight = bottom-top;

	pitch = nwidth;

	double dAvgBrVal=0;
	GausianMeanFilter(pSrc, pDst, CRect(0,0,nwidth,nheight));
	GetAveBrightValue(pSrc, CRect(0,0,nwidth,nheight), 2, dAvgBrVal);

	//black 이미지 따로
	LPBYTE pImgB = new BYTE[nwidth*nheight];
	memset(pImgB, 0, sizeof(BYTE)*nwidth*nheight);
	//white 이미지 따로
	LPBYTE pImgW = new BYTE[nwidth*nheight];
	memset(pImgW, 0, sizeof(BYTE)*nwidth*nheight);

	int nBrightGap = 10;

	for(int j=0; j<nheight; j++)
	{
		for(int i=0; i<nwidth; i++)
		{
			int nVal = *(pSrc+pitch*j+i);

			if((nVal-dAvgBrVal) < 0 && (nVal-dAvgBrVal) < -nBrightGap)
				*(pImgB+pitch*j+i) = 255; // 100
			else if((nVal - dAvgBrVal) > 0 && (nVal - dAvgBrVal) > nBrightGap)
				*(pImgW+pitch*j+i) = 255;
			else
			    *(pDst+pitch*j+i) = 0;
		}
	}
	
	LPBYTE pImgBDst = new BYTE[nwidth*nheight];
	memset(pImgBDst, 0, sizeof(BYTE)*nwidth*nheight);

	LPBYTE pImgWDst = new BYTE[nwidth*nheight];
	memset(pImgWDst, 0, sizeof(BYTE)*nwidth*nheight);

	Erode(pImgB, pImgBDst, CRect(0,0,nwidth,nheight), pitch, 3,3);
	Erode(pImgW, pImgWDst, CRect(0,0,nwidth,nheight), pitch, 3,3);
	Dilate(pImgBDst, pImgB, CRect(0,0,nwidth,nheight), pitch, 3,3);
	Dilate(pImgWDst, pImgW, CRect(0,0,nwidth,nheight), pitch, 3,3);

	for(int j=0; j<nheight; j++)
	{
		for(int i=0; i<nwidth; i++)
		{
			BYTE bBlack = (*(pImgB+pitch*j+i)==255) ? *(pImgB+pitch*j+i)-155 : 0;
			*(pDst+pitch*j+i) = *(pImgW+pitch*j+i) + bBlack;
		}
	}

	delete[] pImgW;
	delete[] pImgB;
	delete[] pImgBDst;
	delete[] pImgWDst;

}

// affine matrix 수평 [1,f,0]
//                    [0,1,0]
//                    [0,0,1]
void _fastcall CInspectFunc::AffineHorizantal(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, float f)
{
	int left,right,top,bottom;
	int pitch;

	left = rcArea.left;
	right = rcArea.right;
	top = rcArea.top;
	bottom = rcArea.bottom;

	int nwidth = right-left;
	int nheight = bottom-top;

	pitch = nwidth;

	int newx, newy;
	for(int j=0; j<nheight; j++)
	{
		for(int i=0; i<nwidth; i++)
		{
			newx = i;
			newy = (int)(f * i +j);
			if(newx<0 || newy < 0) continue;
			if(newx>=nwidth || newy >= nheight) continue;

			*(pDst+pitch*j+i) = *(pSrc+pitch*newy+newx);

		}
	}
}
// affine matrix 수직 [1,0,0]
//                    [f,1,0]
//                    [0,0,1]
void _fastcall CInspectFunc::AffineVerical(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, float f)
{
	int left,right,top,bottom;
	int pitch;

	left = rcArea.left;
	right = rcArea.right;
	top = rcArea.top;
	bottom = rcArea.bottom;

	int nwidth = right-left;
	int nheight = bottom-top;

	pitch = nwidth;

	int newx, newy;
	for(int j=0; j<nheight; j++)
	{
		for(int i=0; i<nwidth; i++)
		{
			newx = int(i + f * j);
			newy = j;
			if(newx<0 || newy < 0) continue;
			if(newx>=nwidth || newy >= nheight) continue;

			*(pDst+pitch*j+i) = *(pSrc+pitch*newy+newx);

		}
	}
}

void _fastcall CInspectFunc::GammaCorrection(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, double dGamma)
{
	int left,right,top,bottom;
	int pitch;

	left = rcArea.left;
	right = rcArea.right;
	top = rcArea.top;
	bottom = rcArea.bottom;

	int nwidth = right-left;
	int nheight = bottom-top;

	pitch = nwidth;

	double dInvGamma = 1/dGamma;

	BYTE srcB=0;
	BYTE dstB=0;
	const int C=1;
	for(int j=0; j<nheight; j++)
	{
		for(int i=0; i<nwidth; i++)
		{
			srcB = *(pSrc+pitch*j+i);
			dstB = C * (pow((double)(srcB / 255.f), dGamma)) * 255 + 0.5;
			
			if(dstB < 0) 
				dstB=0;
			if(dstB > 255) 
				dstB = 255;

			*(pDst+pitch*j+i) = dstB;
		}
	}
}

void _fastcall CInspectFunc::GetBrightMinMax(LPBYTE pSrc, CRect rcArea, int &nMinBright, int &nMaxBright)
{
	int left = rcArea.left; int right = rcArea.right;
	int top = rcArea.top;   int bottom = rcArea.bottom;
	int nPitch = right-left;

	int nMin = 255;
	int nMax = 0;

	for(int j=top; j<bottom; j++)
	{
		for(int i=left; i<right; i++)
		{
			if(*(pSrc + nPitch*j+i) > nMax)
				nMax = *(pSrc + nPitch*j+i);
			if(*(pSrc + nPitch*j+i) < nMin)
				nMin = *(pSrc + nPitch*j+i);	
		}
	}
	nMinBright = nMin;
	nMaxBright = nMax;
}

void _fastcall CInspectFunc::ContrastStrech(LPBYTE pSrc, LPBYTE pDst, CRect rcArea)
{
	int left = rcArea.left; int right = rcArea.right;
	int top = rcArea.top;   int bottom = rcArea.bottom;
	int pitch = right-left;

	
	int nMin, nMax;
	double dAvg=0;
	GetBrightMinMax(pSrc, CRect(left,top,right,bottom),nMin, nMax);
	
	// 노이즈가 스트레칭에 영향없도록 하기 위함
	if(nMin == 0 || nMax == 255)
	{
		int nhist[256]={0,};
		int nMinVal=255;
		int nMaxVal = 0;

		for(int i=rcArea.left; i<rcArea.right; i++)
		{
			for(int j=rcArea.top; j<rcArea.bottom; j++)
			{
				nhist[*(pSrc + pitch * j + i)]++;
			}
		}

		for(int g=0; g<256; g++)
		{
			if(nhist[g]!=0)
			{
				if(g < nMinVal)
					nMinVal = g;
				if(g > nMaxVal)
					nMaxVal = g;
			}
		}

		nMin = (int)(nMinVal + (nMaxVal*0.05));
		nMax = (int)(nMaxVal - (nMaxVal*0.05));

	}

	double dmult = 255.f / (nMax-nMin);

	int grayval=0;
	for(int j=top; j<bottom; j++)
	{
		for(int i=left; i<right; i++)
		{
			grayval = *(pSrc + pitch*j+i) - nMin;
			grayval *= dmult;
			if(grayval < 0) grayval = 0;
			if(grayval > 255) grayval = 255;
			*(pDst+pitch*j+i) = grayval;		
		}
	}
}

void _fastcall CInspectFunc::UnSharpMasking(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int _k)
{
	int left = rcArea.left; int right = rcArea.right;
	int top = rcArea.top;   int bottom = rcArea.bottom;
	int pitch = right-left;
	int nwidth = right-left;
	int nheight = bottom-top;

	LPBYTE pMaskImg = new BYTE[nwidth*nheight];
	memset(pMaskImg, 0, sizeof(BYTE)*nwidth*nheight);

	GausianMeanFilter(pSrc, pDst, rcArea);

	for(int y=top; y<bottom; y++)
	{
		for(int x=left; x<right; x++)
		{
			*(pMaskImg+pitch*y+x) = *(pSrc+pitch*y+x) - *(pDst+pitch*y+x);
		}
	}

	int nVal=0;
	memset(pDst, 0, sizeof(BYTE)*nwidth*nheight);

	for(int y=top; y<bottom; y++)
	{
		for(int x=left; x<right; x++)
		{
			nVal = *(pMaskImg+pitch*y+x);
			*(pDst+pitch*y+x) = *(pSrc+pitch*y+x) + _k*nVal;

		}
	}

	delete[] pMaskImg;
}

void _fastcall CInspectFunc::scramble(PST_FURIEDATA p,  int w, int h, int nChoose, int nselect)
{
	int nWidth = w, nHeight = h;
	int nLength;

	if(nChoose == 1) nLength = nWidth;
	if(nChoose == 2) nLength = nHeight;

	int j=0, m=0;
	double dtemp=0;
	for(int i=0; i<nLength; i++)
	{
		if((i>j && nselect==1) || (i<j && nselect == -1))
		{
			dtemp = p[j].d_Re;
			p[j].d_Re = p[i].d_Re;
			p[i].d_Re = dtemp;

			dtemp = p[j].d_Im;
			p[j].d_Im = p[i].d_Im;
			p[i].d_Im = dtemp;
		}
		m = nLength>>1;
		while((j>=m)&(m>=2))
		{
			j -= m;
			m = m>>1;
		}
		j+=m;
	}
}

void _fastcall CInspectFunc::buterfly(PST_FURIEDATA p,  int w, int h, int nChoose, int nselect)
{
	int N, half_N;
	register int nOffset;
	register int i;
	int nLength, nWidth = w, nHeight = h;
	double temp_re, temp_im;
	double angle, wtemp, w_re, w_im, wp_re, wp_im;
	int j=0;

	if(nChoose == 1) nLength = nWidth;
	if(nChoose == 2) nLength = nHeight;

	N=1;
	
	int nLogVal=0;
	if(nChoose == 1)
	{
		nLogVal = LogB(w,2);
	}
	else
	{
		nLogVal = LogB(h,2);
	}
	

	
	for(int k=0; k<nLogVal; k++) // comment kjw: k 는 이미지 크기의 승수 2^7
	{
		half_N = N;
		N<<=1;
		if(nselect == 1) angle = -2.0 * 3.141592 / N;
		if(nselect == -1) angle = -2.0 * 3.141592 / N * (-1);
		wtemp = sin(0.5 * angle);
		wp_re = -2.0 * wtemp * wtemp;//-(1-cos(angle));//
		wp_im = sin(angle);
		w_re = 1.0;
		w_im = 0.;

		for(nOffset=0; nOffset<half_N; nOffset++)
		{
			for(i=nOffset; i<nLength; i+=N)
			{
				j = i+half_N;
				temp_re = (w_re * p[j].d_Re) - (w_im * p[j].d_Im);
				temp_im = (w_im * p[j].d_Re) + (w_re * p[j].d_Im);
				p[j].d_Re = p[i].d_Re - temp_re;
				p[i].d_Re += temp_re;
				p[j].d_Im = p[i].d_Im - temp_im;
				p[i].d_Im+= temp_im;
			}
			wtemp = w_re;
			w_re = wtemp * wp_re - w_im * wp_im + w_re;
			w_im = w_im * wp_re + wtemp * wp_im + w_im;
		}
	}
	if(nselect == INVERSE)
	{
		for(int i=0; i<nLength; i++)
		{
			p[i].d_Re /= nLength;
			p[i].d_Im /= nLength;
		}
	}

}

void _fastcall CInspectFunc::FFT(PST_FURIEDATA p, int w, int h, int nChoose, int nSelect)
{
	scramble(p, w, h, nChoose, nSelect);
	buterfly(p, w, h, nChoose, nSelect);

}

int _fastcall CInspectFunc::GetFFTSize(int n)
{
	return (int)pow(2.0, int(log((double)(n))/log(2.0) + 0.9999));
}

void _fastcall CInspectFunc::Furie2DTransForm(LPBYTE pSrc, LPBYTE pDst, PPST_FURIE_VAL ppRes, CRect rcArea, int nInverseMode, BOOL bShow)
{
	int left = rcArea.left; int right = rcArea.right;
	int top = rcArea.top;   int bottom = rcArea.bottom;
	int pitch = right-left;
	int nwidth = right-left;
	int nheight = bottom-top;

	int nSelect = nInverseMode;
	int ntemp = 0;

	memset(g_ptagRowData, 0 , sizeof(__ST_FURIE_DATA)*nwidth);
	memset(g_ptagColData, 0 , sizeof(__ST_FURIE_DATA)*nheight);
	//PST_FURIEDATA ptagRowData = new __ST_FURIE_DATA[nwidth];
	//PST_FURIEDATA ptagColData = new __ST_FURIE_DATA[nheight];

	
	if(nSelect == NONE)
	{
		for(int y=0; y < nheight; y++)
		{
			for(int x=0; x<nwidth; x++)
			{
				ppRes[y][x].dRes_Re = (double)(*(pSrc+pitch*(y+top)+(x+left))); 
				ppRes[y][x].dRes_Im = 0.;
			}
		}
	}

	if(nSelect == INVERSE) 
	{
		for(int y=0; y < nheight; y++)
		{
			for(int x=0; x<nwidth; x++)
			{
				ppRes[y][x].dRes_Re = 255; 
				ppRes[y][x].dRes_Im = 0.;
			}
		}
		//CInspectFunc::LowPassFilter(ppRes, nwidth, nheight, 20);
		//CInspectFunc::GausianLowPassFilter(ppRes, nwidth, nheight, 100, false);
		//CInspectFunc::GausseHighPassFilter(ppRes, nwidth, nheight, 0.1);
		NotchFilter(ppRes, nwidth, nheight, 30);

		for(int y=0; y < nheight; y++)
		{
			for(int x=0; x<nwidth; x++)
			{
				*(pDst+pitch*(y+top)+(x+left)) = (BYTE)ppRes[y][x].dRes_Re;
			}
		}
	    return;
	}
	/************************************************************************/
	/* 이미지 -> 주파수 스펙트럼 변환                                       */
	/************************************************************************/
	// 열 변환
	for(int y=0; y<nheight; y++)
	{
		for(int x=0; x<nwidth; x++)
		{
			
			g_ptagColData[x].d_Re = ppRes[y][x].dRes_Re;
			g_ptagColData[x].d_Im = ppRes[y][x].dRes_Im;
		}
		
		FFT(g_ptagColData, GetFFTSize(nwidth), GetFFTSize(nheight), 1, nSelect);

		for(int x=0; x<nwidth; x++)
		{
			ppRes[y][x].dRes_Re = g_ptagColData[x].d_Re;
			ppRes[y][x].dRes_Im = g_ptagColData[x].d_Im;
			
		}
	}

	// 행 변환
	for(int x=0; x<nwidth; x++)
	{
		for(int y=0; y<nheight; y++)
		{
			g_ptagRowData[y].d_Re= ppRes[y][x].dRes_Re;
			g_ptagRowData[y].d_Im= ppRes[y][x].dRes_Im;
		}
	
		FFT(g_ptagRowData, GetFFTSize(nwidth), GetFFTSize(nheight), 2, nSelect);

		for(int y=0; y<nheight; y++)
		{
			ppRes[y][x].dRes_Re = g_ptagRowData[y].d_Re;
			ppRes[y][x].dRes_Im = g_ptagRowData[y].d_Im;
		}
	}
	
	

	if(bShow)
	{
		// 구해진 스펙트럼을 영상으로 변환
		double dmaxV = 0;
		double drange;
		
		int ntransx=0;
		int ntransy=0;
		double dval=0;

		if(nSelect == NONE)
		{
			for(int y=0; y < nheight; y++)
			{
				for(int x=0; x<nwidth; x++)
				{
					dval = sqrt(pow(ppRes[y][x].dRes_Re, 2) + pow(ppRes[y][x].dRes_Im,2));
					ppRes[y][x].dRes = dval; 
					if(ppRes[y][x].dRes > dmaxV)  dmaxV = ppRes[y][x].dRes;

				}
			}

			drange = 255.0 / log(1.0 + dmaxV);

			for(int y=0; y < nheight; y++)
			{
				for(int x=0; x<nwidth; x++)
				{
					ntemp = (int)(log(1.0+ppRes[y][x].dRes)*drange+0.5);

					if(ntemp < 0) 
						*(pDst+pitch*(y+top)+(x+left)) = 0;
					else if(ntemp > 255) 
						*(pDst+pitch*(y+top)+(x+left)) = 255;
					else
						*(pDst+pitch*(y+top)+(x+left)) = (BYTE)ntemp;

				}
			}

		}
	}
	
	
	if(nSelect == INVERSE)
	{
		// 스펙트럼 -> 이미지화
		for(int y=0; y < nheight; y++)
		{
			for(int x=0; x<nwidth; x++)
			{
				ntemp = (int)ppRes[y][x].dRes_Re;

				if(ntemp > 255) ntemp = 255;
				if(ntemp < 0) ntemp = 0;

				BYTE temp = (BYTE)ntemp;
	            
				*(pDst+pitch*(y+top)+(x+left)) = temp;

			}
		}
		
	}

}

BOOL _fastcall CInspectFunc::IsPowerOf2(int n)
{
	int ref = 1;
	while(ref < n)
		ref <<=1;
	if(ref == n)
		return TRUE;
	else
		return FALSE;
}

void _fastcall CInspectFunc::LowPassFilter(PPST_FURIE_VAL ppImg, int W, int H, int nCutOff)
{
	register int i, j;

	int w = W / 2;
	int h = H / 2;
	int x, y;
	for(j=0; j<H; j++)
	{
		for(i=0; i<W; i++)
		{
			x = i + w;
			y = j + h;

			if(x >= W) x -= W;
			if(y >= H) y -= H;

			if(sqrt((double)(x-w)*(x-w) + (y-h)*(y-h)) > nCutOff)
				ppImg[j][i].dRes_Re = ppImg[j][i].dRes_Im = 0;
		}
	}
}

void _fastcall CInspectFunc::GausianLowPassFilter(PPST_FURIE_VAL ppImg, int W, int H, double dD0, BOOL bShift)
{
	// 각 모서리에서 중심으로 갈수록 어두워지는 방식
	// shift 영상으로 본다면 중앙은 밝고 사이드로 가면서 어두워지는 방식
	register int i, j;
	int cx = W/2;
	int cy = H/2;
	int x, y;

	double hVal, dD2;

	if(bShift)
	{
		for(j=0; j<H; j++)
		{
			for(i=0; i<W; i++)
			{	
				dD2 = (cx-i)*(cx-i)+(cy-j)*(cy-j);
				hVal = exp(-1*dD2/(2*dD0*dD0));

				ppImg[j][i].dRes_Re *= hVal;
				ppImg[j][i].dRes_Im *= hVal;
			}
		}
	}
	else
	{
		for(j=0; j<H; j++)
		{
			for(i=0; i<W; i++)
			{
				x = i + cx;
				y = j + cy;

				if(x >= W) x -= W;
				if(y >= H) y -= H;
				dD2 = (cx-x)*(cx-x)+(cy-y)*(cy-y);
				hVal = exp(-1*dD2/(2*dD0*dD0));

				ppImg[j][i].dRes_Re *= hVal;
				ppImg[j][i].dRes_Im *= hVal;

			}
		}
	}
	
}
//Gauss 고역통과 필터
void _fastcall CInspectFunc::GausseHighPassFilter(PPST_FURIE_VAL ppImg, int W, int H, double dCutOff, BOOL bShift)
{
	register int i, j;

	int cx = W/2;
	int cy = H/2;
	int x, y;
	double dist2, hval;

	if(bShift)
	{
		for(j=0; j<H; j++)
		{
			for(i=0; i<W; i++)
			{
				dist2 = (double)(cx-i)*(cx-i) + (cy-j)*(cy-j);
				hval = 1.0 -exp(-dist2 / (2*dCutOff*dCutOff));


				ppImg[j][i].dRes_Re *= hval;
				ppImg[j][i].dRes_Im *= hval;
			}
		}
	}
	else
	{
		for(j=0; j<H; j++)
		{
			for(i=0; i<W; i++)
			{
				x = i + cx;
				y = j + cy;

				if(x >= W) x -= W;
				if(y >= H) y -= H;
				dist2 = (double)(cx-x)*(cx-x) + (cy-y)*(cy-y);
				//dist2 = (double)(x-w)*(x-w) + (y-h)*(y-h);
				hval = 1.0 -exp(-dist2 / (2*dCutOff*dCutOff));


				ppImg[j][i].dRes_Re *= hval;
				ppImg[j][i].dRes_Im *= hval;
			}
		}
	}
	
}

void _fastcall CInspectFunc::NotchFilter(PPST_FURIE_VAL ppImg, int W, int H, double dCutOff)
{
	register int i, j;
	int cx = W/4;
	int cy = H/4;
	int x, y = 0;
	double dist2, hval;

	//leftTop
	for(j=0; j<H/2; j++)
	{
		for(i=0; i<W/2; i++)
		{
			dist2 = (double)(cx-i)*(cx-i) + (cy-j)*(cy-j);
			hval = 1.0 -exp(-dist2 / (2*dCutOff*dCutOff));


			ppImg[j][i].dRes_Re *= hval;
			ppImg[j][i].dRes_Im *= hval;
		}
	}
	//leftBottom
	cx = W/4;
	cy = 3*H/4;
	for(j=H/2; j<H; j++)
	{
		for(i=0; i<W/2; i++)
		{
			dist2 = (double)(cx-i)*(cx-i) + (cy-j)*(cy-j);
			hval = 1.0 -exp(-dist2 / (2*dCutOff*dCutOff));


			ppImg[j][i].dRes_Re *= hval;
			ppImg[j][i].dRes_Im *= hval;
		}
	}

	cx = 3*W/4;
	cy = H/4;
	for(j=0; j<H/2; j++)
	{
		for(i=W/2; i<W; i++)
		{
			dist2 = (double)(cx-i)*(cx-i) + (cy-j)*(cy-j);
			hval = 1.0 -exp(-dist2 / (2*dCutOff*dCutOff));


			ppImg[j][i].dRes_Re *= hval;
			ppImg[j][i].dRes_Im *= hval;
		}
	}

	cx = 3*W/4;
	cy = 3*H/4;
	for(j=H/2; j<H; j++)
	{
		for(i=W/2; i<W; i++)
		{
			dist2 = (double)(cx-i)*(cx-i) + (cy-j)*(cy-j);
			hval = 1.0 -exp(-dist2 / (2*dCutOff*dCutOff));


			ppImg[j][i].dRes_Re *= hval;
			ppImg[j][i].dRes_Im *= hval;
		}
	}

	

}