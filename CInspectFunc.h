#ifndef CINSPECT_H
#define CINSPECT_H

#define OBJECT_VAL 255
#define BACKGROUND_VAL 0
#define PI 3.141592654 
#define _THETA 360

#include <vector>
#include <numeric>

#include "../CRoiClass.h"
#include "../common/CommonDefine.h"


enum{
	// 내부 -> 외부 검출
	INSIDE_LEFT_TOP=0,
	INSIDE_RIGHT_TOP,
	INSIDE_RIGHT_BOTTOM,
	INSIDE_LEFT_BOTTOM,
	
	// 외부 -> 내부 검출
	OUTSIDE_LEFT_TOP,
	OUTSIDE_RIGHT_TOP,
	OUTSIDE_RIGHT_BOTTOM,
	OUTSIDE_LEFT_BOTTOM,
	
};
enum{ //direction
	DirX=0,
	DirY,
	DirXY,
	LTR,
	RTL,
	TTB,
	BTT,
};



namespace CInspectFunc
{
	void  _fastcall CopyImage(LPBYTE psrc, LPBYTE pdst, CRoiClass* roi, int _mode = -1);
	void  _fastcall CopyImage(LPBYTE psrc, LPBYTE pdst, CRect& rect, int _mode=-1);
	// 이진화 및 모폴로지 
	//void _fastcall Binary(LPBYTE pSrc, LPBYTE pDst, CRoiClass* roi, int nTh, int nMode);
	void _fastcall Binary(LPBYTE pSrc, LPBYTE pDst, CRect& rect, int nTh, int nMode);
	static void _fastcall Binary_th(LPBYTE pSrc, LPBYTE pDst, CRoiClass* roi, int nTh, int nMode);
	int _fastcall Erode(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int nSizeX, int nSizeY);
	int _fastcall Dilate(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int nSizeX, int nSizeY);

	// 영상필터 
	void _fastcall GausianMeanFilter(LPBYTE pSrc, LPBYTE pDst, CRect rcArea);
	void _fastcall MidFilter(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int _objectNum);
	void _fastcall SobelSharpness(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int direction);
	void _fastcall LaplacianEdge(LPBYTE pSrc, int* pDst, CRect rcArea);
	// 감마 변환
	void _fastcall GammaCorrection(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, double dGamma);
	// 콘트라스트 스트레치 처리
	void _fastcall ContrastStrech(LPBYTE pSrc, LPBYTE pDst, CRect rcArea);
	void _fastcall UnSharpMasking(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int _k);
	
	// 엣지& 직선피팅관련 
	void _fastcall FindEdgeLine(LPBYTE pSrc, CRect rcArea, int N, int _direction, int _direction2, int nSlope, double* Alpha, double* Beta, double* C);
	int _fastcall LineSlopeToEdgeX_LeftToRight(LPBYTE pSrc, CRect rcArea, int nposy, double* _data, int* dEdgePos, double dSlope);
	int _fastcall LineSlopeToEdgeX_RightToLeft(LPBYTE pSrc, CRect rcArea, int nposy, double* _data, int* dEdgePos, double dSlope);
	int _fastcall LineSlopeToEdgeY_TopToBottom(LPBYTE pSrc, CRect rcArea, int nposx, double* _data, int* dEdgePos, double dSlope);
	int _fastcall LineSlopeToEdgeY_BottomToTop(LPBYTE pSrc, CRect rcArea, int nposx, double* _data, int* dEdgePos, double dSlope);
	BOOL _fastcall LineFitting(int N, int removeCnt, int _dir, double* x, double* y, double* Alpha, double* Beta, double* C);
	int _fastcall ChooseGoodData_FirstBase(int N, double *x, double *y, int space);
	int _fastcall FindCrossPoint(double t1, double a1, double b1, double t2, double a2, double b2, double *cx, double *cy);
	double _fastcall LineToPointDist(double Alpha, double Beta, double C, int x, int y, int direction);
	
	
	// 허프변환 직선검출
	void _fastcall HoughTransLineFitting(LPBYTE pSrc, CRoiClass* roi, int _direction, int _direction2, double _slope, double* A, double* B);
	// 허프변환 원 검출
	void _fastcall HoughTransFindElipse(LPBYTE pSrc, CRect rcArea, int& centerx, int& centery, int& Radius);
	// 블랍처리
	void _fastcall FindBlob(LPBYTE pSrc, CRoiClass* roi, int _th, int _mode, std::vector<CRect>& v);
	
	// 영상 회전
	void _fastcall RotateImage(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int centerx, int centery, double dDegree);
	// Seta 양선형 보간
	void _fastcall BiLinearInterpolation(LPBYTE pSrc, LPBYTE pDst, int nPitch, CRect rcArea, int x, int y, double newx, double newy);
	// 영역 판별 알고리즘
	BOOL _fastcall CnCheckArea(POINT poly[], int nCount, int nX, int nY); 
	// 특징점 추출
	void _fastcall InspectHerrisCorner(LPBYTE fm, CRect rcArea, std::vector<CPoint>& v);
	
	// 곡선 검출 테스트(미완성)
	void _fastcall FindRoundLine(LPBYTE pSrc, CRect rcArea, int _dir, double &A0, double &A1, double &A2);
	void _fastcall CalcThreePoint(CPoint p1, CPoint p2, CPoint p3, double *Alpha, double *Beta, double *C);
	// 곡선 곡률계산
	bool _fastcall get2ndOrderRegression(std::vector<double>* srcX, std::vector<double>* srcY, double *a0, double *a1, double* a2);
	double _fastcall square(double init, double x);
	double _fastcall cubic(double init, double x);
	double _fastcall forth_power(double init, double x);

	void _fastcall GetAveBrightValue(LPBYTE pSrc, CRect rcArea, int nJump, double &dAveVal);
	void _fastcall GetBrightMinMax(LPBYTE pSrc, CRect rcArea, int &nMinBright, int &nMaxBright);

	void _fastcall CheckBlackAndWhite(LPBYTE pSrc, LPBYTE pDst, CRect rcArea);

	// 어파인 기본 변환
	void _fastcall AffineVerical(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, float f); // 수직
	void _fastcall AffineHorizantal(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, float f); // 수평
	
	// 퓨리에 변환
	void _fastcall Furie2DTransForm(LPBYTE fm, LPBYTE pDst, PPST_FURIE_VAL ppRes, CRect rcArea, int nInverseMode, BOOL bShow);
	BOOL _fastcall IsPowerOf2(int n);
	int _fastcall  GetFFTSize(int n);
	void _fastcall FFT(PST_FURIEDATA p, int w, int h, int nChoose, int nSelect);
	void _fastcall scramble(PST_FURIEDATA p,  int w, int h, int nChoose, int nselect);
	void _fastcall buterfly(PST_FURIEDATA p,  int w, int h, int nChoose, int nselect);

	// 주파수 필터
	void _fastcall LowPassFilter(PPST_FURIE_VAL ppImg, int W, int H, int nCutOff);
	void _fastcall GausianLowPassFilter(PPST_FURIE_VAL ppImg, int W, int H, double dD0, BOOL bShift);
	void _fastcall GausseHighPassFilter(PPST_FURIE_VAL ppImg, int W, int H, double dCutOff, BOOL bShift);
	void _fastcall NotchFilter(PPST_FURIE_VAL ppImg, int W, int H, double dCutOff);

	
};

extern PST_FURIEDATA g_ptagRowData;
extern PST_FURIEDATA g_ptagColData;

#endif