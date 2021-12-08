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
	// ���� -> �ܺ� ����
	INSIDE_LEFT_TOP=0,
	INSIDE_RIGHT_TOP,
	INSIDE_RIGHT_BOTTOM,
	INSIDE_LEFT_BOTTOM,
	
	// �ܺ� -> ���� ����
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
	// ����ȭ �� �������� 
	//void _fastcall Binary(LPBYTE pSrc, LPBYTE pDst, CRoiClass* roi, int nTh, int nMode);
	void _fastcall Binary(LPBYTE pSrc, LPBYTE pDst, CRect& rect, int nTh, int nMode);
	static void _fastcall Binary_th(LPBYTE pSrc, LPBYTE pDst, CRoiClass* roi, int nTh, int nMode);
	int _fastcall Erode(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int nSizeX, int nSizeY);
	int _fastcall Dilate(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int nSizeX, int nSizeY);

	// �������� 
	void _fastcall GausianMeanFilter(LPBYTE pSrc, LPBYTE pDst, CRect rcArea);
	void _fastcall MidFilter(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int _objectNum);
	void _fastcall SobelSharpness(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int direction);
	void _fastcall LaplacianEdge(LPBYTE pSrc, int* pDst, CRect rcArea);
	// ���� ��ȯ
	void _fastcall GammaCorrection(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, double dGamma);
	// ��Ʈ��Ʈ ��Ʈ��ġ ó��
	void _fastcall ContrastStrech(LPBYTE pSrc, LPBYTE pDst, CRect rcArea);
	void _fastcall UnSharpMasking(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int _k);
	
	// ����& �������ð��� 
	void _fastcall FindEdgeLine(LPBYTE pSrc, CRect rcArea, int N, int _direction, int _direction2, int nSlope, double* Alpha, double* Beta, double* C);
	int _fastcall LineSlopeToEdgeX_LeftToRight(LPBYTE pSrc, CRect rcArea, int nposy, double* _data, int* dEdgePos, double dSlope);
	int _fastcall LineSlopeToEdgeX_RightToLeft(LPBYTE pSrc, CRect rcArea, int nposy, double* _data, int* dEdgePos, double dSlope);
	int _fastcall LineSlopeToEdgeY_TopToBottom(LPBYTE pSrc, CRect rcArea, int nposx, double* _data, int* dEdgePos, double dSlope);
	int _fastcall LineSlopeToEdgeY_BottomToTop(LPBYTE pSrc, CRect rcArea, int nposx, double* _data, int* dEdgePos, double dSlope);
	BOOL _fastcall LineFitting(int N, int removeCnt, int _dir, double* x, double* y, double* Alpha, double* Beta, double* C);
	int _fastcall ChooseGoodData_FirstBase(int N, double *x, double *y, int space);
	int _fastcall FindCrossPoint(double t1, double a1, double b1, double t2, double a2, double b2, double *cx, double *cy);
	double _fastcall LineToPointDist(double Alpha, double Beta, double C, int x, int y, int direction);
	
	
	// ������ȯ ��������
	void _fastcall HoughTransLineFitting(LPBYTE pSrc, CRoiClass* roi, int _direction, int _direction2, double _slope, double* A, double* B);
	// ������ȯ �� ����
	void _fastcall HoughTransFindElipse(LPBYTE pSrc, CRect rcArea, int& centerx, int& centery, int& Radius);
	// ���ó��
	void _fastcall FindBlob(LPBYTE pSrc, CRoiClass* roi, int _th, int _mode, std::vector<CRect>& v);
	
	// ���� ȸ��
	void _fastcall RotateImage(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, int pitch, int centerx, int centery, double dDegree);
	// Seta �缱�� ����
	void _fastcall BiLinearInterpolation(LPBYTE pSrc, LPBYTE pDst, int nPitch, CRect rcArea, int x, int y, double newx, double newy);
	// ���� �Ǻ� �˰���
	BOOL _fastcall CnCheckArea(POINT poly[], int nCount, int nX, int nY); 
	// Ư¡�� ����
	void _fastcall InspectHerrisCorner(LPBYTE fm, CRect rcArea, std::vector<CPoint>& v);
	
	// � ���� �׽�Ʈ(�̿ϼ�)
	void _fastcall FindRoundLine(LPBYTE pSrc, CRect rcArea, int _dir, double &A0, double &A1, double &A2);
	void _fastcall CalcThreePoint(CPoint p1, CPoint p2, CPoint p3, double *Alpha, double *Beta, double *C);
	// � ������
	bool _fastcall get2ndOrderRegression(std::vector<double>* srcX, std::vector<double>* srcY, double *a0, double *a1, double* a2);
	double _fastcall square(double init, double x);
	double _fastcall cubic(double init, double x);
	double _fastcall forth_power(double init, double x);

	void _fastcall GetAveBrightValue(LPBYTE pSrc, CRect rcArea, int nJump, double &dAveVal);
	void _fastcall GetBrightMinMax(LPBYTE pSrc, CRect rcArea, int &nMinBright, int &nMaxBright);

	void _fastcall CheckBlackAndWhite(LPBYTE pSrc, LPBYTE pDst, CRect rcArea);

	// ������ �⺻ ��ȯ
	void _fastcall AffineVerical(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, float f); // ����
	void _fastcall AffineHorizantal(LPBYTE pSrc, LPBYTE pDst, CRect rcArea, float f); // ����
	
	// ǻ���� ��ȯ
	void _fastcall Furie2DTransForm(LPBYTE fm, LPBYTE pDst, PPST_FURIE_VAL ppRes, CRect rcArea, int nInverseMode, BOOL bShow);
	BOOL _fastcall IsPowerOf2(int n);
	int _fastcall  GetFFTSize(int n);
	void _fastcall FFT(PST_FURIEDATA p, int w, int h, int nChoose, int nSelect);
	void _fastcall scramble(PST_FURIEDATA p,  int w, int h, int nChoose, int nselect);
	void _fastcall buterfly(PST_FURIEDATA p,  int w, int h, int nChoose, int nselect);

	// ���ļ� ����
	void _fastcall LowPassFilter(PPST_FURIE_VAL ppImg, int W, int H, int nCutOff);
	void _fastcall GausianLowPassFilter(PPST_FURIE_VAL ppImg, int W, int H, double dD0, BOOL bShift);
	void _fastcall GausseHighPassFilter(PPST_FURIE_VAL ppImg, int W, int H, double dCutOff, BOOL bShift);
	void _fastcall NotchFilter(PPST_FURIE_VAL ppImg, int W, int H, double dCutOff);

	
};

extern PST_FURIEDATA g_ptagRowData;
extern PST_FURIEDATA g_ptagColData;

#endif