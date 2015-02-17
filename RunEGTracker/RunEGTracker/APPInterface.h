#if !defined(AFX_1100816_H__DF6FB001_509B_11D6_8E99_00E04C498874__INCLUDED_)
#define AFX_1100816_H__DF6FB001_509B_11D6_8E99_00E04C498874__INCLUDED_

#include <cv.h>

#define MAX_TRACKER 100

#define EXPORT
#ifdef EXPORT
#define APPDLL _declspec(dllexport)
#else
#define APPDLL _declspec(dllimport)
#endif 

typedef enum MAIN_ALG_TYPE
{
	ALG_ALARM_AUTOTRACK = 5,    //auto tracking
	ALG_ALARM_CWFGTRACK = 6,    //raw tracking
	ALG_ALARM_IIFTRACK = 7,    //AIF tracker
}MAIN_ALG_TYPE;

typedef enum INTERACTION_TYPE
{
	INTERACTION_InputPt=0,
}INTERACTION_TYPE;

class APPSVInterface
{
public:
	virtual void Initial(void* param) = 0;
	virtual bool LoadParamFromFile(const char* szFilePath) = 0;
	virtual bool SaveParamToFile(const char* szFilePath) = 0;
	virtual bool SetParameter(const unsigned char* lpbyBufIPara, int nBufLen) = 0;
	virtual bool GetParameter(void* lpbyBufIPara) = 0;
	virtual void Process(const unsigned char* lpbyBufData ) = 0;
	virtual void GetResult(void* OutData) = 0;
	virtual void Release() = 0;
	virtual void Draw(const unsigned char* lpbyBufData) = 0;
	virtual void InterAction(INTERACTION_TYPE type,const unsigned char* lpbyBufData) = 0;
};
//every dll must have a gloal function FXCreateObject(IXBaseVSModule** ppBaseModule)

struct InPutAPP_IIFTracker
{
	int nNum;
	IplImage* img;
	CvRect r;
};

struct OutPutAPP_IIFTracker
{
	CvRect r;
	float conv;
};

extern "C" bool APPDLL APPSVCreateObject_IIFTrack(APPSVInterface** ppBaseModule, MAIN_ALG_TYPE AlgType);

#endif
