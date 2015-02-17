#include <cv.h>
#include <vector>
using namespace std;


#if !defined(AFX_Basic_H__DF6FB001_509B_11D6_8E99_00E04C498874__INCLUDED_)
#define AFX_Basic_H__DF6FB001_509B_11D6_8E99_00E04C498874__INCLUDED_


class SVBASIC
{
public:
	virtual void Initial(void* param) = 0;
	virtual void Process(void* InData) = 0;
	virtual void GetResult(void* OutData) = 0;
	virtual void Release() = 0;
};
struct InParamSVFinalBG
{
	int _imgw;
	int _imgh;
	int _iniframe;
	int _thresh;
};
struct InDataSVFinalBG
{
	IplImage* _gray;
	int _renewflag;
};
struct OutPutSVFinalBG
{
	IplImage* _fg;
	IplImage* _bg;
	int _bgok;
};

#define EXPORT

#ifdef EXPORT
#define SVDLL _declspec(dllexport)
#else
#define SVDLL _declspec(dllimport)
#endif 

//#undef SVDLL
//#define SVDLL
extern "C" bool SVDLL SVCreateObject_SVEasyTrack(SVBASIC** ppBaseModule);
extern "C" bool SVDLL SVCreateObject_SVFinalBG(SVBASIC** ppBaseModule);
extern "C" bool SVDLL SVCreateObject_SVStaticsBG(SVBASIC** ppBaseModule);
extern "C" void _declspec(dllimport) SVThreshWave(const int* _data,int* _dst,int _lenth,int _peak);

#endif









































