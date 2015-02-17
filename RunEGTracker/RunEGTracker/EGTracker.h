#include <stdio.h>
#include <string.h>
#include "SVBasic.h"
#include "APPInterface.h"
#include "DAT_CS2.h"
#include "Tracklets.h"
#include "FEAT_PMCSHR.h"
#include "FEAT_motion.h"

//ensure your hardware meets these requires
#define MAX_GroundTruthNum 100000	//maximum of the total boundingboxes in all videos
#define MAX_CamNum 10				//maximum of the camera number
#define MAX_TRACKLET 3000			//maximum of the total tracklet number
#define MAX_TRACKLET_LENGTH 10000	//max length of tracklets
#define MAX_PATH_LENGTH 100			//max length of the file path

struct Boundingbox
{
	int nCam;
	int nFrm;
	int ObjID;
	CvRect r;
};

class CW_EGTracker
{
public: 
	// interface
	bool LoadData(const char* videopath, const char* gthpath, const char* topopath,const char* areapath);
	bool LoadParam(const char* path);
	void Initial(void* param);
	void Process(void* param);
	void GetResult(Boundingbox* output, int* length);
	void Release();

private:
	//Initial
	int c_camnum;
	char** c_path;
	int c_gthnum;
	Boundingbox* c_groundtruth;
	int c_trkletnum;
	Tracklet_Info* c_tracklet;

	CW_MotionOperator c_motion;
	CW_PMCSHROperator c_pmcshr;
	CW_DATbyCS2 c_cs2;

	//Params
	int c_maxflow;
	int c_waitframe;
	float c_mugpairthrd;
	float c_AIFconvthrd;

	//Gauss Background
	int c_bgupdatespeed;
	int c_bggaussthrd;
	int c_fgthrd;

	//read groundtruth
	int ReadGroundtruth(FILE * inputfile,int cam);//obtain c_groundtruth
	int ReadVideo(FILE * datasetfile);//obtain c_path
	void GetTopology(FILE* topofile);//obtain c_topo
	void GetArea(FILE* areafile);//obtain c_area

	//tracklet processing
	void GetTrklet();//obtain c_tracklet and c_trkletnum
	void GetTrkletbyFrm(IplImage* frame,IplImage* fgimage,APPSVInterface** ppBaseModule,int* tracker_idx,Boundingbox* cur_output,int cur_output_num, int nFrm, int Cam);
	void NewTrklet(IplImage* frame,IplImage* fgimage, CvRect r, int ID, int nFrm, int Cam);
	void AddtoTrklet(IplImage* frame,IplImage* fgimage, CvRect r, int ID, int nFrm, int Cam, float conv);
	void NewTracker(IplImage* frame,APPSVInterface** ppBaseModule,int* tracker_idx,CvRect r,int ID,int new_idx);//new_idx dicides the position of trackers
	void EndTrklet(int ID);
	void SelectTrklet(int length_thrd);
	void CpyTrklet(Tracklet_Info* tracklet1,Tracklet_Info* tracklet2);
	void DelTrklet(Tracklet_Info* tracklet2);
	void OrderTrklet();
	void SwitchTrklet(Tracklet_Info* tracklet1,Tracklet_Info* tracklet2);

	//Graph Preparation
	bool* c_topo;
	float c_sim[10];
	void Compute_Sims(int waitframe);

	//Output data for test
	void TrkletToOutput4Test(char* path,int waitframe);

	//DataAssociation
	void OutputToDMX(char* dmx_path, int* n, int* e, int maxf,int waitframe);

	//Update Tracklet Results
	void UpdateTrklet(int* links);
	int FindPreTracklet(int i,int* result);

	//others
	void GetFG(IplImage* src1, IplImage* src2, IplImage* dst);
	void pyrfilter(IplImage* src,IplImage* dst);
	double TwoRectOverlap(CvRect r1,CvRect r2,int flag);
};
