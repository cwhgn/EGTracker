#include "Tracklets.h"

#define MAX_CLUSTERS_mcsh 30//30 classes
#define MCSH_COLOR_DISTANCE_THRE 0.06//the distance threshold bewteen colors

class CW_PMCSHROperator
{
public: 
	//extract features
	void cal_img_mcsh(IplImage* frame,IplImage* frame_bk, CvRect rect, mcsh * mcsh_region);
	void cal_tl_pmcsh(Tracklet_Info* tl_i);
	void cal_tl_mcsh(Tracklet_Info* tl_i);

	//compute similarity
	float dist_feature( mcsh * h1, mcsh * h2 );//two msch
	float dist_mug_mcsh(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2);//mug between two sets of mcsh
	float dist_mug_tl(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2);//mug between two pmcsh (tracklets)
	float dist_avg_tl(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2);//average between two pmcsh (tracklets)

	//select suitable pair
	float dist_mug(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2);

private:
	bool  del_color(int r, int g,int b);//delete bad pixels
	float MCSH_Color_Distance(float p1[MAX_NUM_mcsh][4],float p2[MAX_NUM_mcsh][4], int i,int j);
};