#include <math.h>
#include "Tracklets.h"

#define Assign1(a,b,c,d) (pPoint1[0]=a,pPoint1[1]=b,pPoint1[2]=c,pPoint1[3]=d);
#define Assign2(a,b,c,d) (pPoint2[0]=a,pPoint2[1]=b,pPoint2[2]=c,pPoint2[3]=d);

struct MatArea
{
	CvPoint LeftUp;
	CvPoint RightUp;
	CvPoint RightDown;
	CvPoint LeftDown;
	CvPoint Center;
};

struct Relationship
{
	bool connect;
	MatArea area1;
	MatArea area2;
};

class CW_MotionOperator
{
private:
	int step1;
	int step2;
	int cam_num;
	Relationship* rl;

public:
	void initial(int camnum);
	bool load_areas(FILE* areafile);
	float dist_motion(Tracklet_Info* tracklet1,Tracklet_Info* tracklet2,float* distance);
	void release();
};