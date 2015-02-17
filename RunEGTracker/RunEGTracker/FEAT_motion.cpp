#include "FEAT_motion.h"

void CW_MotionOperator::initial(int camnum)
{
	step1=25;
	step2=25;
	cam_num=camnum;
	rl=(Relationship*)malloc(camnum*camnum*sizeof(Relationship));
	memset(rl,0,camnum*camnum*sizeof(Relationship));
}

bool CW_MotionOperator::load_areas(FILE* areafile)
{
	int n=0;
	while(feof(areafile)==0)
	{
		char flag[100]="";
		fscanf(areafile,"%s",flag);
		if(!strcmp(flag,"Relationship:"))
		{
			int cam1=0;
			int cam2=0;
			fscanf(areafile,"%d %d", &cam1, &cam2);
			if(cam1<=0||cam1>cam_num||cam2<=0||cam2>cam_num)
			{
				printf("Cam name in Area.txt is wrong!\n");
				return 0;
			}
			cam1--;
			cam2--;
			char areaflag[100];
			fscanf(areafile,"%s",areaflag);
			if(!strcmp(areaflag,"Area:"))
			{
				n++;

				rl[cam1*cam_num+cam2].connect=1;
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area1.LeftUp.x, &rl[cam1*cam_num+cam2].area1.LeftUp.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area1.RightUp.x, &rl[cam1*cam_num+cam2].area1.RightUp.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area1.RightDown.x, &rl[cam1*cam_num+cam2].area1.RightDown.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area1.LeftDown.x, &rl[cam1*cam_num+cam2].area1.LeftDown.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area1.Center.x, &rl[cam1*cam_num+cam2].area1.Center.y);
				
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area2.LeftUp.x, &rl[cam1*cam_num+cam2].area2.LeftUp.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area2.RightUp.x, &rl[cam1*cam_num+cam2].area2.RightUp.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area2.RightDown.x, &rl[cam1*cam_num+cam2].area2.RightDown.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area2.LeftDown.x, &rl[cam1*cam_num+cam2].area2.LeftDown.y);
				fscanf(areafile,"%d %d", &rl[cam1*cam_num+cam2].area2.Center.x, &rl[cam1*cam_num+cam2].area2.Center.y);
				
				rl[cam2*cam_num+cam1].connect=1;
				rl[cam2*cam_num+cam1].area1.LeftUp=rl[cam1*cam_num+cam2].area2.LeftUp;
				rl[cam2*cam_num+cam1].area1.RightDown=rl[cam1*cam_num+cam2].area2.RightDown;
				rl[cam2*cam_num+cam1].area1.RightUp=rl[cam1*cam_num+cam2].area2.RightUp;
				rl[cam2*cam_num+cam1].area1.LeftDown=rl[cam1*cam_num+cam2].area2.LeftDown;
				rl[cam2*cam_num+cam1].area1.Center=rl[cam1*cam_num+cam2].area2.Center;

				rl[cam2*cam_num+cam1].area2.LeftUp=rl[cam1*cam_num+cam2].area1.LeftUp;
				rl[cam2*cam_num+cam1].area2.RightDown=rl[cam1*cam_num+cam2].area1.RightDown;
				rl[cam2*cam_num+cam1].area2.RightUp=rl[cam1*cam_num+cam2].area1.RightUp;
				rl[cam2*cam_num+cam1].area2.LeftDown=rl[cam1*cam_num+cam2].area1.LeftDown;
				rl[cam2*cam_num+cam1].area2.Center=rl[cam1*cam_num+cam2].area1.Center;
			}
		}
	}
	return 1;
}

float CW_MotionOperator::dist_motion(Tracklet_Info* tracklet1,Tracklet_Info* tracklet2,float* distance)
{
	int t,i,j;
	int cam_a=tracklet1->nCam;
	int cam_b=tracklet2->nCam;

	if (tracklet1->length<=step1)step1=tracklet1->length-1;
	if (tracklet2->length<=step2)step2=tracklet2->length-1;
	int inv_t=tracklet2->frm[0]-tracklet1->frm[tracklet1->length-1];
	float dist=0;
	float dist1=0;
	float dist2=0;
	if (cam_a==cam_b)
	{
		CvPoint2D32f p1=cvPoint2D32f(tracklet1->r[tracklet1->length-1-step1].x+tracklet1->r[tracklet1->length-1-step1].width/2.0,
			tracklet1->r[tracklet1->length-1-step1].y+tracklet1->r[tracklet1->length-1-step1].height/2.0);
		CvPoint2D32f p2=cvPoint2D32f(tracklet1->r[tracklet1->length-1].x+tracklet1->r[tracklet1->length-1].width/2.0,
			tracklet1->r[tracklet1->length-1].y+tracklet1->r[tracklet1->length-1].height/2.0);
		CvPoint2D32f p3=cvPoint2D32f(tracklet2->r[0].x+tracklet2->r[0].width/2.0,tracklet2->r[0].y+tracklet2->r[0].height/2.0);
		CvPoint2D32f p4=cvPoint2D32f(tracklet2->r[step2].x+tracklet2->r[step2].width/2.0,tracklet2->r[step2].y+tracklet2->r[step2].height/2.0);
		CvPoint2D32f p_1=cvPoint2D32f(p2.x+inv_t*(p2.x-p1.x)/(step1+INF_MIN),p2.y+inv_t*(p2.y-p1.y)/(step1+INF_MIN));
		CvPoint2D32f p_2=cvPoint2D32f(p3.x-inv_t*(p4.x-p3.x)/(step2+INF_MIN),p3.y-inv_t*(p4.y-p3.y)/(step2+INF_MIN));
		dist1 = cvSqrt((p_1.x-p3.x)*(p_1.x-p3.x)+(p_1.y-p3.y)*(p_1.y-p3.y));
		dist2 = cvSqrt((p_2.x-p2.x)*(p_2.x-p2.x)+(p_2.y-p2.y)*(p_2.y-p2.y));
	}
	else
	{
		CvPoint2D32f p1=cvPoint2D32f(tracklet1->r[tracklet1->length-1-step1].x+tracklet1->r[tracklet1->length-1-step1].width/2.0,
			tracklet1->r[tracklet1->length-1-step1].y+tracklet1->r[tracklet1->length-1-step1].height/2.0);
		CvPoint2D32f p2=cvPoint2D32f(tracklet1->r[tracklet1->length-1].x+tracklet1->r[tracklet1->length-1].width/2.0,
			tracklet1->r[tracklet1->length-1].y+tracklet1->r[tracklet1->length-1].height/2.0);
		CvPoint2D32f p3=cvPoint2D32f(tracklet2->r[0].x+tracklet2->r[0].width/2.0,tracklet2->r[0].y+tracklet2->r[0].height/2.0);
		CvPoint2D32f p4=cvPoint2D32f(tracklet2->r[step2].x+tracklet2->r[step2].width/2.0,tracklet2->r[step2].y+tracklet2->r[step2].height/2.0);

		CvMat* matdetect1 = (CvMat*)malloc(sizeof(CvMat*));
		CvMat* matdetect2 = (CvMat*)malloc(sizeof(CvMat*));
		CvPoint pt1,pt2;
		CvPoint pPoint1[4];
		CvPoint pPoint2[4];
		matdetect1=cvCreateMatHeader(1,4, CV_32SC2);
		matdetect2=cvCreateMatHeader(1,4, CV_32SC2);
		
		int c1=cam_a-1;
		int c2=cam_b-1;
		Assign1(rl[c1*cam_num+c2].area1.LeftUp,rl[c1*cam_num+c2].area1.RightUp,rl[c1*cam_num+c2].area1.RightDown,rl[c1*cam_num+c2].area1.LeftDown);
		pt1=rl[c1*cam_num+c2].area1.Center;
		Assign2(rl[c1*cam_num+c2].area2.LeftUp,rl[c1*cam_num+c2].area2.RightUp,rl[c1*cam_num+c2].area2.RightDown,rl[c1*cam_num+c2].area2.LeftDown);
		pt2=rl[c1*cam_num+c2].area2.Center;
		if(rl[c1*cam_num+c2].connect!=1)
			printf("No relationship between %d and %d\n",cam_a,cam_b);

		cvSetData(matdetect1,pPoint1,sizeof(CvPoint)*4);
		cvSetData(matdetect2,pPoint2,sizeof(CvPoint)*4);

		if (cvPointPolygonTest(matdetect1,p2,0)>=0)
		{
			dist1=0;
		}
		else
		{
			float distmin=INF;
			for (t=1;t<=inv_t;t++)
			{
				CvPoint2D32f p_1=cvPoint2D32f(p2.x+t*(p2.x-p1.x)/(step1+INF_MIN),p2.y+t*(p2.y-p1.y)/(step1+INF_MIN));
				float dist_tmp = cvSqrt((p_1.x-pt1.x)*(p_1.x-pt1.x)+(p_1.y-pt1.y)*(p_1.y-pt1.y));
				if (dist_tmp<distmin)distmin=dist_tmp;
			}
			dist1=distmin;
		}
		if (cvPointPolygonTest(matdetect2,p3,0)>=0)
		{
			dist2=0;
		}
		else
		{
			float distmin=INF;
			for (t=1;t<=inv_t;t++)
			{
				CvPoint2D32f p_2=cvPoint2D32f(p3.x-t*(p4.x-p3.x)/(step2+INF_MIN),p3.y-t*(p4.y-p3.y)/(step2+INF_MIN));
				float dist_tmp = cvSqrt((p_2.x-pt2.x)*(p_2.x-pt2.x)+(p_2.y-pt2.y)*(p_2.y-pt2.y));
				if (dist_tmp<distmin)distmin=dist_tmp;
			}
			dist2=distmin;
		}
	}
	dist=(dist1+dist2)/2;

	distance[0]=dist;

	float e=2.718281828459;
	float k=-0.01;
	float value=pow(e,k*dist);

	return value;
}

void CW_MotionOperator::release()
{
	free(rl);
}