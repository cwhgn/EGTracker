#include "EGTracker.h"
#include <Windows.h>
#include<tchar.h>

bool CW_EGTracker::LoadParam(const char* path)
{
	FILE *paramfile = fopen(path, "r");
	if (!paramfile) {
		fprintf(stderr, "Parammeters file (Param.ini) not available. Stopping.\n");
		exit(-1);
	}
	fclose(paramfile);

	LPTSTR lpPath = new char[MAX_PATH];	
	strcpy(lpPath, path);

	//CamInfo
	char strSection[100] = "CamInfo";
	char* strSectionKey = "";
	char* strValue = "";
	char inBuf[1024];

	strSectionKey = "camnum";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_camnum = atoi(strValue);

	//EGTracker
	sprintf(strSection, "EGTracker");

	strSectionKey = "maxflow";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_maxflow = atoi(strValue);

	strSectionKey = "waitframe";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_waitframe = atoi(strValue);

	strSectionKey = "mug_pairthrd";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_mugpairthrd = atof(strValue);

	//Gauss Background
	sprintf(strSection, "GaussBackground");

	strSectionKey = "bg_updatespeed";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_bgupdatespeed = atoi(strValue);

	strSectionKey = "bg_gaussthrd";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_bggaussthrd = atoi(strValue);

	strSectionKey = "fg_thrd";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_fgthrd = atoi(strValue);

	//AIF Tracker
	sprintf(strSection, "AIFTracker");

	strSectionKey = "conv_thrd";
	GetPrivateProfileString (strSection,strSectionKey, NULL, inBuf, 1024, lpPath); 
	strValue = inBuf;
	c_AIFconvthrd = atof(strValue);

	return 1;
}

bool CW_EGTracker::LoadData(const char* videopath, const char* gthpath, const char* topopath, const char* areapath)
{
	//Read groundtruth
	FILE *datasetfile = fopen(videopath, "r");
	if (!datasetfile) {
		fprintf(stderr, "Video list file (Video.txt) not available. Stopping.\n");
		exit(-1);
	}
	int camnum = ReadVideo(datasetfile);
	if(camnum!=c_camnum)
		printf("Video.txt contains different video number from camera number.\n");
	fclose(datasetfile);

	//Read video
	FILE *inputfile = fopen(gthpath, "r");
	if (!inputfile) {
		fprintf(stderr, "Input groundtruth file (Groundtruth.txt) not available. Stopping.\n");
		exit(-1);
	}
	c_gthnum = ReadGroundtruth(inputfile,0);
	fclose(inputfile);

	//Read topology information
	FILE *topofile = fopen(topopath, "r");
	if (!topofile) {
		fprintf(stderr, "Topology information file (Topology.txt) not available. Stopping.\n");
		exit(-1);
	}
	GetTopology(topofile);
	fclose(topofile);

	//Read exit/enter areas bewteen cameras
	c_motion.initial(c_camnum);
	FILE *areafile = fopen(areapath, "r");
	if (!areafile) {
		fprintf(stderr, "exit/enter area information file (Area.txt) not available. Stopping.\n");
		exit(-1);
	}
	GetArea(areafile);
	fclose(areafile);

	return 1;
}

void CW_EGTracker::Initial(void* param)
{
	int i;
	c_tracklet = (Tracklet_Info*)malloc(sizeof(Tracklet_Info)*MAX_TRACKLET);
	for (i=0;i<MAX_TRACKLET;i++)
	{
		c_tracklet[i].r = (CvRect*)malloc(sizeof(CvRect)*MAX_TRACKLET_LENGTH);
		c_tracklet[i].frm = (int*)malloc(sizeof(int)*MAX_TRACKLET_LENGTH);
		c_tracklet[i].nCam = -1;
		c_tracklet[i].id = -1;
		c_tracklet[i].length = -1;
		c_tracklet[i].conv = -1;
		c_tracklet[i].t = -1;
		c_tracklet[i].p = -1;
		c_tracklet[i].f_free = -1;
		//c_tracklet[i].f = (mcsh*)malloc(sizeof(mcsh)*MAX_TRACKLET_LENGTH);
	}

	c_camnum=0;
	c_gthnum=0;
	c_trkletnum=0;
	c_path = (char**)malloc(sizeof(char*)*MAX_CamNum);
	for (i = 0 ; i < MAX_CamNum ; i ++)
	{
		c_path[i] = (char*)malloc(sizeof(char)*100);
	}	
	c_groundtruth = (Boundingbox*)malloc(sizeof(Boundingbox)*MAX_GroundTruthNum);
	memset(c_groundtruth,0,MAX_GroundTruthNum*sizeof(Boundingbox));

	c_topo = (bool*)malloc(MAX_CamNum*MAX_CamNum*sizeof(bool));
	memset(c_topo,0,MAX_CamNum*MAX_CamNum*sizeof(bool));
	memset(c_sim,0,10*sizeof(float));
}

void CW_EGTracker::Process(void* param)
{
	int nodes,edges,maxf,waitframe;
	maxf = c_maxflow;
	waitframe = c_waitframe;

	//Obtain tracklets,PreProcess tracklet,Computing Sim
	GetTrklet();
	SelectTrklet(0);//the length of tracklets should be more than 0
	printf("Total tracklet number is£º%d\n",c_trkletnum);
	OrderTrklet();//Order tracklets by time

	//Output for test by step
	//char path[100]="F:\\Project\\UGMCT\\Results\\";
	//TrkletToOutput4Test(path,waitframe);

	//Compute the mean and variance
	Compute_Sims(waitframe);
	printf("SCT:mean=%.4f,var=%.4f\nICT:mean=%.4f,var=%.4f\n",c_sim[1],c_sim[3],c_sim[2],c_sim[4]);

	//compute the result
	char dmx_path[100]="sample.dmx";
	OutputToDMX(dmx_path,&nodes,&edges,maxf,waitframe);
	c_cs2.Initial(dmx_path,nodes,edges,maxf);
	c_cs2.LoadGraph();
	c_cs2.Process();
	int connects=0;
	int* links=(int*)malloc(nodes*sizeof(int));
	memset(links,0,nodes*sizeof(int));
	c_cs2.GetResult(links,&connects);

	UpdateTrklet(links);
}

void CW_EGTracker::GetResult(Boundingbox* output, int* length)
{		
	int i,j;
	int l=0;
	for (i=0;i<c_trkletnum;i++)
	{
		for (j=0;j<c_tracklet[i].length;j++)
		{
			output[l].nCam=c_tracklet[i].nCam;
			output[l].nFrm=c_tracklet[i].frm[j];
			output[l].ObjID=c_tracklet[i].id;
			output[l].r=c_tracklet[i].r[j];
			l++;
		}
	}
	length[0]=l;
}

void CW_EGTracker::Release()
{
	int i;
	free(c_groundtruth);
	for (i=0;i<MAX_TRACKLET;i++)
	{
		free(c_tracklet[i].r);
		free(c_tracklet[i].frm);
		//free(tracklet[i].f);
		if(c_tracklet[i].t>0)free(c_tracklet[i].pf);
	}
	free(c_tracklet);
	for (i=0;i<MAX_CamNum;i++)
	{
		free(c_path[i]);
	}	
	free(c_path);
	free(c_topo);
	c_motion.release();
}

////////////////////////////////////////////PRIVATE////////////////////////////////////////////////////////

int CW_EGTracker::ReadVideo(FILE * datasetfile)
{
	int i;
	int cams=0;
	for(i=0;feof(datasetfile)==0;i++)
	{
		int flag=fscanf(datasetfile,"%s",c_path[i]);
		if(flag>0)cams++;
	}
	return cams;
}

int CW_EGTracker::ReadGroundtruth(FILE * inputfile,int cam)//×¢:cam=0 indicates all cameras,cam=n means data of camera n.
{
	int i;
	int length=0;
	for(i=0;feof(inputfile)==0;i++)
	{
		int flag=fscanf(inputfile,"%d,%d,%d,%d,%d,%d,%d",&c_groundtruth[i].nCam,&c_groundtruth[i].nFrm,&c_groundtruth[i].ObjID,
			&c_groundtruth[i].r.x,&c_groundtruth[i].r.y,&c_groundtruth[i].r.width,&c_groundtruth[i].r.height);
		c_groundtruth[i].ObjID=-1;//for exoeriment two£¬all ID should be -1
		if(flag>0)
		{
			length++;
			if(cam>0&&c_groundtruth[i].nCam!=cam)
			{
				i--;
				length--;
			}
		}
	}
	return length;
}

void CW_EGTracker::pyrfilter(IplImage* src,IplImage* dst)
{   
	CvSize size=cvSize(src->width,src->height);
	IplImage* pyr=cvCreateImage(cvSize((size.width&(-2))/2,(size.height&(-2))/2),8,1);
	cvPyrDown(src,pyr,7);
	cvDilate(pyr,pyr,0,1);
	cvErode(pyr,pyr,NULL,1);
	cvPyrUp(pyr,dst,7);
	cvReleaseImage(&pyr);
}

void CW_EGTracker::GetFG(IplImage* src1, IplImage* src2, IplImage* dst)
{
	//CvRect rect = cvRect(0,0,dst->width,1);
	cvAbsDiff(src1,src2,dst);
	pyrfilter(dst,dst);
	//cvSetImageROI(dst,rect);
	//cvZero(dst);
	//cvResetImageROI(dst);
	cvSmooth(dst,dst);
	cvThreshold(dst,dst,c_fgthrd,255,CV_THRESH_BINARY);
}

void CW_EGTracker::GetTrklet()
{
	int i,j;
	IplImage* frame;
	int sum_gth=0;
	for (int cam=1;cam<=c_camnum;cam++)
	{
		CvCapture* pCapture = NULL;
		pCapture = cvCaptureFromFile(c_path[cam-1]);
		frame = cvQueryFrame( pCapture );
		cvFlip(frame,frame);frame->origin=0;

		//Gauss
		SVBASIC *staticbg;
		InParamSVFinalBG ipbg;
		InDataSVFinalBG idbg;
		OutPutSVFinalBG opbg_statics;
		//Gauss initial
		ipbg._imgw = frame->width;ipbg._imgh = frame->height;
		SVCreateObject_SVStaticsBG(&staticbg);
		ipbg._thresh=c_bggaussthrd;
		ipbg._iniframe=c_bgupdatespeed;
		staticbg->Initial(&ipbg);
		//others
		int nFramNum=0;
		IplImage* gray = cvCreateImage(cvGetSize(frame),frame->depth,1);
		IplImage* bgimage = cvCreateImage(cvGetSize(frame),frame->depth,1);
		IplImage* fgimage = cvCreateImage(cvGetSize(frame),frame->depth,1);
		cvCvtColor(frame,gray,CV_RGB2GRAY);
		cvZero(bgimage);

		//IIFTracker initial
		int tracker_idx[MAX_TRACKER];
		APPSVInterface* ppBaseModule[MAX_TRACKER];
		for (i=0;i<MAX_TRACKER;i++)
		{
			tracker_idx[i]=-1;
			APPSVCreateObject_IIFTrack(&ppBaseModule[i], ALG_ALARM_IIFTRACK);	
		}

		bool do_while=1;
		while(do_while)
		{
			//Create the foreground image		
			idbg._gray = gray;
			staticbg->Process(&idbg);
			staticbg->GetResult(&opbg_statics);
			if (opbg_statics._bgok)
			{
				cvCopyImage(opbg_statics._bg,bgimage);
			}
			GetFG(bgimage,gray,fgimage);

			int cur_output_num=0;
			Boundingbox cur_output[MAX_TRACKER];
			for (i=0;i<c_gthnum;i++)
			{
				if (c_groundtruth[i].nCam == cam && c_groundtruth[i].nFrm==nFramNum)
				{
					memcpy(&cur_output[cur_output_num],&c_groundtruth[i],sizeof(Boundingbox));
					cur_output_num++;
				}
			}

			GetTrkletbyFrm(frame,fgimage,ppBaseModule,tracker_idx,cur_output,cur_output_num,nFramNum,cam);

			if (frame = cvQueryFrame( pCapture ))
			{
				cvFlip(frame,frame);frame->origin=0;
				cvCvtColor(frame,gray,CV_RGB2GRAY);
				nFramNum++;
			}
			else do_while=0;
		}

		staticbg->Release();
		cvReleaseImage(&gray);
		cvReleaseImage(&bgimage);
		cvReleaseImage(&fgimage);
		cvReleaseCapture(&pCapture);
		printf("cam[%d]finished\n",cam);
	}

	for (i=0;i<c_trkletnum;i++)
	{
		if (c_tracklet[i].f!=NULL)
		{
			c_pmcshr.cal_tl_pmcsh(&c_tracklet[i]);
			c_pmcshr.cal_tl_mcsh(&c_tracklet[i]);
			free(c_tracklet[i].f);
			c_tracklet[i].f=NULL;
			c_tracklet[i].f_free=9;

			printf("free mcsh: ID-%d,idx-%d,lastframe-%d\n",c_tracklet[i].id,i,c_tracklet[i].frm[0]+c_tracklet[i].length);
		}
	}

}

void CW_EGTracker::GetTrkletbyFrm(IplImage* frame,IplImage* fgimage,APPSVInterface** ppBaseModule,int* tracker_idx,Boundingbox* cur_output,int cur_output_num, int nFrm, int Cam)
{
	int i,j;
	OutPutAPP_IIFTracker outparam[MAX_TRACKER];
	bool* cur_state=(bool*)malloc(cur_output_num*sizeof(bool));
	memset(cur_state,0,cur_output_num*sizeof(bool));
	bool* trk_state=(bool*)malloc(MAX_TRACKER*sizeof(bool));
	memset(trk_state,0,MAX_TRACKER*sizeof(bool));
	float* overlap=(float*)malloc(MAX_TRACKER*cur_output_num*sizeof(float));
	memset(overlap,0,MAX_TRACKER*cur_output_num*sizeof(float));

	//Build the relationship between current_rects and tracker_results
	int trackers=0;
	for (i=0;i<MAX_TRACKER;i++)
	{
		if (tracker_idx[i]>=0)
		{
			trackers++;
			ppBaseModule[i]->Process((unsigned char *)frame->imageData);
			ppBaseModule[i]->GetResult(&outparam[i]);
			if (outparam[i].conv<c_AIFconvthrd)continue;
			for (j=0;j<cur_output_num;j++)
			{
				overlap[i*cur_output_num+j]=TwoRectOverlap(cur_output[j].r,outparam[i].r,2);
			}
		}
	}
	//printf("trackers=%d,currentRs=%d\n",trackers,cur_output_num);

	//do while
	int max_trk=-1,max_cur=-1;
	while (1)
	{
		float max_overlap=0;
		for (i=0;i<MAX_TRACKER;i++)
		{
			if (tracker_idx[i]>=0)
			{
				for (j=0;j<cur_output_num;j++)
				{
					if (cur_state[j]!=1&&trk_state[i]!=1&&max_overlap<overlap[i*cur_output_num+j])
					{
						max_trk=i;
						max_cur=j;
						max_overlap=overlap[i*cur_output_num+j];
					}
				}
			}
		}
		if(max_overlap<=0)break;

		cur_state[max_cur]=1;
		trk_state[max_trk]=1;
		overlap[max_trk*cur_output_num+max_cur]=0;
		int ID=tracker_idx[max_trk];
		CvRect r=cur_output[max_cur].r;
		ppBaseModule[max_trk]->Release();tracker_idx[max_trk]=-1;
		NewTracker(frame,ppBaseModule,tracker_idx,r,ID,max_trk);//Rect=outparam[max_trk].r,ID=tracker_idx[max_trk];
		AddtoTrklet(frame,fgimage,r,ID,nFrm,Cam,outparam[max_trk].conv);//ID=tracker_idx[i],Rect=r;

		//draw
		CvRect tmp_r=r;
		cvRectangle(frame,cvPoint(tmp_r.x,tmp_r.y),cvPoint(tmp_r.x+tmp_r.width,tmp_r.y+tmp_r.height),CV_RGB(0,0,255));
		CvFont font;
		char text[100];
		sprintf(text,"%d",ID);
		cvInitFont( &font,CV_FONT_VECTOR0,0.5,0.5,0,1,8);
		cvPutText(frame,text,cvPoint(tmp_r.x+tmp_r.width/2,tmp_r.y+tmp_r.height/2), &font,CV_RGB(0,0,255));
	}

	//free the left tracker
	for (i=0;i<MAX_TRACKER;i++)
	{
		if (trk_state[i]!=1&&tracker_idx[i]>=0)
		{
			EndTrklet(tracker_idx[i]);
			ppBaseModule[i]->Release();tracker_idx[i]=-1;
		}
	}

	//new the left cur_output
	for (i=0;i<cur_output_num;i++)
	{
		if (cur_state[i]==0)
		{
			NewTracker(frame,ppBaseModule,tracker_idx,cur_output[i].r,c_trkletnum,-1);//Rect=cur_output[i].r,ID=c_trkletnum;
			NewTrklet(frame,fgimage,cur_output[i].r,c_trkletnum,nFrm,Cam);//ID=c_trkletnum,c_trkletnum++,Rect=cur_output[i].r,Frm=nFrm,Cam=Cam;

			//draw
			CvRect tmp_r=cur_output[i].r;
			cvRectangle(frame,cvPoint(tmp_r.x,tmp_r.y),cvPoint(tmp_r.x+tmp_r.width,tmp_r.y+tmp_r.height),CV_RGB(255,0,0));
			CvFont font;
			char text[100];
			sprintf(text,"%d",c_trkletnum-1);
			cvInitFont( &font,CV_FONT_VECTOR0,0.5,0.5,0,1,8);
			cvPutText(frame,text,cvPoint(tmp_r.x+tmp_r.width/2,tmp_r.y+tmp_r.height/2), &font,CV_RGB(255,0,0));
		}
	}

	cvNamedWindow("frame");cvShowImage("frame",frame);
	if (cvWaitKey(10)=='s')cvWaitKey(0);

	free(cur_state);
	free(trk_state);
	free(overlap);
}

double CW_EGTracker::TwoRectOverlap(CvRect r1,CvRect r2,int flag)//The denominator is: smaller rect(flag=0)£¬larger rect (flag=1)£¬r1 (flag=2)£¬r2 (flag=3)
{
	int w1,w2,h1,h2,w,h,wmin,hmin,area,coarea;
	double overlapK=0;

	w1 = abs(r2.x+r2.width-r1.x);
	w2 = abs(r1.x+r1.width-r2.x);
	h1 = abs(r2.y+r2.height-r1.y);
	h2 = abs(r1.y+r1.height-r2.y);
	w = r1.width+r2.width;
	h = r1.height+r2.height;

	if(w1<w&&w2<w&&h1<h&&h2<h)
	{
		wmin = min(min(w1,w2),min(r1.width,r2.width));
		hmin = min(min(h1,h2),min(r1.height,r2.height));
		area = wmin*hmin;

		if(flag==0)coarea = min(r1.height*r1.width,r2.height*r2.width);
		else if(flag==1)coarea = max(r1.height*r1.width,r2.height*r2.width);
		else if(flag==2)coarea = r1.height*r1.width;
		else if(flag==3)coarea = r2.height*r2.width;

		overlapK = ((double)area/coarea);
		//printf("overlapK=%f\n",overlapK);
	}
	return overlapK;
}

void CW_EGTracker::AddtoTrklet(IplImage* frame,IplImage* fgimage, CvRect r, int ID, int nFrm, int Cam, float conv)
{
	int i;
	for (i=0;i<c_trkletnum;i++)
	{
		if (c_tracklet[i].id==ID)
		{
			int l=c_tracklet[i].length;
			c_tracklet[i].frm[l] = nFrm;
			c_tracklet[i].r[l] = r;
			c_tracklet[i].conv = (c_tracklet[i].conv*c_tracklet[i].length + conv)/(c_tracklet[i].length+1);
			c_tracklet[i].length++;
			c_pmcshr.cal_img_mcsh(frame,fgimage,r,&c_tracklet[i].f[l]);
			return;
		}
	}
}

void CW_EGTracker::NewTrklet(IplImage* frame,IplImage* fgimage, CvRect r, int ID, int nFrm, int Cam)
{
	c_tracklet[c_trkletnum].length = 0;
	c_tracklet[c_trkletnum].frm[0] = nFrm;
	c_tracklet[c_trkletnum].r[0] = r;
	c_tracklet[c_trkletnum].nCam = Cam;
	c_tracklet[c_trkletnum].id = ID;
	c_tracklet[c_trkletnum].length++;
	c_tracklet[c_trkletnum].f = (mcsh*)malloc(sizeof(mcsh)*MAX_TRACKLET_LENGTH);
	c_tracklet[c_trkletnum].conv = 1;
	c_pmcshr.cal_img_mcsh(frame,fgimage,r,&c_tracklet[c_trkletnum].f[0]);
	c_trkletnum++;
}

void CW_EGTracker::NewTracker(IplImage* frame,APPSVInterface** ppBaseModule,int* tracker_idx,CvRect r,int ID,int new_idx)
{
	int i;
	InPutAPP_IIFTracker inparam;

	if(new_idx>0)
	{
		if(tracker_idx[new_idx]!=-1)
			printf("NewTracker BUG!\n");
		i=new_idx;
	}
	else
	{
		for (i=0;i<MAX_TRACKER;i++)
		{
			if (tracker_idx[i]<0)
			{
				break;
			}
		}
	}

	//APPSVCreateObject_IIFTrack(&ppBaseModule[i], ALG_ALARM_IIFTRACK);						
	inparam.img=frame;
	inparam.nNum=60;
	inparam.r=r;
	ppBaseModule[i]->Initial(&inparam);
	tracker_idx[i]=ID;
}

void CW_EGTracker::EndTrklet(int ID)
{
	int i;
	for (i=0;i<c_trkletnum;i++)
	{
		if (c_tracklet[i].id==ID)
		{
			c_pmcshr.cal_tl_pmcsh(&c_tracklet[i]);
			c_pmcshr.cal_tl_mcsh(&c_tracklet[i]);
			free(c_tracklet[i].f);
			c_tracklet[i].f=NULL;
			c_tracklet[i].f_free=9;
			return;
		}
	}
}

void CW_EGTracker::SelectTrklet(int length_thrd)
{
	int i;
	for (i=0;i<c_trkletnum;i++)
	{
		if (c_tracklet[i].length<=length_thrd)
		{
			if (c_trkletnum-1==i)
			{
				DelTrklet(&c_tracklet[i]);
			}
			else
			{
				CpyTrklet(&c_tracklet[i],&c_tracklet[c_trkletnum-1]);
			}
			i--;
			c_trkletnum--;
		}
	}
}

void CW_EGTracker::CpyTrklet(Tracklet_Info* tracklet1,Tracklet_Info* tracklet2)
{
	memcpy(tracklet1->r,tracklet2->r,MAX_TRACKLET_LENGTH*sizeof(CvRect));
	memset(tracklet2->r,0,MAX_TRACKLET_LENGTH*sizeof(CvRect));
	memcpy(tracklet1->frm,tracklet2->frm,MAX_TRACKLET_LENGTH*sizeof(int));
	memset(tracklet2->frm,0,MAX_TRACKLET_LENGTH*sizeof(int));
	if (tracklet1->pf!=0)
	{
		free(tracklet1->pf);tracklet1->pf=0;
	}
	if (tracklet2->pf!=0)
	{
		tracklet1->pf=(mcsh*)malloc(sizeof(mcsh)*tracklet2->p);
		memcpy(tracklet1->pf,tracklet2->pf,tracklet2->p*sizeof(mcsh));
		free(tracklet2->pf);tracklet2->pf=0;
	}
	memcpy(&tracklet1->m,&tracklet2->m,sizeof(mcsh));
	tracklet1->id=tracklet2->id;
	tracklet2->id=-1;
	tracklet1->length=tracklet2->length;
	tracklet2->length=-1;
	tracklet1->nCam=tracklet2->nCam;
	tracklet2->nCam=-1;
	tracklet1->p=tracklet2->p;
	tracklet2->p=-1;
	tracklet1->t=tracklet2->t;
	tracklet2->t=-1;
	tracklet1->conv=tracklet2->conv;
	tracklet2->conv=-1;
	tracklet1->f_free=tracklet2->f_free;
	tracklet2->f_free=-1;
}

void CW_EGTracker::DelTrklet(Tracklet_Info* tracklet2)
{
	memset(tracklet2->r,0,MAX_TRACKLET_LENGTH*sizeof(CvRect));
	memset(tracklet2->frm,0,MAX_TRACKLET_LENGTH*sizeof(int));
	memset(&tracklet2->m,0,sizeof(mcsh));
	if (tracklet2->pf!=0)
	{
		free(tracklet2->pf);tracklet2->pf=0;
	}
	tracklet2->id=-1;
	tracklet2->length=-1;
	tracklet2->nCam=-1;
	tracklet2->p=-1;
	tracklet2->t=-1;
	tracklet2->conv=-1;
	tracklet2->f_free=-1;
}

void CW_EGTracker::OrderTrklet()
{
	int i,j;
	for (i=0;i<c_trkletnum-1;i++)
	{
		for (j=0;j<c_trkletnum-1-i;j++)
		{
			if (c_tracklet[j].frm[0]>c_tracklet[j+1].frm[0])
			{
				SwitchTrklet(&c_tracklet[j],&c_tracklet[j+1]);
			}
		}
	}
}

void CW_EGTracker::SwitchTrklet(Tracklet_Info* tracklet1,Tracklet_Info* tracklet2)
{
	Tracklet_Info tmptrk;
	tmptrk.r = (CvRect*)malloc(sizeof(CvRect)*MAX_TRACKLET_LENGTH);
	tmptrk.frm = (int*)malloc(sizeof(int)*MAX_TRACKLET_LENGTH);
	tmptrk.pf = (mcsh*)malloc(sizeof(mcsh)*9);
	CpyTrklet(&tmptrk,tracklet1);
	CpyTrklet(tracklet1,tracklet2);
	CpyTrklet(tracklet2,&tmptrk);
	free(tmptrk.r);
	free(tmptrk.frm);
}

void CW_EGTracker::GetTopology(FILE* topofile)
{
	int i,a,b;
	int cam_num=c_camnum;
	for(a=0;a<cam_num;a++)
	{
		for(b=0;b<cam_num;b++)
		{
			int flag=fscanf(topofile,"%d",&c_topo[a*cam_num+b]);
			printf("a=%d,b=%d,topo=%d\n",a,b,c_topo[a*cam_num+b]);
			if(flag<0)
				printf("Topology.txt contains not enough information!\n");
		}
	}
}

void CW_EGTracker::GetArea(FILE* areafile)
{
	c_motion.load_areas(areafile);
}

void CW_EGTracker::Compute_Sims(int waitframe)
{
	int i,j;
	float sum_conv=0;
	float sum_sim_s=0;
	float sum_sim_c=0;
	int sim_pair_s=0;
	int sim_pair_c=0;
	for (i=0;i<c_trkletnum;i++)
	{
		sum_conv += c_tracklet[i].conv;
		for (j=0;j<c_trkletnum;j++)
		{
			int inv_t=c_tracklet[j].frm[0]-c_tracklet[i].frm[c_tracklet[i].length-1];
			int a=c_tracklet[i].nCam;
			int b=c_tracklet[j].nCam;
			float flag_mug=c_pmcshr.dist_mug(&c_tracklet[i],&c_tracklet[j]);
			if(inv_t>0&&inv_t<waitframe&&c_topo[(a-1)*c_camnum+(b-1)]==1&&flag_mug<=c_mugpairthrd)
			{
				float tmp_sim=c_pmcshr.dist_avg_tl(&c_tracklet[i],&c_tracklet[j]);
				//appearance
				if (a==b) {sum_sim_s += tmp_sim; sim_pair_s++;}
				else      {sum_sim_c += tmp_sim; sim_pair_c++;}
			}
		}
	}
	c_sim[0]=sum_conv/c_trkletnum;
	float sim_s=sum_sim_s/sim_pair_s;
	c_sim[1] = sim_s;
	float sim_c=sum_sim_c/sim_pair_c;
	c_sim[2] = sim_c;
	//printf("appearance pair:s-%d,c-%d; motion pair:s-%d,c-%d\n",sim_pair_s,sim_pair_c,motion_pair_s,motion_pair_c);

	float sum_var_s=0;
	float sum_var_c=0;
	for (i=0;i<c_trkletnum;i++)
	{
		for (j=0;j<c_trkletnum;j++)
		{
			int inv_t=c_tracklet[j].frm[0]-c_tracklet[i].frm[c_tracklet[i].length-1];
			int a=c_tracklet[i].nCam;
			int b=c_tracklet[j].nCam;
			float flag_mug=c_pmcshr.dist_mug(&c_tracklet[i],&c_tracklet[j]);
			if(inv_t>0&&inv_t<waitframe&&c_topo[(a-1)*c_camnum+(b-1)]==1&&flag_mug<=c_mugpairthrd)
			{
				float tmp_sim=c_pmcshr.dist_avg_tl(&c_tracklet[i],&c_tracklet[j]);
				if (a==b) sum_var_s += (tmp_sim-sim_s)*(tmp_sim-sim_s);
				else      sum_var_c += (tmp_sim-sim_c)*(tmp_sim-sim_c);
			}
		}
	}
	float var_s=sum_var_s/sim_pair_s;
	c_sim[3] = cvSqrt(var_s);
	float var_c=sum_var_c/sim_pair_c;
	c_sim[4] = cvSqrt(var_c);
}

void CW_EGTracker::TrkletToOutput4Test(char* path,int waitframe)
{
	int i,j;

	//Output Tracklets
	char trklet_path[100];
	strcat(strcpy(trklet_path,path),"tracklet.txt");
	FILE *trklet_opf = fopen(trklet_path, "w");
	for (i=0;i<c_trkletnum;i++)
	{
		for (j=0;j<c_tracklet[i].length;j++)
		{
			fprintf(trklet_opf,"%d,%d,%d,%d,%d,%d,%d,%d\n",c_tracklet[i].nCam,c_tracklet[i].frm[j],c_tracklet[i].id,i,
				c_tracklet[i].r[j].x,c_tracklet[i].r[j].y,c_tracklet[i].r[j].width,c_tracklet[i].r[j].height);
		}
	}
	fclose(trklet_opf);

	//Output Other Features
	int m_length=c_trkletnum*c_trkletnum;
	char conv_path[100];
	char a_path[100];
	char non_a_path[100];
	char m_path[100];
	char t_path[100];
	char data_path[100];
	char scene_path[100];
	char mug_path[100];
	strcat(strcpy(conv_path,path),"conv_m.txt");
	strcat(strcpy(a_path,path),"a_m.txt");
	strcat(strcpy(non_a_path,path),"non_a_m.txt");
	strcat(strcpy(m_path,path),"m_m.txt");
	strcat(strcpy(t_path,path),"t_m.txt");
	strcat(strcpy(data_path,path),"data_m.txt");
	strcat(strcpy(scene_path,path),"scene_m.txt");
	strcat(strcpy(mug_path,path),"mug_m.txt");
	float* a_m=(float*)malloc(m_length*sizeof(float));
	float* non_a_m=(float*)malloc(m_length*sizeof(float));
	float* m_m=(float*)malloc(m_length*sizeof(float));
	float* t_m=(float*)malloc(m_length*sizeof(float));
	int* data_m=(int*)malloc(m_length*sizeof(int));
	int* scene_m=(int*)malloc(m_length*sizeof(int));
	float* mug_m=(float*)malloc(m_length*sizeof(float));
	memset(a_m,0,m_length*sizeof(float));
	memset(non_a_m,0,m_length*sizeof(float));
	memset(m_m,0,m_length*sizeof(float));
	memset(t_m,0,m_length*sizeof(float));
	memset(scene_m,0,m_length*sizeof(int));
	memset(data_m,0,m_length*sizeof(int));
	memset(mug_m,0,m_length*sizeof(float));

	int near_thrd=5;
	for (i=1;i<=c_trkletnum;i++)
	{
		for (j=1;j<=c_trkletnum;j++)
		{
			int inv_t=c_tracklet[j-1].frm[0]-c_tracklet[i-1].frm[c_tracklet[i-1].length-1];
			int a=c_tracklet[i-1].nCam;
			int b=c_tracklet[j-1].nCam;
			if (inv_t>0&&c_topo[(a-1)*c_camnum+(b-1)]==1)
			{
				mug_m[(i-1)*c_trkletnum + j-1] = c_pmcshr.dist_mug(&c_tracklet[i-1],&c_tracklet[j-1]);

				float mean = 0;
				float var = 1;
				if(a==b)
				{
					mean = c_sim[1]-c_sim[2];
					var = c_sim[3]/c_sim[4];
				}

				non_a_m[(i-1)*c_trkletnum + j-1] = c_pmcshr.dist_avg_tl(&c_tracklet[i-1],&c_tracklet[j-1]);
				a_m[(i-1)*c_trkletnum + j-1] = ( non_a_m[(i-1)*c_trkletnum + j-1] - mean )/var;
				float dist=0;
				m_m[(i-1)*c_trkletnum + j-1] = c_motion.dist_motion(&c_tracklet[i-1],&c_tracklet[j-1],&dist);
				t_m[(i-1)*c_trkletnum + j-1]= (inv_t<=near_thrd)?1.0:1-min(1,(inv_t-near_thrd)/(float)(waitframe-near_thrd+1));

				data_m[(i-1)*c_trkletnum + j-1]=1;
				scene_m[(i-1)*c_trkletnum + j-1]=a*10+b;
			}
		}
	}

	FILE* conv_f=fopen(conv_path,"w");
	FILE* a_f=fopen(a_path,"w");
	FILE* non_a_f=fopen(non_a_path,"w");
	FILE* m_f=fopen(m_path,"w");
	FILE* t_f=fopen(t_path,"w");
	FILE* data_f=fopen(data_path,"w");
	FILE* scene_f=fopen(scene_path,"w");
	FILE* mug_f=fopen(mug_path,"w");
	for(i=0;i<c_trkletnum;i++)
	{
		fprintf(conv_f,"%.4f ",c_tracklet[i].conv);
		for(j=0;j<c_trkletnum;j++)
		{
			fprintf(a_f,"%.4f ",a_m[i*c_trkletnum+j]);
			fprintf(non_a_f,"%.4f ",non_a_m[i*c_trkletnum+j]);
			fprintf(m_f,"%.4f ",m_m[i*c_trkletnum+j]);
			fprintf(t_f,"%.4f ",t_m[i*c_trkletnum+j]);
			fprintf(data_f,"%d ",data_m[i*c_trkletnum+j]);
			fprintf(scene_f,"%d ",scene_m[i*c_trkletnum+j]);
			fprintf(mug_f,"%f ",mug_m[i*c_trkletnum+j]);
		}fprintf(a_f,"\n");fprintf(non_a_f,"\n");fprintf(m_f,"\n");fprintf(t_f,"\n");fprintf(scene_f,"\n");fprintf(conv_f,"\n");fprintf(data_f,"\n");fprintf(mug_f,"\n");
	}
	fclose(conv_f);
	fclose(a_f);
	fclose(non_a_f);
	fclose(m_f);
	fclose(t_f);
	fclose(data_f);
	fclose(scene_f);
	fclose(mug_f);

	free(a_m);
	free(non_a_m);
	free(m_m);
	free(t_m);
	free(scene_m);
	free(data_m);
	free(mug_m);
}

void CW_EGTracker::OutputToDMX(char* dmx_path,int* n, int* e, int maxf,int waitframe)
{
	int i,j;
	FILE *dmx_f = fopen(dmx_path, "w");

	int edges=0;
	for (i=1;i<=c_trkletnum;i++)
	{
		edges=edges+3;
		for (j=1;j<=c_trkletnum;j++)
		{
			int inv_t=c_tracklet[j-1].frm[0]-c_tracklet[i-1].frm[c_tracklet[i-1].length-1];
			int a=c_tracklet[i-1].nCam;
			int b=c_tracklet[j-1].nCam;
			if (inv_t>0&&(inv_t<waitframe)&&c_topo[(a-1)*c_camnum+(b-1)]==1)edges++;
		}
	}
	int nodes=c_trkletnum*2+2;
	e[0]=edges;
	n[0]=nodes;

	fprintf(dmx_f,"p min %d %d\n", nodes, edges);
	fprintf(dmx_f,"n %d %d\n", 1, maxf);
	fprintf(dmx_f,"n %d %d\n", nodes, 0-maxf);

	float t_sim=0;
	float a_sim=0;
	float non_a_sim=0;
	float m_sim=0;
	float pexit=0;// or using c_maxflow/trklet_num;
	int near_thrd=5;
	for (i=1;i<=c_trkletnum;i++)
	{
		fprintf(dmx_f,"a %d %d %d %d %f\n",1,i*2,0,1,pexit);
		fprintf(dmx_f,"a %d %d %d %d %f\n",i*2+1,nodes,0,1,pexit);
		float tmp=c_tracklet[i-1].conv<=0?INF_MIN:c_tracklet[i-1].conv;
		tmp=c_tracklet[i-1].conv>=1?1-INF_MIN:c_tracklet[i-1].conv;
		fprintf(dmx_f,"a %d %d %d %d %f\n",i*2,i*2+1,0,1,1.0*log10((1-tmp)/tmp));
		for (j=1;j<=c_trkletnum;j++)
		{
			int inv_t=c_tracklet[j-1].frm[0]-c_tracklet[i-1].frm[c_tracklet[i-1].length-1];
			int a=c_tracklet[i-1].nCam;
			int b=c_tracklet[j-1].nCam;
			if (inv_t>0&&(inv_t<waitframe)&&c_topo[(a-1)*c_camnum+(b-1)]==1)
			{
				float mean = 0;
				float var = 1;
				if(a==b)
				{
					mean = c_sim[1]-c_sim[2];
					var = c_sim[3]/c_sim[4];
				}

				non_a_sim = c_pmcshr.dist_avg_tl(&c_tracklet[i-1],&c_tracklet[j-1]);
				a_sim = ( non_a_sim - mean )/var;
				float dist=0;
				m_sim = c_motion.dist_motion(&c_tracklet[i-1],&c_tracklet[j-1],&dist);
				t_sim= (inv_t<=near_thrd)?1.0:1-min(1,(inv_t-near_thrd)/(float)(waitframe-near_thrd+1));

				float k1=0.5;
				float k2=0.5;
				float k3=0.0;

				non_a_sim = non_a_sim<=0?INF_MIN:non_a_sim;non_a_sim = non_a_sim>=1?1:non_a_sim;
				a_sim = a_sim<=0?INF_MIN:a_sim;a_sim = a_sim>=1?1:a_sim;
				m_sim = m_sim<=0?INF_MIN:m_sim;m_sim = m_sim>=1?1:m_sim;
				t_sim = t_sim<=0?INF_MIN:t_sim;t_sim = t_sim>=1?1:t_sim;

				int miusN=6;
				float final_sim=-k1*log10(a_sim)-k2*log10(m_sim)-k3*log10(t_sim);

				fprintf(dmx_f,"a %d %d %d %d %f\n",i*2+1,j*2,0,1,final_sim);
			}
		}
	}
	fclose(dmx_f);
}

void CW_EGTracker::UpdateTrklet(int* links)
{
	int i,j;
	int* best_result=(int*)malloc(c_trkletnum*sizeof(int));
	for (i=0;i<c_trkletnum;i++)best_result[i]=-1;
	int connects=0;

	for (i=1;i<=c_trkletnum;i++)
	{
		if (links[i*2]!=0)
		{
			int Endn=(links[i*2]+1)/2;
			if((links[i*2]+1)%2!=0)
				printf("A wrong link appears! BUG!\n");
			if(best_result[Endn-1]!=-1&&Endn!=c_trkletnum+1)
				printf("Too many edges to one node! BUG!\n");
			best_result[Endn-1]=i-1;
			connects++;
		}
	}
	printf("connects:%d\n",connects);

	int new_id_num=0;
	int* new_id=(int*)malloc(c_trkletnum*c_trkletnum*sizeof(int));
	int* new_id_idx=(int*)malloc(c_trkletnum*sizeof(int));
	for (i=0;i<c_trkletnum;i++)
	{
		if (best_result[i]==-1)
		{
			new_id[new_id_num*c_trkletnum]=i;
			new_id_idx[new_id_num]=1;
			c_tracklet[i].id=new_id_num;
			new_id_num++;
		}
	}
	for (i=0;i<c_trkletnum;i++)
	{
		if (best_result[i]!=-1)
		{
			int k=FindPreTracklet(i,best_result);
			for (j=0;j<new_id_num;j++)
			{
				if (new_id[j*c_trkletnum]==k)
				{
					new_id[j*c_trkletnum+new_id_idx[j]]=i;
					new_id_idx[j]++;
					c_tracklet[i].id=j;
				}
			}
		}
	}

}

int CW_EGTracker::FindPreTracklet(int i,int* result)
{
	int k=i;
	int pre_i=result[i];
	if (pre_i>=0)
	{		
		k=FindPreTracklet(pre_i,result);
	}
	return k;
}