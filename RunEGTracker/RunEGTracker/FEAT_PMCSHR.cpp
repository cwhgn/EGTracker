#include "FEAT_PMCSHR.h"

void CW_PMCSHROperator::cal_img_mcsh(IplImage* frame,IplImage* frame_bk, CvRect rect, mcsh * mcsh_region)
{
	float * mr;

	IplImage * tmp , * tmp_bk;

	cvSetImageROI( frame, rect);
	tmp = cvCreateImage( cvGetSize(frame), 8, 3 );
	cvCopy( frame, tmp, NULL );
	cvResetImageROI( frame );	

	cvSetImageROI( frame_bk, rect);
	tmp_bk = cvCreateImage( cvGetSize(frame_bk), 8, 1 );
	cvCopy( frame_bk, tmp_bk, NULL );
	cvResetImageROI( frame_bk );

	//frame1
	unsigned long int size_object=0;

	for (int i=0;i<tmp_bk->height;i++)
	{
		for (int k=0;k<tmp_bk->width;k++)
		{
			if ((unsigned char)tmp_bk->imageData[i*tmp_bk->widthStep+k]==255)
			{
				int b = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3];
				int g = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+1];
				int r = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+2];
				if(del_color(r,g,b))
					size_object++;
			}
		}
	}
	int hsvflag=1;
	if (size_object==0)
	{
		hsvflag=0;
		for (int i=0;i<tmp_bk->height;i++)
		{
			for (int k=0;k<tmp_bk->width;k++)
			{
				if ((unsigned char)tmp_bk->imageData[i*tmp_bk->widthStep+k]==255)
				{
					int b = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3];
					int g = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+1];
					int r = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+2];
					size_object++;
				}
			}
		}
	}
	if (size_object ==0)
	{
		mcsh_region->n = 0;
		return;
	}

	CvMat *clusters;//matrix after classification
	clusters = cvCreateMat (size_object, 1, CV_32SC1);
	CvMat *points;//matrix before classification
	points = cvCreateMat (size_object, 1, CV_32FC3);

	unsigned int j=0;
	for (int i=0;i<tmp->height;i++)
	{
		for (int k=0;k<tmp->width;k++)
		{
			if ((unsigned char)tmp_bk->imageData[i*tmp_bk->widthStep+k]==255)
			{
				int b = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3];
				int g = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+1];
				int r = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+2];
				if(del_color(r,g,b)||hsvflag==0)
				{
					points->data.fl[j*3] = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3];
					points->data.fl[j*3 + 1] = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+1];
					points->data.fl[j*3 + 2] = (unsigned char) tmp->imageData[tmp->widthStep*i+k*3+2];
					j++;
				}
			}
		} 
	}//3 channel data

	cvKMeans2 (points, MAX_CLUSTERS_mcsh, clusters,
		//cvTermCriteria (CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 10, 1.0));
		cvTermCriteria (CV_TERMCRIT_EPS, 10, 0.5));
	CvMat *color = cvCreateMat (MAX_CLUSTERS_mcsh, 1, CV_32FC3);
	CvMat *count = cvCreateMat (MAX_CLUSTERS_mcsh, 1, CV_32SC1);//as counting matrix
	cvSetZero (color);
	cvSetZero (count);

	for (int i = 0; i < size_object; i++)
	{
		int idx = clusters->data.i[i];
		int j = ++count->data.i[idx];

		color->data.fl[idx * 3 ] = color->data.fl[idx * 3 ] * (j - 1) / j + points->data.fl[i * 3 ] / j;
		color->data.fl[idx * 3 + 1] = color->data.fl[idx * 3 + 1] * (j - 1) / j + points->data.fl[i * 3 + 1] / 	j;
		color->data.fl[idx * 3 + 2] = color->data.fl[idx * 3 + 2] * (j - 1) / j + points->data.fl[i * 3 + 2] / 	j;
	}
	mcsh_region->n = MAX_NUM_mcsh*4;
	mr = mcsh_region->m;
	memset( mr, 0, mcsh_region->n * sizeof(float) );

	float result[MAX_CLUSTERS_mcsh][4],temp[4];

	for (int idx=0;idx<MAX_CLUSTERS_mcsh;idx++)
	{

		result[idx][0]= (float) count->data.i[idx]/size_object;
		result[idx][1]= (float) color->data.fl[idx * 3+2 ];//note, the order is 3-2-1
		result[idx][2]= (float) color->data.fl[idx * 3+1 ];
		result[idx][3]= (float) color->data.fl[idx * 3 ];
	}
	//from large to small
	for   (int i  =  0;  i  <  MAX_CLUSTERS_mcsh;   i++)   
	{   
		for   (int j=0;j<MAX_CLUSTERS_mcsh-1-i;j++)   
		{   
			if   (result[j][0]   <   result[j+1][0])   
			{   
				for (int k=0;k<4;k++)
				{
					temp[k]=result[j][k];
					result[j][k]   =   result[j+1][k];   
					result[j+1][k] =  temp[k]; 
				}
			}   
		}   
	} 

	//get the MAX_NUM_mcsh largest colors
	for (int i=0;i<MAX_NUM_mcsh;i++)
	{
		for (int j=0;j<4;j++)
		{
			mr[i*4+j] = result[i][j];
		}
	}
	cvReleaseMat(&count);
	cvReleaseMat(&color);
	cvReleaseMat(&points);
	cvReleaseMat(&clusters);
	cvReleaseImage(&tmp_bk);
	cvReleaseImage(&tmp);
}

void CW_PMCSHROperator::cal_tl_mcsh(Tracklet_Info* tl_i)
{
	int	size_object = MAX_NUM_mcsh*tl_i->length;
	int mcsh_num=tl_i->length;

	CvMat *clusters;//matrix after classification
	clusters = cvCreateMat (size_object, 1, CV_32SC1);
	CvMat *points;//matrix before classification
	points = cvCreateMat (size_object, 1, CV_32FC3);
	CvMat *idxps;//smaple numbers before classification
	idxps = cvCreateMat (size_object, 1, CV_32FC1);

	unsigned int k=0;
	for (int i=0;i<tl_i->length;i++)
	{
		for (int s=0;s<MAX_NUM_mcsh;s++)
		{
			points->data.fl[k*3] = tl_i->f[i].m[s*4+1];
			points->data.fl[k*3 + 1] = tl_i->f[i].m[s*4+2];
			points->data.fl[k*3 + 2] = tl_i->f[i].m[s*4+3];
			idxps->data.fl[k] = tl_i->f[i].m[s*4];
			k++;
		}
	}

	cvKMeans2 (points, MAX_CLUSTERS_mcsh, clusters,
		//cvTermCriteria (CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 10, 1.0));
		cvTermCriteria (CV_TERMCRIT_EPS, 10, 0.5));
	CvMat *color = cvCreateMat (MAX_CLUSTERS_mcsh, 1, CV_32FC3);
	CvMat *count = cvCreateMat (MAX_CLUSTERS_mcsh, 1, CV_32FC1);//as counting
	cvSetZero (color);
	cvSetZero (count);

	for (int i = 0; i < size_object; i++)
	{
		int idx = clusters->data.i[i];
		count->data.fl[idx] = count->data.fl[idx] + idxps->data.fl[i];

		color->data.fl[idx * 3 ] = color->data.fl[idx * 3 ] * (count->data.fl[idx] - idxps->data.fl[i])/count->data.fl[idx] + points->data.fl[i * 3 ] * idxps->data.fl[i]/count->data.fl[idx];
		color->data.fl[idx * 3 + 1] = color->data.fl[idx * 3 + 1] * (count->data.fl[idx] - idxps->data.fl[i])/count->data.fl[idx] + points->data.fl[i * 3 + 1] * idxps->data.fl[i]/count->data.fl[idx];
		color->data.fl[idx * 3 + 2] = color->data.fl[idx * 3 + 2] * (count->data.fl[idx] - idxps->data.fl[i])/count->data.fl[idx] + points->data.fl[i * 3 + 2] * idxps->data.fl[i]/count->data.fl[idx];
	}
	tl_i->m.n = MAX_NUM_mcsh*4;
	float *mr = tl_i->m.m;
	memset( mr, 0, tl_i->m.n * sizeof(float) );

	float result[MAX_CLUSTERS_mcsh][4],temp[4];

	for (int idx=0;idx<MAX_CLUSTERS_mcsh;idx++)
	{

		result[idx][0]= (float) count->data.fl[idx]/mcsh_num;
		result[idx][1]= (float) color->data.fl[idx * 3+2 ];//note£¬the order is 3-2-1
		result[idx][2]= (float) color->data.fl[idx * 3+1 ];
		result[idx][3]= (float) color->data.fl[idx * 3 ];
	}
	//from large to small
	for   (int i  =  0;  i  <  MAX_CLUSTERS_mcsh;   i++)   
	{   
		for   (int h=0;h<MAX_CLUSTERS_mcsh-1-i;h++)   
		{   
			if   (result[h][0]   <   result[h+1][0])   
			{   
				for (int k=0;k<4;k++)
				{
					temp[k]=result[h][k];
					result[h][k]   =   result[h+1][k];   
					result[h+1][k] =  temp[k]; 
				}
			}   
		}   
	} 

	//get the MAX_NUM_mcsh largest colors
	for (int i=0;i<MAX_NUM_mcsh;i++)
	{
		for (int k=0;k<4;k++)
		{
			mr[i*4+k] = result[i][k];
		}
	}
	cvReleaseMat(&count);
	cvReleaseMat(&color);
	cvReleaseMat(&points);
	cvReleaseMat(&clusters);
	cvReleaseMat(&idxps);
}

void CW_PMCSHROperator::cal_tl_pmcsh(Tracklet_Info* tl_i)
{
	//compute the period p
	float maxsim=0;
	int t_min=20;
	int t_max=40;
	printf("id=%d\n",tl_i->id);
	if (tl_i->length>t_max)
	{
		for (int u=t_min;u<=t_max;u=u+1)
		{
			float sum_sim=0;
			int j;
			for (j=0;j<tl_i->length-u;j++)
			{
				float sim=dist_feature( &tl_i->f[j], &tl_i->f[j+u]);
				sum_sim += sim;
			}
			float usim= sum_sim/j;
			if (usim>maxsim)
			{
				maxsim=usim;
				tl_i->t=u;
			}
		}
		tl_i->p=tl_i->length/tl_i->t;
		if(tl_i->length%tl_i->t>0)tl_i->p++;
	}
	else 
	{
		tl_i->t=tl_i->length;
		tl_i->p=1;
	}

	tl_i->pf = (mcsh*)malloc(sizeof(mcsh)*tl_i->p);//release outside
	for (int j=0;j<tl_i->p;j++)
	{	
		int size_object=0;
		int mcsh_num=(j+1)*tl_i->t>tl_i->length?tl_i->length-j*tl_i->t:tl_i->t;
		size_object = MAX_NUM_mcsh*mcsh_num;

		CvMat *clusters;//matrix after classification
		clusters = cvCreateMat (size_object, 1, CV_32SC1);
		CvMat *points;//matrix before classification
		points = cvCreateMat (size_object, 1, CV_32FC3);
		CvMat *idxps;//·Ösample numbers before classification
		idxps = cvCreateMat (size_object, 1, CV_32FC1);

		unsigned int k=0;
		for (int i=tl_i->t*j;i<tl_i->t*(j+1)&&i<tl_i->length;i++)
		{
			for (int s=0;s<MAX_NUM_mcsh;s++)
			{
				points->data.fl[k*3] = tl_i->f[i].m[s*4+1];
				points->data.fl[k*3 + 1] = tl_i->f[i].m[s*4+2];
				points->data.fl[k*3 + 2] = tl_i->f[i].m[s*4+3];
				idxps->data.fl[k] = tl_i->f[i].m[s*4];
				k++;
			}
		}

		cvKMeans2 (points, MAX_CLUSTERS_mcsh, clusters,
			//cvTermCriteria (CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 10, 1.0));
			cvTermCriteria (CV_TERMCRIT_EPS, 10, 0.5));
		CvMat *color = cvCreateMat (MAX_CLUSTERS_mcsh, 1, CV_32FC3);
		CvMat *count = cvCreateMat (MAX_CLUSTERS_mcsh, 1, CV_32FC1);//as counting
		cvSetZero (color);
		cvSetZero (count);

		for (int i = 0; i < size_object; i++)
		{
			int idx = clusters->data.i[i];
			count->data.fl[idx] = count->data.fl[idx] + idxps->data.fl[i];

			color->data.fl[idx * 3 ] = color->data.fl[idx * 3 ] * (count->data.fl[idx] - idxps->data.fl[i])/count->data.fl[idx] + points->data.fl[i * 3 ] * idxps->data.fl[i]/count->data.fl[idx];
			color->data.fl[idx * 3 + 1] = color->data.fl[idx * 3 + 1] * (count->data.fl[idx] - idxps->data.fl[i])/count->data.fl[idx] + points->data.fl[i * 3 + 1] * idxps->data.fl[i]/count->data.fl[idx];
			color->data.fl[idx * 3 + 2] = color->data.fl[idx * 3 + 2] * (count->data.fl[idx] - idxps->data.fl[i])/count->data.fl[idx] + points->data.fl[i * 3 + 2] * idxps->data.fl[i]/count->data.fl[idx];
		}
		tl_i->pf[j].n = MAX_NUM_mcsh*4;
		float *mr = tl_i->pf[j].m;
		memset( mr, 0, tl_i->pf[j].n * sizeof(float) );

		float result[MAX_CLUSTERS_mcsh][4],temp[4];

		for (int idx=0;idx<MAX_CLUSTERS_mcsh;idx++)
		{

			result[idx][0]= (float) count->data.fl[idx]/mcsh_num;
			result[idx][1]= (float) color->data.fl[idx * 3+2 ];//note, the order is 3-2-1
			result[idx][2]= (float) color->data.fl[idx * 3+1 ];
			result[idx][3]= (float) color->data.fl[idx * 3 ];
		}
		//from large to small
		for   (int i  =  0;  i  <  MAX_CLUSTERS_mcsh;   i++)   
		{   
			for   (int h=0;h<MAX_CLUSTERS_mcsh-1-i;h++)   
			{   
				if   (result[h][0]   <   result[h+1][0])   
				{   
					for (int k=0;k<4;k++)
					{
						temp[k]=result[h][k];
						result[h][k]   =   result[h+1][k];   
						result[h+1][k] =  temp[k]; 
					}
				}   
			}   
		} 

		//get the MAX_NUM_mcsh largest colors
		for (int i=0;i<MAX_NUM_mcsh;i++)
		{
			for (int k=0;k<4;k++)
			{
				mr[i*4+k] = result[i][k];
			}
		}
		cvReleaseMat(&count);
		cvReleaseMat(&color);
		cvReleaseMat(&points);
		cvReleaseMat(&clusters);
		cvReleaseMat(&idxps);
	}
}

float CW_PMCSHROperator::dist_feature( mcsh * h1, mcsh * h2 )
{
	//compute the similarity
	if (h1->n == 0 || h2->n == 0)
	{
		return -1;
	}

	float score_mcsh=0;
	float score12=0;
	float score21=0;
	float sim_max=0;
	float sim_min=0;

	float distance[2];
	int c12[MAX_NUM_mcsh],c21[MAX_NUM_mcsh];
	float MCSH_Blob1[MAX_NUM_mcsh][4], MCSH_Blob2[MAX_NUM_mcsh][4];
	for (int i=0; i<MAX_NUM_mcsh; i++)
	{
		for (int j=0; j<4; j++)
		{
			MCSH_Blob1[i][j] = h1->m[i*4 + j];
			MCSH_Blob2[i][j] = h2->m[i*4 + j];
		}
	}

	distance[1]=0;
	distance[0]=0;

	for (int i=1;i<=MAX_NUM_mcsh;i++)//just the first MAX_NUM_mcsh colors
	{
		for (int j=1;j<=MAX_NUM_mcsh;j++)
		{
			distance[1]=MCSH_Color_Distance(MCSH_Blob1,MCSH_Blob2,i-1,j-1);
			if (distance[1] < MCSH_COLOR_DISTANCE_THRE)
			{
				distance[0] += MCSH_Blob2[j-1][0];
			}

		}

		score12+=min(MCSH_Blob1[i-1][0] , distance[0]);
		distance[1]=distance[0]=0;
	}

	distance[1]=0;
	distance[0]=0;

	for (int i=1;i<=MAX_NUM_mcsh;i++)//just the first MAX_NUM_mcsh colors
	{
		for (int j=1;j<=MAX_NUM_mcsh;j++)
		{
			distance[1]=MCSH_Color_Distance(MCSH_Blob2,MCSH_Blob1,i-1,j-1);
			if (distance[1] < MCSH_COLOR_DISTANCE_THRE)
			{
				distance[0] += MCSH_Blob1[j-1][0];
			}

		}
		score21+=min(MCSH_Blob2[i-1][0] , distance[0]);
		distance[1]=distance[0]=0;
	}
	float total1=0, total2=0;

	for (int i=0;i<MAX_NUM_mcsh;i++)
	{
		total1+=MCSH_Blob1[i][0];
		total2+=MCSH_Blob2[i][0];

	}
	if (total1==0||total2==0)
	{
		return -1;
	}

	sim_min=min(score12/total1,score21/total2);//the larger sim_min value£¬the more similar
	score_mcsh = sim_min;
	sim_max=max(score12,score21)/total1;

	//if(sim_max-sim_min>0.02)
	//score_mcsh=1-(sim_max-sim_min)/(sim_min+sim_max);

	return score_mcsh; //the larger, the more similar
}


float CW_PMCSHROperator::dist_mug_mcsh(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2)
{
	int i,j;
	float maxsim=0;
	float minsim=1;
	for (i=0;i<tracklet1->length;i++)
	{
		for (j=0;j<tracklet2->length;j++)
		{
			float sim=dist_feature(&tracklet1->f[i],&tracklet2->f[j]);
			if (sim<minsim)minsim=sim;
			if (sim>maxsim)maxsim=sim;
		}
	}

	return 1-(maxsim-minsim)/(maxsim+minsim);
}


float CW_PMCSHROperator::dist_mug_tl(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2)
{
	int i,j;
	float maxsim=0;
	float minsim=1;
	int count=0;
	float sim=-1;
	for (i=0;i<tracklet1->p;i++)
	{
		for (j=0;j<tracklet2->p;j++)
		{
			sim=dist_feature(&tracklet1->pf[i],&tracklet2->pf[j]);
			if (sim<0) continue;
			if (sim<minsim)minsim=sim;
			if (sim>maxsim)maxsim=sim;
			count++;
		}
	}

	if (count<=1)return sim;
	if ((maxsim-minsim)<=0.2)return maxsim;
	else return maxsim*(1-(maxsim-minsim)/(maxsim+minsim));
}

float CW_PMCSHROperator::dist_avg_tl(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2)
{
	int i,j;
	float sumsim=0;
	for (i=0;i<tracklet1->p;i++)
	{
		for (j=0;j<tracklet2->p;j++)
		{
			float sim=dist_feature(&tracklet1->pf[i],&tracklet2->pf[j]);
			sumsim += sim;
		}
	}

	if (tracklet1->p<=1)tracklet1->p=1;
	if (tracklet2->p<=1)tracklet2->p=1;
	return sumsim/(tracklet1->p*tracklet2->p);
}

float CW_PMCSHROperator::dist_mug(Tracklet_Info* tracklet1, Tracklet_Info* tracklet2)//for filtering, not for comparison
{
	int i,j;
	float maxsim=0;
	float minsim=1;
	int count=0;
	float sim=-1;
	for (i=0;i<tracklet1->p;i++)
	{
		for (j=0;j<tracklet2->p;j++)
		{
			sim=dist_feature(&tracklet1->pf[i],&tracklet2->pf[j]);
			if (sim<0) continue;
			if (sim<minsim)minsim=sim;
			if (sim>maxsim)maxsim=sim;
			count++;
		}
	}

	if (count<=1)return 0;
	return maxsim-minsim;
}

////////////////////////////////////////////PRIVATE////////////////////////////////////////////////////////

bool CW_PMCSHROperator::del_color(int r, int g,int b)
{
	bool gd=1;
	int maxc=max(max(r,g),b);
	int minc=min(min(r,g),b);
	float v=maxc/255.0;
	float s=(maxc-minc)/(float)maxc;
	if (v<0.1||s<0.1)
	{
		gd=0;
	}
	return gd;
}

float CW_PMCSHROperator::MCSH_Color_Distance(float p1[MAX_NUM_mcsh][4],float p2[MAX_NUM_mcsh][4], int i,int j)
{
	float distance=0;
	float d1,d2;
	float c1=(float)(p1[i][1]-p2[j][1]);
	c1=c1*c1;
	float c2=(float)(p1[i][2]-p2[j][2]);
	c2=c2*c2;
	float c3=(float)(p1[i][3]-p2[j][3]);
	c3=c3*c3;	
	d1=sqrt(c3+c2+c1);
	d2=sqrt(p1[i][1]*p1[i][1]+p1[i][2]*p1[i][2]+p1[i][3]*p1[i][3])+sqrt(p2[j][1]*p2[j][1]+p2[j][2]*p2[j][2]+p2[j][3]*p2[j][3]);

	distance = d1/d2;

	return distance;

}