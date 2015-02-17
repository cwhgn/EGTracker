#pragma once

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <stdio.h>

#define INF 1000000
#define INF_MIN 0.000001

#define MAX_NUM_mcsh 25

struct mcsh {
	float m[MAX_NUM_mcsh * 4];
	int n;                     
};

struct Tracklet_Info
{
	CvRect* r;
	int* frm;
	int nCam;
	int id;
	int length;
	float conv;

	int t;//the length of period
	int p;//the number of periods
	mcsh* pf;
	mcsh* f;
	int f_free;
	mcsh m;
};