// RunEGTracker.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "EGTracker.h"

#define GroundTruthNum 100000

void WriteOutput(Boundingbox* output, int length, FILE *fp)
{
	int i;
	for(i=0;i<length;i++)
	{
		fprintf(fp,"%d,%d,%d,%d,%d,%d,%d\n",output[i].nCam,output[i].nFrm,output[i].ObjID,output[i].r.x,output[i].r.y,output[i].r.width,output[i].r.height);
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	if (argc!=7)
	{
		printf("There should be six input params.\nFirst five inputs are the pathes of Param.ini,Video.txt,Groundtruth.txt, Topology.txt and Area.txt.\nThe last path is the output path.\nNote: they should be inputed orderly!");
		return -1;
	}
	printf("%d\n%s\n%s\n%s\n%s\n%s\n%s\n",argc,argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
	CW_EGTracker tracker;
    tracker.Initial(0);
	tracker.LoadParam(argv[1]);
	tracker.LoadData(argv[2],argv[3],argv[4],argv[5]);
	tracker.Process(0);
	int length = 0;
	Boundingbox* output;
	output = (Boundingbox*)malloc(sizeof(Boundingbox)*GroundTruthNum);
	tracker.GetResult(output,&length);
	tracker.Release();

	char optpath[MAX_PATH_LENGTH];
	sprintf(optpath,"%s/output.txt",argv[6]);
	FILE *outputfile = fopen(optpath, "w");
	WriteOutput(output, length, outputfile);
	fclose(outputfile);
	free(output);

	return 0;
}

