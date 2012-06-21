int testImage1, testImageD;

void mouseHandler(int event, int x, int y, int flags, void* param)
{
	IplImage* image = (IplImage*) param;
	IplImage* img = cvCloneImage( image );
	if (event == CV_EVENT_MOUSEMOVE)
	{
		//printf("\nClicked at y=%d x=%d \n", y, x);
		int x1 = y;
		int y1 = x;
		for (int x1=0;x1<listViewPoints[testImage1]->height; x1++)
		{
			int x2 = listViewPoints[testImage1]->originalMatchingResultX[testImageD][x1][y1];
			int y2 = listViewPoints[testImage1]->originalMatchingResultY[testImageD][x1][y1];
			if (x2>=0) 
			{
				cvLine(img, cvPoint(y1,x1), cvPoint(y2+listViewPoints[testImage1]->width,x2), 
				cvScalar(
				128*(x1%2),
				128*((x1+1)%2)
				,0),1);
			}
		}
		x1 = y;
		int x2 = listViewPoints[testImage1]->originalMatchingResultX[testImageD][x1][y1];
		int y2 = listViewPoints[testImage1]->originalMatchingResultY[testImageD][x1][y1];
		if (x2>=0) 
			cvLine(img, cvPoint(y1,x1), cvPoint(y2+listViewPoints[testImage1]->width,x2), 
				cvScalar(0,0,255),1);

		cvShowImage("mainWin", img );
	}
	cvReleaseImage(&img );

}

void testMatchingResult(int imageId1, int d)
{
	int imageId2 = listViewPoints[imageId1]->stereoPairs[d];
	testImage1 = imageId1;
	testImageD = d;
	cvNamedWindow("mainWin", CV_WINDOW_AUTOSIZE); 

	int height = listViewPoints[imageId1]->height;
	int width = listViewPoints[imageId1]->width;
	IplImage* img=cvCreateImage(cvSize(width*2,height),IPL_DEPTH_8U,3);
	CvScalar s;

	for (int i=0;i<height;i++)
	{
		for (int j=0;j<width;j++) 
		{
			s=cvGet2D(listViewPoints[imageId1]->image,i,j);
			s.val[1] = s.val[2] = s.val[0];
			cvSet2D(img,i,j,s); // set the (i,j) pixel value
		}

		for (int j=width;j<width+width;j++) 
		{
			s=cvGet2D(listViewPoints[imageId2]->image,i,j-width);
			s.val[1] = s.val[2] = s.val[0];
			cvSet2D(img,i,j,s); // set the (i,j) pixel value
		}
	}



	// show the image
	cvShowImage("mainWin", img );

	cvSetMouseCallback("mainWin",mouseHandler,(void*) img);

	// wait for a key
	cvWaitKey(0);

	// release the image
	cvReleaseImage(&img );

	cvDestroyWindow("mainWin");


}