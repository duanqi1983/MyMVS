extern viewPoint** listViewPoints;
extern int nViews;
extern ModelingOptions MOptions;

void readMiddleBuryData()
{
	FILE * myfile;
	myfile = fopen ((path+string("/")+inputFileName).c_str(),"r");

	if (myfile == NULL) cout << "Unable to open file: "<<path<<"/"<<inputFileName;

	fscanf(myfile,"%d",&nViews);
	printf("nView = %d\n",nViews);

	listViewPoints =  new viewPoint*[nViews];

	for (int i=0;i<nViews;i++)
	{		 
		char name[80];
		fscanf(myfile,"%s",name);		
		cout<<endl<<"Name: "<<name<<endl;

		double K[9],R[9],T[3];
		for (int j=0;j<9;j++) fscanf(myfile,"%lf",&K[j]);
		for (int j=0;j<9;j++) 
		{
			fscanf(myfile,"%lf",&R[j]);
			//R[j] = -R[j];
		}
		for (int j=0;j<3;j++) 
		{
			fscanf(myfile,"%lf",&T[j]);
			//T[j] = - T[j];
		}
		listViewPoints[i] = new viewPoint(i,(path+string("/")+string(name)),K,R,T);
	}	 

	// read silhouette
	for (int i=0;i<nViews;i++)
	{
		char name[80];
		if (fscanf(myfile,"%s",name)==-1) listViewPoints[i]->readSilhouette(NULL); else
			listViewPoints[i]->readSilhouette((path+(string)"/"+(string)name).c_str());
	}

	fclose(myfile);


}

void readLiJianguoData()
{
	FILE * myfile;
	myfile = fopen ((path+string("/")+inputFileName).c_str(),"r");

	if (myfile == NULL) cout << "Unable to open file: "<<path<<"/"<<inputFileName;

	fscanf(myfile,"%d",&nViews);
	printf("nView = %d\n",nViews);

	listViewPoints =  new viewPoint*[nViews];

	for (int i=0;i<nViews;i++)
	{		 
		char name[80];
		fscanf(myfile,"%s",name);		
		cout<<endl<<"Name: "<<name<<endl;

		int w,h;
		fscanf(myfile,"%d", &w);
		fscanf(myfile,"%d", &h);
		double k[5];
		for (int j=0; j<5; j++) fscanf(myfile,"%lf",&k[j]);
		double K[9],R[9],T[3],P[12];
		for (int j=0;j<9;j++) K[j]=0;
		K[0] = K[4] = k[0];
		K[2] = k[3];
		K[5] = k[4];
		K[8] = 1;
		for (int j=0;j<9;j++) fscanf(myfile,"%lf",&R[j]);
		for (int j=0;j<3;j++) fscanf(myfile,"%lf",&T[j]);
		for (int j=0;j<12;j++) fscanf(myfile,"%lf",&P[j]);
		//for (int j=0;j<8;j++) P[j]/=4;
		//getKfromPRT(R,T,P,K);

		//listViewPoints[i] = new viewPoint(i,string(name),K,R,T,P);
		listViewPoints[i] = new viewPoint(i,(path+string("/")+string(name)),K,R,T);
	}	 
	
	// read silhouette
	for (int i=0;i<nViews;i++)
	{
		char name[80];
		if (fscanf(myfile,"%s",name)==-1) listViewPoints[i]->readSilhouette(NULL); else
			listViewPoints[i]->readSilhouette((path+(string)"/"+(string)name).c_str());
	}
	fclose(myfile);


}

void readMiddleBuryData2(string foldername)
{
	if (foldername.find_last_of("\\") == foldername.length()-1) {
		foldername.erase(foldername.end()-1);
	}
	string par_filename = foldername; par_filename += "\\";
	par_filename += foldername.substr(foldername.find_last_of("\\")+1);
	par_filename += "_par.txt";
	FILE * myfile;
	myfile = fopen (par_filename.c_str(),"r");

	if (myfile == NULL) cout << "Unable to open file: "<<par_filename << endl;

	fscanf(myfile,"%d",&nViews);
	printf("nView = %d\n",nViews);

	listViewPoints =  new viewPoint*[nViews];
	vector<std::string> maskfilelist;
	GetDirFileList((foldername+"\\masks\\"), maskfilelist, MOptions.MaskFileTypes);

	for (int i=0;i<nViews;i++)
	{		
		cout << ".";
		char name[80];
		fscanf(myfile,"%s",name);	
		string img_name = foldername + "\\" + string(name);
		//cout<<endl<<"Name: "<<img_name<<endl;

		double K[9],R[9],T[3];
		for (int j=0;j<9;j++) fscanf(myfile,"%lf",&K[j]);
		for (int j=0;j<9;j++) 
		{
			fscanf(myfile,"%lf",&R[j]);
			//R[j] = -R[j];
		}
		for (int j=0;j<3;j++) 
		{
			fscanf(myfile,"%lf",&T[j]);
			//T[j] = - T[j];
		}
		listViewPoints[i] = new viewPoint(i,img_name,K,R,T);
		if (maskfilelist.size() == nViews) {
			listViewPoints[i]->readSilhouette(maskfilelist[i].c_str());
		} else {
			listViewPoints[i]->readSilhouette(NULL);
		}
	}	 
	cout << "Done"<<endl;
	fclose(myfile);


}