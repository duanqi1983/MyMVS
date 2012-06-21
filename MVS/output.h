void outputVerticesWithNormals(vector<track> tracks, string filename)
{
	//string filename = path +"WithNormals.ply";
	printf("\nWriting to %s file\n",filename.c_str());
	FILE * pFile;
	pFile = fopen(filename.c_str(),"w");
	fprintf(pFile,"ply\n");
	fprintf(pFile,"format ascii 1.0\n");
	int validPoints = 0;
	for (int i=0;i<tracks.size();i++) {
		if (tracks[i].valid) 
			validPoints++;
	}
	fprintf(pFile,"element vertex %d\n", validPoints);
	fprintf(pFile,"property float x\n");
	fprintf(pFile,"property float y\n");
	fprintf(pFile,"property float z\n");
	fprintf(pFile,"property float nx\n");
	fprintf(pFile,"property float ny\n");
	fprintf(pFile,"property float nz\n");
	fprintf(pFile,"end_header\n");	
	for (int i=0;i<tracks.size();i++) {
		if (tracks[i].valid)
			fprintf(pFile,"%9.9lf %9.9lf %9.9lf %9.9lf %9.9lf %9.9lf \n",
			tracks[i].X[0],tracks[i].X[1],tracks[i].X[2],
			sin(tracks[i].X[3])*cos(tracks[i].X[4]),
			sin(tracks[i].X[3])*sin(tracks[i].X[4]),
			cos(tracks[i].X[3]));
	}  
	fclose(pFile);
}


void writeToOBJ(vector<track> tracks, string filename)
{
	printf("\nWriting to %s file\n",filename.c_str());
	FILE * pFile;
	pFile = fopen(filename.c_str(),"w");
	for (int i=0;i<tracks.size();i++) 
		if (tracks[i].valid)
			fprintf(pFile,"v %9.9lf %9.9lf %9.9lf \n",
			tracks[i].X[0],tracks[i].X[1],tracks[i].X[2]);
	fclose(pFile);
}


void writeToNPTS(vector<track> tracks, string filename)
{
	printf("\nWriting to %s file\n",filename.c_str());
	FILE * pFile;
	pFile = fopen(filename.c_str(),"w");
	for (int i=0;i<tracks.size();i++) {
		if (tracks[i].valid)
			fprintf(pFile,"%9.9lf %9.9lf %9.9lf %9.9lf %9.9lf %9.9lf \n",
			tracks[i].X[0],tracks[i].X[1],tracks[i].X[2],
			sin(tracks[i].X[3])*cos(tracks[i].X[4]),
			sin(tracks[i].X[3])*sin(tracks[i].X[4]),
			cos(tracks[i].X[3]));
	}  
	fclose(pFile);
}

void writeToNPTS2(vector<track> tracks, string filename)
{
	cout << "Writing to " << filename;
	fstream fout(filename, ios::out);
	if (!fout) {
		cout << "Can not open " << filename << " to save data..." << endl;
		return;
	}
	fout << setprecision (8);
	for (int i=0;i<tracks.size();i++) {
		if (tracks[i].valid) {
			fout << tracks[i].X[0] << " " << tracks[i].X[1] << " " << tracks[i].X[2] << " "
				<< sin(tracks[i].X[3])*cos(tracks[i].X[4]) << " " << sin(tracks[i].X[3])*sin(tracks[i].X[4]) << " "
				<< cos(tracks[i].X[3]) << endl;
		}	
	}  
	fout.close();
	cout << "Done..." << endl;
}

void outputPara(char * fileName)
{
	printf("\nWrting parameters to %s file\n", fileName);
	FILE *pFile;
	pFile = fopen((path + "/" + (string)fileName).c_str(), "w");
	fprintf(pFile,"%d\n", nViews);
	for (int i=0;i<nViews;i++)
	{
		fprintf(pFile,"%s\n",listViewPoints[i]->imageName.c_str());
		fprintf(pFile,"%d %d\n",listViewPoints[i]->width,listViewPoints[i]->height);
		fprintf(pFile,"%.9lf %.9lf 0 %.9lf %.9lf\n",
			listViewPoints[i]->K[0],listViewPoints[i]->K[4],
			listViewPoints[i]->K[2],listViewPoints[i]->K[5]);

		for (int j=0;j<9;j++) fprintf(pFile,"%.9lf ", listViewPoints[i]->R[j]);
		fprintf(pFile,"\n");
		for (int j=0;j<3;j++) fprintf(pFile,"%.9lf ", listViewPoints[i]->T[j]);
		fprintf(pFile,"\n");
		for (int j=0;j<12;j++) fprintf(pFile,"%.9lf ", listViewPoints[i]->P[j]);
		fprintf(pFile,"\n\n");

	}
	fclose(pFile);

}

void outputPara2(char * fileName, double scale)
{
	printf("\nWrting parameters to %s file\n", fileName);
	FILE *pFile;
	pFile = fopen((path + "/" + (string)fileName).c_str(), "w");
	fprintf(pFile,"%d\n", nViews);
	for (int i=0;i<nViews;i++)
	{
		fprintf(pFile,"%s ",listViewPoints[i]->imageName.c_str());
		//fprintf(pFile,"%.9lf %.9lf 0 %.9lf %.9lf\n",listViewPoints[i]->K.at<double>(0,0),listViewPoints[i]->K.at<double>(1,1),listViewPoints[i]->K.at<double>(0,2),listViewPoints[i]->K.at<double>(1,2));
		for (int ii=0;ii<2;ii++)
		{
			for (int jj=0;jj<3;jj++) fprintf(pFile,"%.9lf ", scale*listViewPoints[i]->K[ii*3+jj]);
			//fprintf(pFile,"\n");
		}
		fprintf(pFile,"%.9lf %.9lf %.9lf ",0.,0.,1.);
		for (int ii=0;ii<3;ii++)
		{
			for (int jj=0;jj<3;jj++) fprintf(pFile,"%.9lf ", listViewPoints[i]->R[ii*3+jj]);
			//fprintf(pFile,"\n");
		}
		//fprintf(pFile,"\n");
		for (int j=0;j<3;j++) fprintf(pFile,"%.9lf ", listViewPoints[i]->T[j]);
		fprintf(pFile,"\n");
		/*
		for (int ii=0;ii<3;ii++)
		{
			for (int jj=0;jj<4;jj++) fprintf(pFile,"%.9lf ", listViewPoints[i]->P.at<double>(ii,jj));
			//fprintf(pFile,"\n");
		}
		fprintf(pFile,"\n\n");
		*/

	}
	fclose(pFile);

}

bool PoissonReconstruction(string inputfile, string &output_prefix, int depth_value = 12, int samplenode = 4)
{
	string PRExeFile = "PoissonRecon.64.exe";
	if (!FileExisted(PRExeFile.c_str())) {
		std::cerr << "Please put the PoissonRecon.exe file into the current folder ..." << std::endl;
		return false;
	}
	char buffer[10];
	output_prefix += "_";itoa(depth_value, buffer, 10); output_prefix += buffer;
	output_prefix += "_";itoa(samplenode, buffer, 10); output_prefix += buffer;
	output_prefix += ".ply";
	string cmd = PRExeFile + " --in ";
	cmd += inputfile;
	cmd += " --out ";
	cmd += output_prefix;
	cmd += " --depth ";
	itoa(depth_value, buffer, 10); cmd += buffer;
	cmd += " --samplesPerNode ";
	itoa(samplenode, buffer, 10); cmd += buffer;
	std::cout << std::endl << cmd << endl;
	WinExec(cmd.c_str(),0);
	while(!FileExisted(output_prefix.c_str()))
		::Sleep(5000);
}