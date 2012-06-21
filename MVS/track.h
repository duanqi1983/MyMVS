#ifndef _TRACK_HPP
#define _TRACK_HPP
#include "ModelingOptions.h"

struct Visible_Info
{
	double view_angle;
	double input_angle;
	double output_angle;
	int view_id;
	MyColor view_color;

	double New_CosAngle;
	double Lambda;
	double LightDirection[3];
	double LightLuminance[3];
	Visible_Info():view_angle(0.0),view_id(-1),Lambda(-1) {
		view_color(0) = view_color(1) = view_color(2) = 0;
		view_angle = input_angle = output_angle = 0.0;
	}
	Visible_Info(double _angle, int _id, int _r, int _g, int _b):view_angle(_angle),view_id(_id),Lambda(-1) {
		view_color(0) = _r; view_color(1) = _g; view_color(2) = _b; input_angle = view_angle/2; output_angle = view_angle/2;
	}
	Visible_Info(double _angle, double _iangle, double _oangle, int _id, int _r, int _g, int _b):view_angle(_angle),view_id(_id),Lambda(-1) {
		view_color(0) = _r; view_color(1) = _g; view_color(2) = _b; input_angle = _iangle; output_angle = _oangle;
	}
};

bool VIFLess_InputAngle(const Visible_Info &t1, const Visible_Info &t2)
{
	return t1.input_angle < t2.input_angle;
}

bool VIFLess_ViewAngle(const Visible_Info &t1, const Visible_Info &t2)
{
	return t1.view_angle < t2.view_angle;
}

class track
{
public:
	double X[5]; //3D coordinate, theta, phi spherical coordinate; r = 1
	int size; // number of pixel
	vector<int> views; // list of views
	vector<int> px; // 2D coordinate X
	vector<int> py; // 2D coordinate Y
	vector<bool> reliable;
	vector<double> view_intensity;
	double avg_view_intensity;
	int reliableViews;
	bool valid;
	double fx;

	double Norm[3];
	double PhotometricNorm[3];
	vector<Visible_Info> vec_visible_info;
	double SurfaceLambda; // that is just the lambda in lambertain model
	double Fitted_SurfaceLambda;
	double *vec_visible_color;


	track()
	{
		size = 0;
		reliableViews = 0;
		valid = true;
		views.clear();
		px.clear();
		py.clear();
		reliable.clear();
		view_intensity.clear();
		vec_visible_info.clear();
		vec_visible_color = NULL;
	}

	

	void addView(int viewId, int x, int y)
	{
		size++;
		views.push_back(viewId);
		px.push_back(x);
		py.push_back(y);
		reliable.push_back(false);
		
	}

	void addViewColor(int viewId, int r, int g, int b, double _vangle = 0.0, double _iangle = 0.0, double _oangle = 0.0)
	{
		for (int i = 0; i < this->vec_visible_info.size(); i ++) {
			if (this->vec_visible_info[i].view_id == viewId) {
				this->vec_visible_info[i].view_angle	= _vangle;
				this->vec_visible_info[i].input_angle	= _iangle; 
				this->vec_visible_info[i].output_angle	= _oangle;;
				this->vec_visible_info[i].view_color(0) = r;
				this->vec_visible_info[i].view_color(1) = g;
				this->vec_visible_info[i].view_color(2) = b;
				return;
			}
		}
		Visible_Info vif(_vangle, _iangle, _oangle, viewId, r, g, b);
		vec_visible_info.push_back(vif);
	}

	void AddViewLambda(int viewId, double _Lambda)
	{
		for (int i = 0; i < this->vec_visible_info.size(); i ++) {
			if (this->vec_visible_info[i].view_id == viewId) {
				this->vec_visible_info[i].Lambda = _Lambda;
				return;
			}
		}
		Visible_Info vif; vif.view_id = viewId; vif.Lambda = _Lambda;
		this->vec_visible_info.push_back(vif);
	}

	void AddViewLightInfo(int viewId, double* _lightdirection, double* _lightluminance)
	{
		for (int i = 0; i < this->vec_visible_info.size(); i ++) {
			if (this->vec_visible_info[i].view_id == viewId) {
				matrixCopy(_lightdirection, this->vec_visible_info[i].LightDirection, 3);
				matrixCopy(_lightluminance, this->vec_visible_info[i].LightLuminance, 3);
				return;
			}
		}
		Visible_Info vif; vif.view_id = viewId;
		matrixCopy(_lightdirection, vif.LightDirection, 3);
		matrixCopy(_lightluminance, vif.LightLuminance, 3);
		this->vec_visible_info.push_back(vif);
	}

	void removeMismatchedPixels()
	{
		int i = 0;
		while (i<size)
		{
			while (i<size && reliable[i]) i++;
			if (i>=size) break;
			while (!reliable[size-1])
			{
				size--;
				reliable.erase(reliable.begin()+size);
				views.erase(views.begin()+size);
				px.erase(px.begin()+size);
				py.erase(py.begin()+size);
			}
			if (i<size-1)
			{
				reliable[i] = reliable[size-1];
				views[i] = views[size-1];
				px[i] = px[size-1];
				py[i] = py[size-1];
				reliable[size-1] = false;
			}			
			i++;
		}		
		
	}



	~track()
	{
		views.clear();
		px.clear();
		py.clear();
		reliable.clear();
		view_intensity.clear();
		vec_visible_info.clear();
		delete vec_visible_color;
	}
};


#endif