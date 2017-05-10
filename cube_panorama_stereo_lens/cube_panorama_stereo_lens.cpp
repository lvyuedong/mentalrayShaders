/**************************************************
Source code originally comes from Kok Hwee
modified by lvyuedong
***************************************************/

#include <math.h>
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2008\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2008\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2009\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2009\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2010\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2010\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2011\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2011\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2012\shader.h"

#include "D:\projects\mentalray\mentalrayForMayaLibs\64bit\2013\shader.h"

#pragma warning (disable: 4244)

struct cube_panorama_stereo_lens {
	miInteger whichEye;
	miScalar interocular;
	miScalar zeroParallax;
	miScalar cube_width;
	miScalar cube_height;
	miScalar cube_depth;
	miBoolean   lens_correction;
	miBoolean	front_en;
	miBoolean	left_en;
	miBoolean	right_en;
	miBoolean	top_en;
	miBoolean	bottom_en;
	miBoolean	back_en;
};

extern "C" DLLEXPORT int cube_panorama_stereo_lens_version(void) {return 1;}

float cube_panorama_stereo_lens_max(float a, float b, float c){
	return ( a>=b?(a>=c?a:c):(b>=c?b:c) );
}

extern "C" DLLEXPORT miBoolean cube_panorama_stereo_lens(
  miColor *result,
  miState *state,
  struct cube_panorama_stereo_lens *paras
)
{
	miInteger eye = *mi_eval_integer(&paras->whichEye);
    miScalar ocular = *mi_eval_scalar(&paras->interocular) / 2.0;
    miScalar ocularPower = ocular * ocular;
    miScalar zero = *mi_eval_scalar(&paras->zeroParallax);

	miScalar w = *mi_eval_scalar(&paras->cube_width);
	miScalar h = *mi_eval_scalar(&paras->cube_height);
	miScalar d = *mi_eval_scalar(&paras->cube_depth);
	
	if(w<=0 || h<=0 || d<=0) return miFALSE;
	
	float w_x = w / (w+w+d+d);
	float d_x = d / (w+w+d+d);
	float d_y = d / (h+d+d);
	float h_y = h / (h+d+d);

	float max_side = cube_panorama_stereo_lens_max(w,h,d);
	float scale_w = w/max_side;
	float scale_h = h/max_side;
	float scale_d = d/max_side;
	
	miBoolean correct = *mi_eval_boolean(&paras->lens_correction);
	miBoolean f_en = *mi_eval_boolean(&paras->front_en);
	miBoolean l_en = *mi_eval_boolean(&paras->left_en);
	miBoolean r_en = *mi_eval_boolean(&paras->right_en);
	miBoolean t_en = *mi_eval_boolean(&paras->top_en);
	miBoolean bm_en = *mi_eval_boolean(&paras->bottom_en);
	miBoolean bk_en = *mi_eval_boolean(&paras->back_en);

	float x_ratio = state->raster_x / state->camera->x_resolution;
	float y_ratio = state->raster_y / state->camera->y_resolution;
	
	int side = -1;
	if( x_ratio<d_x ){
		if( y_ratio>=d_y && y_ratio<=(1-d_y) ) side = 1;
	}else if( x_ratio>=d_x && x_ratio<=(d_x+w_x) ) {
		if( y_ratio<d_y ) side = 4;
		else if( y_ratio>=d_y && y_ratio<=(1-d_y) ) side = 0;
		else if( y_ratio>(1-d_y) ) side = 3;
	}else if( x_ratio>(d_x+w_x) && x_ratio<=(1-w_x) ) {
		if( y_ratio>=d_y && y_ratio<=(1-d_y) ) side = 2;
	}else if( x_ratio>(1-w_x) ){
		if( y_ratio>=d_y && y_ratio<=(1-d_y) ) side = 5;
	}
	
	if( side == -1 ){
		result->r = result->g = result->b = result->a = 0;
		return miTRUE;
	}

	miVector raydir;
	// we use standard cube in 1 unit
	float radius = sqrt( 0.75 );
	float uval = 0;
	float vval = 0;
	switch( side ){
		case 0:	// front
				raydir.x = ((x_ratio-d_x)/w_x - 0.5f);
				raydir.y = ((y_ratio-d_y)/h_y - 0.5f);
				raydir.z = -0.5f;
				uval = M_PI + atan2(raydir.x, 0.5f);
				vval = M_PI_2 + atan2(raydir.y, sqrt(0.25f+raydir.x*raydir.x));
				/*
				if(correct){
					uval = M_PI + raydir.x * M_PI_2;
					vval = M_PI_2 + raydir.y * M_PI_2;
					rayy = -cos(vval);
					tmpy = sqrt(radius - rayy * rayy);
					rayx = -sin(uval) * tmpy;
					rayz = cos(uval) * tmpy;
					raydir.x = rayx; raydir.y = rayy; raydir.z = rayz;
				}
				*/
				break;
		case 1:	// left
				raydir.x = -0.5f;
				raydir.y = ((y_ratio-d_y)/h_y - 0.5f);
				raydir.z = -(x_ratio/d_x - 0.5f);
				uval = M_PI_2 - atan2(raydir.z, 0.5f);
				vval = M_PI_2 + atan2(raydir.y, sqrt(0.25f+raydir.z*raydir.z));
				/*
				if(correct){
					uval = M_PI_2 - raydir.z * M_PI_2;
					vval = M_PI_2 + raydir.y * M_PI_2;
					rayy = -cos(vval);
					tmpy = sqrt(radius - rayy * rayy);
					rayx = -sin(uval) * tmpy;
					rayz = cos(uval) * tmpy;
					raydir.x = rayx; raydir.y = rayy; raydir.z = rayz;
				}
				*/
				break;
		case 2:	// right
				raydir.x = 0.5f;
				raydir.y = ((y_ratio-d_y)/h_y - 0.5f);
				raydir.z = ((x_ratio-(d_x+w_x))/d_x - 0.5f);
				uval = M_PI + M_PI_2 + atan2(raydir.z, 0.5f);
				vval = M_PI_2 + atan2(raydir.y, sqrt(0.25f+raydir.z*raydir.z));
				/*
				if(correct){
					uval = M_PI + M_PI_2 + raydir.z * M_PI_2;
					vval = M_PI_2 + raydir.y * M_PI_2;
					rayy = -cos(vval);
					tmpy = sqrt(radius - rayy * rayy);
					rayx = -sin(uval) * tmpy;
					rayz = cos(uval) * tmpy;
					raydir.x = rayx; raydir.y = rayy; raydir.z = rayz;
				}
				*/
				break;
		case 3:	// top
				raydir.x = ((x_ratio-d_x)/w_x - 0.5f);
				raydir.y = 0.5f;
				raydir.z = ((y_ratio-(d_y+h_y))/d_y - 0.5f);
				uval = -atan2( raydir.x, raydir.z );
				if(uval<0.0f) uval = M_PI + M_PI + uval;
				vval = atan2(0.5f, sqrt(raydir.x*raydir.x + raydir.z*raydir.z)) + M_PI_2;
				break;
		case 4:	// bottom
				raydir.x = ((x_ratio-d_x)/w_x - 0.5f);
				raydir.y = -0.5f;
				raydir.z = -(y_ratio/d_y - 0.5f);
				uval = -atan2( raydir.x, raydir.z );
				if(uval<0.0f) uval = M_PI + M_PI + uval;
				vval = atan2(-0.5f, sqrt(raydir.x*raydir.x + raydir.z*raydir.z)) + M_PI_2;
				break;
		case 5:	// back
				raydir.x = -((x_ratio-(d_x+d_x+w_x))/w_x - 0.5f);
				raydir.y = ((y_ratio-d_y)/h_y - 0.5f);
				raydir.z = 0.5f;
				uval = -atan2(raydir.x, 0.5f);
				if(uval<0.0f) uval = M_PI + M_PI + uval;
				vval = M_PI_2 + atan2(raydir.y, sqrt(0.25f+raydir.x*raydir.x));
				/*
				if(correct){
					uval = -raydir.x * M_PI_2;
					if(uval<0.0f) uval = M_PI + M_PI + uval;
					vval = M_PI_2 + raydir.y * M_PI_2;
					rayy = -cos(vval);
					tmpy = sqrt(radius - rayy * rayy);
					rayx = -sin(uval) * tmpy;
					rayz = cos(uval) * tmpy;
					raydir.x = rayx; raydir.y = rayy; raydir.z = rayz;
				}
				*/
				break;
	}

	miVector raydir_t;
	mi_vector_from_camera(state, &raydir_t, &raydir);

	miVector rayorig, rayorig_t;
	rayorig.x = 0; rayorig.y = 0; rayorig.z = 0;
	miScalar virtualOcular = 0;
	if(zero != 0) virtualOcular = sqrt( ocularPower + pow(ocularPower/zero,2) );
	virtualOcular = ocular;
	if( eye == 0) {	//left eye
		rayorig.x = -sin(uval - M_PI_2) * virtualOcular + virtualOcular;
		rayorig.z = cos(uval - M_PI_2) * virtualOcular;
	}else if( eye == 1 ){		//right eye
	  	rayorig.x = -sin(uval + M_PI_2) * virtualOcular - virtualOcular;
		rayorig.z = cos(uval + M_PI_2) * virtualOcular;
	}
	mi_point_from_camera(state, &rayorig_t, &rayorig);
	
	if( (side==0&&f_en) ||
		(side==1&&l_en) ||
		(side==2&&r_en) ||
		(side==3&&t_en) ||
		(side==4&&bm_en) ||
		(side==5&&bk_en) ){
		return (mi_trace_eye(result, state, &rayorig_t, &raydir_t));
	}else {
		result->r = result->g = result->b = result->a = 0;
		return miTRUE;
	}

	//return miFALSE;
}
