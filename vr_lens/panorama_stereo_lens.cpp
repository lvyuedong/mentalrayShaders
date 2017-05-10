/**************************************************
panorama lens source code derives from Kok Hwee
modified to stereo panorama lens by lvyuedong
***************************************************/

#include <math.h>
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2008\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2008\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2009\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2009\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2011\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2011\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\32bit\2012\shader.h"
//#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2012\shader.h"
#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2014\shader.h"

//#include "E:\projects\mentalray\mentalrayForMaxLibs\64bit\2009\shader.h"
//#include "E:\projects\mentalray\mentalrayForMaxLibs\32bit\2009\shader.h"
//#include "E:\projects\mentalray\mentalrayForMaxLibs\64bit\2011\shader.h"
//#include "E:\projects\mentalray\mentalrayForMaxLibs\32bit\2011\shader.h"

struct panorama_stereo_lens {
	miInteger whichEye;
	miScalar interocular;
	miScalar zeroParallax;
	miScalar horizontal_angle;
	miScalar vertical_angle;
	miBoolean	useAbsoluteAngle;
	miScalar horiz_angle_start;
	miScalar horiz_angle_end;
	miScalar vert_angle_start;
	miScalar vert_angle_end;
};

extern "C" DLLEXPORT int panorama_stereo_lens_version(void) {return 1;}

extern "C" DLLEXPORT miBoolean panorama_stereo_lens(
  miColor *result,
  miState *state,
  struct panorama_stereo_lens *paras
)
{
  miScalar uval, vval, ufactor, vfactor, rayx, rayy, rayz, tmpy;
  miVector raydir, rayorig, rayorig_t, raydir_t;
  
  miInteger eye = *mi_eval_integer(&paras->whichEye);
  miScalar ocular = *mi_eval_scalar(&paras->interocular) / 2.0;
  miScalar ocularPower = ocular * ocular;
  miScalar zero = *mi_eval_scalar(&paras->zeroParallax);
  miScalar hangle = *mi_eval_scalar(&paras->horizontal_angle);
  miScalar vangle = *mi_eval_scalar(&paras->vertical_angle);
  miBoolean abAngle = *mi_eval_boolean(&paras->useAbsoluteAngle);
  miScalar h_start = *mi_eval_scalar(&paras->horiz_angle_start);
  miScalar h_end = *mi_eval_scalar(&paras->horiz_angle_end);
  miScalar v_start = *mi_eval_scalar(&paras->vert_angle_start);
  miScalar v_end = *mi_eval_scalar(&paras->vert_angle_end);
  miScalar virtualOcular = 0;

  if(!abAngle){
	ufactor = (hangle / 360.0) * (M_PI + M_PI);
	vfactor = (vangle / 180.0) * M_PI;
	uval = ufactor * (state->raster_x / state->camera->x_resolution) + (M_PI + M_PI - ufactor)*0.5;
	vval = vfactor * (state->raster_y / state->camera->y_resolution) + (M_PI - vfactor)*0.5;
  }else {
			if( h_end <= h_start || v_end <= v_start) return(miFALSE);
			ufactor = ( (h_end - h_start) / 360.0) * (M_PI + M_PI);
			vfactor = ( (v_end - v_start) / 180.0) * M_PI;
			uval = ufactor * (state->raster_x / state->camera->x_resolution) + h_start / 180.0 * M_PI;
			vval = vfactor * (state->raster_y / state->camera->y_resolution) + v_start / 180.0 * M_PI;
	}

	rayy = -cos(vval);				//right hand coordination system in mental ray
	tmpy = sqrt(1 - rayy * rayy);
	rayx = -sin(uval) * tmpy;
	rayz = cos(uval) * tmpy;
  raydir.x = rayx; raydir.y = rayy; raydir.z = rayz;
  
  if(zero != 0)
	  virtualOcular = sqrt( ocularPower + pow(ocularPower/zero,2) );
  
  rayorig.x = 0; rayorig.y = 0; rayorig.z = 0;
  if( eye == 0) {	//left eye
			rayorig.x = -sin(uval - M_PI_2) * virtualOcular + virtualOcular;
			rayorig.z = cos(uval - M_PI_2) * virtualOcular;
  }else if( eye == 1 ){		//right eye
  		rayorig.x = -sin(uval + M_PI_2) * virtualOcular - virtualOcular;
			rayorig.z = cos(uval + M_PI_2) * virtualOcular;
  }
  
  mi_vector_from_camera(state, &raydir_t, &raydir);
  mi_point_from_camera(state, &rayorig_t, &rayorig);

  return (mi_trace_eye(result, state, &rayorig_t, &raydir_t));
}
