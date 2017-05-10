#include "E:\projects\mentalray\mentalrayForMayaLibs\64bit\2014\shader.h"
#include "math.h"

static miScalar radius;

struct fisheye_lens {
	miInteger method;			//0 for linear, 1 for cosine
	miInteger lens_type;	//0 for circular, 1 for full frame
	miScalar angle_of_view;
	miScalar distort;
};

extern "C" DLLEXPORT int fisheye_lens_version(void) {return 1;}

extern "C" DLLEXPORT void fisheye_lens_init(
	miState *state,
	struct fisheye_lens *paras,
	miBoolean *inst_req
)
{
	if (!paras) *inst_req = miTRUE;
	else {
		if(*mi_eval_integer(&paras->lens_type))	
			radius = sqrt(pow((miScalar)state->camera->x_resolution,2)+pow((miScalar)state->camera->y_resolution,2))/2.0;
		else {
						if(state->camera->aspect>=1) radius = float(state->camera->y_resolution) / 2.0;
						else radius = float(state->camera->x_resolution) / 2.0;
				 }
	}
}

extern "C" DLLEXPORT miBoolean fisheye_lens(
	miColor *result,
	miState *state,
	struct fisheye_lens *paras
)
{
	miScalar focal = *mi_eval_scalar(&paras->angle_of_view);
	miScalar focal_rad_2 = focal / 180.0 * M_PI_2;
	miScalar skew = *mi_eval_scalar(&paras->distort);
	miVector org,dir,org_t,dir_t;
	miScalar x_r,y_r;
	miScalar dist,c_dist,div;
	
	x_r = state->camera->x_resolution;
	y_r = state->camera->y_resolution;
	
	org.x = 0; org.y = 0; org.z = 0;
	dir.x = state->raster_x-x_r/2.0; dir.y = state->raster_y-y_r/2.0; dir.z = 0;
	dist = mi_vector_dist(&dir,&org);
	if(dist>radius){
		result->r = result->g = result->b = result->a = 0;
		return miTRUE;
	}
	
	if(*mi_eval_integer(&paras->method)){
		c_dist = cos((1.0-dist/radius)*M_PI_2)*radius;
		if(c_dist!=0) div = dist/c_dist;
		else div = 0;
		mi_vector_mul(&dir,pow(div,skew));
		dir.z = -(x_r/2.0)/tan(focal/360.0*M_PI);
	}else {
					if(dist<=0) dir.z = -radius;
					else if(dist>=radius) dir.z = 0;
					else dir.z = -tan( pow((1-dist/radius),skew)*focal_rad_2 + M_PI_2 - focal_rad_2 ) * dist;
				}
	mi_vector_normalize(&dir);
	mi_vector_from_camera(state,&dir_t,&dir);
	mi_point_from_camera(state,&org_t,&org);
	
	return(mi_trace_eye(result,state,&org_t,&dir_t));
}
