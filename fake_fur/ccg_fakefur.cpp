#include <math.h>
#include "..\..\mentalrayForMayaLibs\64bit\2014\shader.h"
#include "..\..\mentalrayForMayaLibs\64bit\2014\mayaapi.h"
#include "crystalcg.h"

miBoolean ccg_illegalTex(
	miVector t1,
	miVector t2,
	miVector t3
	)
{
	miScalar tex1,tex2,tex3;

	if(t1.x==t1.y&&t1.x==t1.z&&t1.x==0) tex1 = 0;
		else tex1 = 1;
	if(t2.x==t2.y&&t2.x==t2.z&&t2.x==0) tex2 = 0;
		else tex2 = 1;
	if(t3.x==t3.y&&t3.x==t3.z&&t3.x==0) tex3 = 0;
		else tex3 = 1;
			
	if(tex1==0&&tex2==0) return miFALSE;
		else if(tex2==0&&tex3==0)	return miFALSE;
			else if(tex1==0&&tex3==0)	return miFALSE;
				else return miTRUE;
}


miBoolean ccg_getTangentUV(
	miState		*state,
	miInteger	idx,
	miVector 	t1,
	miVector 	t2,
	miVector 	t3,
	miVector 	p1,
	miVector 	p2,
	miVector 	p3,
	miVector 	*out_u,
	miVector 	*out_v
	)
{
	miVector	vc,vcA,vcB;
	miScalar	fu21,fv21,fu31,fv31;
	int num = 0;
  
	mi_query( miQ_NUM_BUMPS, state, miNULLTAG, &num);

	if ( num > idx && state->bump_x_list != NULL )
	{
		*out_u = state->bump_x_list[idx];
		*out_v = state->bump_y_list[idx];
		mi_vector_from_object(state, out_u, out_u);
		mi_vector_from_object(state, out_v, out_v);
		mi_vector_normalize(out_u);
		mi_vector_normalize(out_v);
		return miTRUE;
	}

	fu21 = t2.x - t1.x;
	fv21 = t2.y - t1.y;
	fu31 = t3.x - t1.x;
	fv31 = t3.y - t1.y;

	vcA.x = p2.x - p1.x;	vcA.y = fu21;		vcA.z = fv21;
	vcB.x = p3.x - p1.x;	vcB.y = fu31;		vcB.z = fv31;
	mi_vector_prod(&vc, &vcA, &vcB);
	if(vc.x!=0.0)
	{
		out_u->x = -vc.y/vc.x;
		out_v->x = -vc.z/vc.x;
	}

	vcA.x = p2.y - p1.y;	vcA.y = fu21;		vcA.z = fv21;
	vcB.x = p3.y - p1.y;	vcB.y = fu31;		vcB.z = fv31;
	mi_vector_prod(&vc, &vcA, &vcB);
	if(vc.x!=0.0)
	{
		out_u->y = -vc.y/vc.x;
		out_v->y = -vc.z/vc.x;
	}

	vcA.x = p2.z - p1.z;	vcA.y = fu21;		vcA.z = fv21;
	vcB.x = p3.z - p1.z;	vcB.y = fu31;		vcB.z = fv31;
	mi_vector_prod(&vc, &vcA, &vcB);
	if(vc.x!=0.0)
	{
		out_u->z = -vc.y/vc.x;
		out_v->z = -vc.z/vc.x;
	}

	mi_vector_normalize(out_u);
	mi_vector_normalize(out_v);

	return miTRUE;
}



miScalar ccg_SC_specular_blinn(
	miScalar ecc,
	miScalar roll,
	miVector N,
	miVector L,
	miVector V,
	miScalar LN,
	miScalar VN ) 
{
	miVector H;
    float NH, NH2, NHSQ, Dd, Gg, VH, Ff, tmp;

	mi_vector_normalize(&L);
    mi_vector_add(&H, &V, &L);
	mi_vector_normalize(&H);
	//LN = mi_vector_dot(&L, &N);
    NH = mi_vector_dot(&N, &H);
    NHSQ = NH*NH;
    NH2 = NH * 2.0f;
    Dd = (ecc+1.0f) / (NHSQ + ecc);
    Dd *= Dd;
    VH = mi_vector_dot(&V, &H);
    if( VN < LN )
    {
		if( VN * NH2 < VH )
			Gg = NH2 / VH;
		else
			Gg = 1.0f / VN;
    }
    else
    {
		if( LN * NH2 < VH )
			Gg = (LN * NH2) / (VH * VN);
		else
			Gg = 1.0f / VN;
    }
    /* poor man's Fresnel */
    tmp = pow((1.0f - VH), 3.0f);
    Ff = tmp + (1.0f - tmp) * roll;
	return (Dd * Gg * Ff);
}


struct ccg_fakefur_t
{
     miColor	    skinColor;
	 miColor		skinDiffuseMapMult;
     miColor	    hairColor;
	 miColor		hairDiffuseMapMult;
     miColor		incandescence;
     miVector 		normalMapping;
     //hair shading
     miColor      	Kajiya_specColor1;
     miScalar  		Kajiya_specExp1;
     miColor      	Kajiya_specColor2;
     miScalar  		Kajiya_specExp2;
	 miColor		Goldman_specColor;
	 miScalar		Goldman_specPower;
	 miColor		hairSpecularMapMult;
	 miScalar		hairDirection;
	 miInteger		hairSpecularMode;
	 miInteger		hairDiffuseMode;
	 //hair parameters
     miScalar       hairDensity;
	 miScalar		hairDensityMapMult;
     miScalar  		backwardScattering; // hair backward scattering
     miScalar  		forwardScattering; // hair forward scattering
     miScalar 		hairOverHairShadow; // density of the hair-over-hair shadows
     miScalar  		shadowAttenuation; // hair over surface shadow attenuation
     miScalar     	startShadowCos; // cosine of the starting shadow termination
     miScalar       endShadowCos; // cosine of the ending shadow termination
     miScalar       hairLength; // hair length
     miScalar     	hairBaseRadius; // hair radius at base
     miScalar      	hairTipRadius; // hair radius at tip
														     //miScalar  tangentBreakupAmt;
														     //miScalar     tangentBreakup;
     //skin shading
     miColor		skinSpecularColor;
     miScalar		eccentricity;
     miScalar		specularRolloff;
	 miColor		skinSpecularMapMult;
     miBoolean		disableSkinSpecular;
     miBoolean		directBlinnSpecular;
	 
	 //sss components
	 miColor		sss;
	 miColor		sssMapMult;
     
	 miBoolean		surface;
     int            mode;
     int            i_lights;
     int	        n_lights;
     miTag          lights[1];
};


extern "C" DLLEXPORT int ccg_fakefur_version(void) { return 1;}


miScalar sinBetweenTwoVector(miVector a, miVector b)
{
	miVector c;
	mi_vector_prod(&c, &a, &b);
	return (mi_vector_norm(&c));
}

// fakefurOpacity 
// Computes the mean opacity of a patch of fur as viewed from a given angle 
// Probability of a random ray striking the single hair from direction V
miScalar fakefur_opacity (miVector *V, miScalar An, miScalar D,
						miVector *T, miVector *NN)
{
   miScalar g, alpha;
   
   // we assume that V and T are unit vectors 
   miVector B;
   mi_vector_prod(&B, V, T);
   alpha = mi_vector_dot(V, NN);
   g = mi_vector_norm(&B);
   if ( alpha > miSCALAR_EPSILON ) g /= alpha;
   
   // Opacity is an inverted Gaussian
   alpha = exp( D * An * g );
   if ( alpha <= miSCALAR_EPSILON ) return 0.0f;
   alpha = 1 - ( 1 / alpha ); 
   return alpha;
}

/*
vector shiftTangent( const vector& T, const vector& n, const miScalar amt )
{
   vector shiftedT = T + n * amt;
   shiftedT.normalize();
   return shiftedT;
}
*/


miScalar fakefur_smoothstep(miScalar min, miScalar max, miScalar x)
{
	if(x<min) return 0.0f;
	else if(x>max) return 1.0f;
	else {
		miScalar t = (x - min) / (max - min);
		return (t*t*(3.0f - 2.0f*t));
	}
}

miScalar fakefur_strandSpecular( miVector *T, miVector *H, miScalar exponent)
{
   miScalar TdH = mi_vector_dot(T, H);
   miScalar sinTH = sqrt( 1.0f - TdH * TdH );
   miScalar dirAtten = fakefur_smoothstep(-1.0f, 0.0f, TdH);
   return ( dirAtten * pow(sinTH, exponent ) );
}


extern "C" DLLEXPORT miBoolean ccg_fakefur(
	miColor      *result,
	miState     *state,
	ccg_fakefur_t*  p)
{
   if ( state->type == miRAY_SHADOW ) return miFALSE;

   //miScalar v = state->bary[1];
   float hair_density =   *mi_eval_scalar( &p->hairDensity );
   miScalar	hair_density_map = *mi_eval_scalar( &p->hairDensityMapMult );
   hair_density *= hair_density_map;
   //const miScalar bu           =   mr_eval( p->tangentBreakupAmt );
   //const miScalar bup          =   mr_eval( p->tangentBreakup ) * v;

   miColor skinColor      =   *mi_eval_color( &p->skinColor );
   miColor	skinDiffMap		=  *mi_eval_color( &p->skinDiffuseMapMult );
   miColor hairColor        =   *mi_eval_color( &p->hairColor );
   miColor	hairDiffMap		=	*mi_eval_color( &p->hairDiffuseMapMult );
   miColor specColor1          =   *mi_eval_color( &p->Kajiya_specColor1 );
   miScalar exp1         =   *mi_eval_scalar( &p->Kajiya_specExp1 );
   miColor specColor2          =   *mi_eval_color( &p->Kajiya_specColor2 );
   miScalar exp2         =   *mi_eval_scalar( &p->Kajiya_specExp2 );
   miColor	specColor3		= *mi_eval_color(&p->Goldman_specColor);
   miScalar specPower	= *mi_eval_scalar(&p->Goldman_specPower);
   miColor	hair_spec_map = *mi_eval_color( &p->hairSpecularMapMult );
   miInteger hairShadingMode	= *mi_eval_integer(&p->hairSpecularMode);
   miInteger hairDiffMode	= *mi_eval_integer(&p->hairDiffuseMode);
   miScalar p_reflect    =   *mi_eval_scalar( &p->backwardScattering );
   miScalar p_transmit   =   *mi_eval_scalar( &p->forwardScattering );
   miScalar p_surface    =   *mi_eval_scalar( &p->shadowAttenuation );
   miScalar s            =   *mi_eval_scalar( &p->hairOverHairShadow );
   miScalar w_min        =   *mi_eval_scalar( &p->startShadowCos );
   miScalar w_max        =   *mi_eval_scalar( &p->endShadowCos );
   miScalar len          =   *mi_eval_scalar( &p->hairLength );
   miScalar rb           =   *mi_eval_scalar( &p->hairBaseRadius );
   miScalar rt           =   *mi_eval_scalar( &p->hairTipRadius );
   
  miBoolean	*mayaStateDiffuse = NULL, *mayaStateSpecular = NULL;
	miBoolean	emitDiffuse = miTRUE;		/* default emits diffuse */
	miBoolean	emitSpecular = miTRUE;	/* default emits specular */
  
   miVector normalbend = *mi_eval_vector(&p->normalMapping);

   // normalized eye direction vector
   miVector V = state->dir;
   mi_vector_neg(&V);
   mi_vector_normalize(&V);

   //
   // TODO: get hair normals/tangents from Maya's Fur export data files
   //

   // normalized hair tangent vector
   miVector T;
   //if(*mi_eval_boolean(&p->surface))
		//T = state->bump_x_list[0];
   //else T = state->tex_list[0];
		const miVector *t[3], *po[3];
		miVector	u_dir, v_dir;
		if(mi_tri_vectors(state, 't', 0, &t[0], &t[1], &t[2]))
		{
			if(ccg_illegalTex(*t[0],*t[1],*t[2])&&mi_tri_vectors(state, 'p', 0, &po[0], &po[1], &po[2]))
			{
				ccg_getTangentUV(state, 0, *t[0], *t[1], *t[2], *po[0], *po[1], *po[2], &u_dir, &v_dir);
				miScalar xangle = *mi_eval_scalar(&p->hairDirection);
				if(fmod(xangle,360)!=0.0)
				{
					miMatrix rotationU;
					mi_matrix_rotate_axis(rotationU, &state->normal, xangle/360.0f*float(CCG_PI));
					mi_vector_transform(&u_dir, &u_dir, rotationU);
					mi_vector_transform(&v_dir, &v_dir, rotationU);
					mi_vector_normalize(&u_dir);
					mi_vector_normalize(&v_dir);
				}
			}
		}else {
				ccg_vector_init(&u_dir, 0);
				ccg_vector_init(&v_dir, 0);
			  }
	 T = u_dir;
   //mi_vector_normalize(&T);
   
   miVector Nn = state->normal;
   mi_vector_normalize(&Nn);
   miVector hairN = Nn;
   miVector skinN = Nn;

   // hack: shift T a little bit
   //T = shiftTangent( T, hairN, bu * snoise(bup) );
   
   float hair_area = ( len * ( rb + rt ) ) / 2.0f;
   
   // Cross product of the tangent along the hair and view vector
   miVector TcV;
   mi_vector_prod(&TcV, &T, &V);
   float sinTV = mi_vector_norm(&TcV);

	// compute the hair/skin visibility ratio
	float hair_skin_visibility = fakefur_opacity( &V, hair_area, hair_density, &T, &hairN );
	 
	miColor skinSpecColor = *mi_eval_color(&p->skinSpecularColor);
	float eccen = *mi_eval_scalar(&p->eccentricity);
	if (fabs(eccen) == 1.0f) eccen = 0.999999f;
	eccen = 1.0f / (pow(eccen, 2) - 1.0f);
		
	 float rolloff = *mi_eval_scalar(&p->specularRolloff);
	 miColor	skin_spec_map = *mi_eval_color( &p->skinSpecularMapMult );
	 miBoolean disableSkinSpec = *mi_eval_boolean(&p->disableSkinSpecular);
     
   int m = *mi_eval_integer(&p->mode);
   int n_l = *mi_eval_integer(&p->n_lights);
   int i_l = *mi_eval_integer(&p->i_lights);
   miTag *light = mi_eval_tag(&p->lights);
   if (m == 1)   /* modify light list (inclusive mode) */
			mi_inclusive_lightlist(&n_l, &light, state);
		else if (m == 2)  /* modify light list (exclusive mode) */
			mi_exclusive_lightlist(&n_l, &light, state);
		else if(m == 4){
			n_l = 0;
			light = 0;
		}

   miColor Ci;
   Ci.r = Ci.g = Ci.b = Ci.a = 0;
   miColor Cl;
   Cl.r = Cl.g = Cl.b = Cl.a = 0;
   miScalar NdL = 0;
   miVector L;
   L.x = L.y = L.z = 0;

   if (m == 4 || n_l){
		for (mi::shader::LightIterator iter(state, light, n_l);!iter.at_end(); ++iter) {
			miColor Cs;
			Cs.r = Cs.g = Cs.b = Cs.a = 0;
			int	samples = 0;
			while (iter->sample()) {
				//get custom maya light properties
				if (mayabase_stateitem_get(state,
						MBSI_LIGHTDIFFUSE, &mayaStateDiffuse,
						MBSI_LIGHTSPECULAR, &mayaStateSpecular,
						MBSI_NULL)) {
					emitDiffuse = *mayaStateDiffuse;
					emitSpecular = *mayaStateSpecular;
				}
						
				iter->get_contribution(&Cl);
				NdL = iter->get_dot_nl();
				L = iter->get_direction();

				miVector Ln = L;
				mi_vector_normalize(&Ln);
				
				// scattering of light by fur
				miVector TcL;
				 mi_vector_prod(&TcL, &T, &Ln); 
				 float sinTL = mi_vector_norm(&TcL);
				 float hairNdL = mi_vector_dot(&hairN, &Ln);
				 float skinNdL = hairNdL;
			    
				float opacL = fakefur_opacity( &L, hair_area, hair_density, &T, &hairN );
				// using the fakefur opacity function, compute the hair-over-hair
				 // shadow attenuation
				 float hair_over_hair = 1 - s * opacL;
				 // scattering of light by fur
				 float k = mi_vector_dot(&TcL, &TcV) / ( sinTL * sinTV );
				 // directional attenuation factor
				 float f_dir = ((1 + k) / 2 * p_reflect) + ((1 - k) / 2 * p_transmit);
				// Surface normal factor as a quick and dirty way to adjust shadowing
				 float f_surface = 1 + p_surface * ( fakefur_smoothstep( w_min, w_max, hairNdL) - 1 );
				// using the fakefur opacity function, compute the
				// hair-over-skin shadow factor
				float hair_over_skin_shadow = 1 - opacL;
				float alpha_hair = hair_skin_visibility * hair_over_hair * f_dir * f_surface;
				float alpha_skin = (1 - hair_skin_visibility) * hair_over_skin_shadow;

				//specular
				miColor ks;
				ks.r = ks.g = ks.b = 0;
				
				if(emitSpecular)
				{
					 miColor khs;		//hair specular factor
					 khs.r = khs.g = khs.b = 0;
					 miVector H;
					 mi_vector_add(&H, &V, &L);
					 mi_vector_normalize(&H);
					 
					 	if(hairShadingMode==0) {
					 		float tmp1 = fakefur_strandSpecular(&T, &H, exp1);
					 		float tmp2 = fakefur_strandSpecular(&T, &H, exp2);
							khs.r = tmp1 * specColor1.r + tmp2 * specColor2.r;
							khs.g = tmp1 * specColor1.g + tmp2 * specColor2.g;
							khs.b = tmp1 * specColor1.b + tmp2 * specColor2.b;
						}
						else if(hairShadingMode==1) {
										float tmp = pow( (mi_vector_dot(&T, &Ln)*mi_vector_dot(&T, &V)+sinTL*sinTV), specPower);
										khs.r = tmp * specColor3.r;
										khs.g = tmp * specColor3.g;
										khs.b = tmp * specColor3.b;
									}
						else {
										float tmp1 = fakefur_strandSpecular(&T, &H, exp1);
					 					float tmp2 = fakefur_strandSpecular(&T, &H, exp2);
					 					float tmp = pow( (mi_vector_dot(&T, &Ln)*mi_vector_dot(&T, &V)+sinTL*sinTV), specPower);
										khs.r = tmp1 * specColor1.r + tmp2 * specColor2.r;
										khs.g = tmp1 * specColor1.g + tmp2 * specColor2.g;
										khs.b = tmp1 * specColor1.b + tmp2 * specColor2.b;
										khs.r += tmp * specColor3.r;
										khs.g += tmp * specColor3.g;
										khs.b += tmp * specColor3.b;
									}
						khs.r *= hair_spec_map.r;
						khs.g *= hair_spec_map.g;
						khs.b *= hair_spec_map.b;
						
						miColor kss;		//skin specular factor
						kss.r = kss.g = kss.b = 0;
						if(!disableSkinSpec){
								float VdotN = mi_vector_dot(&V,&Nn);
								float blinn_f = ccg_SC_specular_blinn(eccen, rolloff, Nn, L, V, NdL, VdotN);
								if(blinn_f > 1e-5){
									kss.r = blinn_f * skinSpecColor.r;
									kss.g = blinn_f * skinSpecColor.g;
									kss.b = blinn_f * skinSpecColor.b;
								}
						}
						kss.r *= skin_spec_map.r;
						kss.g *= skin_spec_map.g;
						kss.b *= skin_spec_map.b;
						
						if(!*mi_eval_boolean(&p->directBlinnSpecular)){
								ks.r = (alpha_hair*khs.r + alpha_skin*kss.r)*Cl.r;
								ks.g = (alpha_hair*khs.g + alpha_skin*kss.g)*Cl.g;
								ks.b = (alpha_hair*khs.b + alpha_skin*kss.b)*Cl.b;
						}else {
											ks.r = (alpha_hair*khs.r + kss.r)*Cl.r;
											ks.g = (alpha_hair*khs.g + kss.g)*Cl.g;
											ks.b = (alpha_hair*khs.b + kss.b)*Cl.b;
									}
				}
				
				//diffuse
				miColor kd;
				kd.r = kd.g = kd.b = 0;
				
				if(emitDiffuse)
				{
						miColor khd;	//hair diffuse
						khd.r = khd.g = khd.b = 0;
						if(hairDiffMode==0)
						{
							khd.r = hairNdL * hairColor.r;
							khd.g = hairNdL * hairColor.g;
							khd.b = hairNdL * hairColor.b;
						}
						else {
								khd.r = sinTL * hairColor.r;
								khd.g = sinTL * hairColor.g;
								khd.b = sinTL * hairColor.b;
							}
						khd.r *= hairDiffMap.r;
						khd.g *= hairDiffMap.g;
						khd.b *= hairDiffMap.b;
						
						miColor ksd;	//skin diffuse
						ksd.r = ksd.g = ksd.b = 0;
						ksd.r = skinNdL * skinColor.r;
						ksd.g = skinNdL * skinColor.g;
						ksd.b = skinNdL * skinColor.b;
						ksd.r *= skinDiffMap.r;
						ksd.g *= skinDiffMap.g;
						ksd.b *= skinDiffMap.b;
						
						kd.r = (alpha_hair*khd.r + alpha_skin*ksd.r)*Cl.r;
						kd.g = (alpha_hair*khd.g + alpha_skin*ksd.g)*Cl.g;
						kd.b = (alpha_hair*khd.b + alpha_skin*ksd.b)*Cl.b;
				}

				 Cs.r += ks.r + kd.r;
				 Cs.g += ks.g + kd.g;
				 Cs.b += ks.b + kd.b;
			}
			samples = iter->get_number_of_samples();
			if (samples) {
				Cs.r /= (miScalar)samples;
				Cs.g /= (miScalar)samples;
				Cs.b /= (miScalar)samples;
				Ci.r += Cs.r;
				Ci.g += Cs.g;
				Ci.b += Cs.b;
			}
		}
	}
	
	// sss component
	miColor sssCol = *mi_eval_color( &p->sss );
	miColor sssMap = *mi_eval_color( &p->sssMapMult );

	Ci.r += sssCol.r * sssMap.r;
	Ci.g += sssCol.g * sssMap.g;
	Ci.b += sssCol.b * sssMap.b;

   // handle opacity
   *result = Ci;
   result->a = 1.0f;
   
   miColor incan = *mi_eval_color(&p->incandescence);
   result->r += incan.r;
   result->g += incan.g;
   result->b += incan.b;

   return miTRUE;
}

