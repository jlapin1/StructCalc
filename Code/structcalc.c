/*Whats new
- Different sigma values for each residue
	- sigma tensor is set before every calculation in plane_orientation_find, measurables, csdc
	- sigma values are read in from file sigin.csv, and printed out to file sigout.csv
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sort.h>
#include <dirent.h>

int
	kb_max   = 1,
	dry_runs = 1,
	CH       = 1,
	RAMA     = 1,
	stride   = 2,		// number of planes evaluated
	simits   = 1000,	// At least this many simplex iterations to find a low score
	stepbk   = 5,
	file     = 1;
double
	tol       = 1E3,	// Proceed to next residue only if tol is passed; do not set above 1E10
	tol2      = 1E0,	// To speed up algorithm, stop search once score is below tol2
    ramcut    = 2E3,
	noisiness = 0,
    wt[]      = {1,1,1},
	PHI0,PSI0,dPHI,dPSI,
	PHI_max,PHI_min,PSI_max,PSI_min,
	***surf_RAM,
	elapsed;
int
    sz,
	RAM_crit[5];
double
	pi = 3.141592653589793238,
	d2r = 0.01745329251994329547,  // rad deg^-1
	car = 50.699,
	a_nc = 151.8,
	a_n = 12.5,
	tetra = 110.5,
	tetid,
	HNCa = 118.2,
	NCCa = 115.6,
	HNCo = 119.5,
	gamma1 = 18.5,
	gamma2 = 90,
	chi0 = 10545.76120, // Hz
	chi1 = 23334.73464537244, // Hz
	omega,
	angle1,
	angle2,
	sixty;

char * confpath = "config.conf";

gsl_complex
	i;

gsl_matrix_complex
	*Dcs, *Dcsdag, *Mcs, *T_Q, *Tdag,
	*inner1, *inner2, *inner3, *inner1gly, *inner4, *inner1Q, *inner2Q, *inner3Q, *inner3PQ, *inner4Q, *inner1glyQ, *ytetra, *eang12;

gsl_rng
	* r;
	
const gsl_multimin_fminimizer_type * TT;

clock_t
	begin,
	end;

typedef struct{
	int
		* resnum,
		lower,
		upper,
		ind,
		stride;
	double
		*** params,
		** freq,
		** target,
		** sigma,
		answers[3],
		THETA0,
		PHI0;
	char
		* restyp;
	gsl_vector_complex
		* Q0;
} robj;
typedef struct{
	int 
		ind[3],
		typ[3];
	double
		ans[10],
		low[10],
		err;
	gsl_vector
		* guess;
} aobj;
void readconf(char *fpath){
	// This subroutine sets a bunch of global variables
	// atoll(): convert string to integer
	// atof(): convert string to float
	FILE * u = fopen(fpath, "r");
	char * nm = (char*)malloc(20*sizeof(char)),
		 * opt = (char*)malloc(20*sizeof(char)),
         buf[4][10];
	int fout, m, cnt1, cnt2;
	while(feof(u)==0){
		fout = fscanf(u,"%s", nm); // Read first string
		if(*nm!='#'){ // If not a comment line
			fout = fscanf(u,"%s\n",opt); // Read second string
			if (strcmp(nm, "kb_max")==0)
				kb_max = atoll(opt);
			else if(strcmp(nm,"CH")==0)
				CH = atoll(opt);
			else if(strcmp(nm,"RAMA")==0)
				RAMA = atoll(opt);
            else if(strcmp(nm,"ramcut")==0)
                ramcut = atof(opt);
			else if(strcmp(nm,"stride")==0)
				stride = atoll(opt);
			else if(strcmp(nm,"simits")==0)
				simits = atoll(opt);
			else if(strcmp(nm,"stepbk")==0)
				stepbk = atoll(opt);
			else if(strcmp(nm,"file")==0)
				file = atoll(opt);
			else if(strcmp(nm,"tol")==0)
				tol = atof(opt);
			else if(strcmp(nm,"tol2")==0)
				tol2 = atof(opt);
			else if(strcmp(nm,"noisiness")==0)
				noisiness = atof(opt);
			else if(strcmp(nm,"car")==0)
				car = atof(opt);
            else if(strcmp(nm,"wt")==0){
                cnt1=0;cnt2=0;
                for(m=0;m<strlen(opt);m++){
                    if(opt[m]==','){
                        buf[cnt1][cnt2] = '\0';
                        cnt1++;
                        if(cnt1==3)
                            perror("Error: Too many values for wt");
                        cnt2=0;
                    }
                    else if(opt[m]=='\n')
                        break;
                    else{
                        buf[cnt1][cnt2] = opt[m];
                        cnt2++;
                    }
                }
                for(m=0;m<3;m++)
                    wt[m] = atof(buf[m]);
            }
            else if(strcmp(nm,"ramrng")==0){
                cnt1=0;cnt2=0;
                for(m=0;m<strlen(opt);m++){
                    if(opt[m]==','){
                        buf[cnt1][cnt2] = '\0';
                        cnt1++;
                        if(cnt1==4)
                            perror("Error: Too many values for ramrng");
                        cnt2=0;
                    }
                    else if(opt[m]=='\n')
                        break;
                    else{
                        buf[cnt1][cnt2] = opt[m];
                        cnt2++;
                    }
                }
                PHI0 = atof(buf[0]);
                PSI0 = atof(buf[1]);
                dPHI = atof(buf[2]);
                dPSI = atof(buf[3]);
            }
			else if(strcmp(nm,"a_nc")==0)
				a_nc = atof(opt);
			else if(strcmp(nm,"a_n")==0)
				a_n = atof(opt);
            else if(strcmp(nm,"tetra")==0)
                tetra = atof(opt);
            else if(strcmp(nm,"HNCa")==0)
                HNCa = atof(opt);
            else if(strcmp(nm,"NCCa")==0)
                NCCa = atof(opt);
            else if(strcmp(nm,"HNCo")==0)
                HNCo = atof(opt);
            else if(strcmp(nm,"gamma1")==0)
                gamma1 = atof(opt);
            else if(strcmp(nm,"gamma2")==0)
                gamma2 = atof(opt);
            else if(strcmp(nm,"chi0")==0)
                chi0 = atof(opt);
            else if(strcmp(nm,"chi1")==0)
                chi1 = atof(opt);
		}
		else // If a comment (#), skip everything that is not a newline and then read the newline
			fout = fscanf(u,"%*[^\n]\n");
	}
	fclose(u);
	free(nm);free(opt);
}
void wigner(double one, double two, double three, gsl_matrix_complex * out){
	double
		cos_b,
		sin_b_sq_2;
	gsl_complex
		ex_ia,
		ex_ig,
		D33,D31,D23,D32;
	cos_b = cos(two);
	sin_b_sq_2 = sin(two)/sqrt(2);
	ex_ia = gsl_complex_exp(gsl_complex_mul_real(i,one));
	ex_ig = gsl_complex_exp(gsl_complex_mul_real(i,three));
	D33 = gsl_complex_mul(ex_ia,gsl_complex_mul_real(ex_ig,((1+cos_b)/2)));
	D31 = gsl_complex_mul(ex_ia,gsl_complex_mul_real(gsl_complex_conjugate(ex_ig),((1-cos_b)/2)));
	D23 = gsl_complex_mul_real(ex_ig,sin_b_sq_2);
	D32 = gsl_complex_mul_real(ex_ia,sin_b_sq_2);
	gsl_matrix_complex_set(out,0,0,gsl_complex_conjugate(D33));
	gsl_matrix_complex_set(out,0,1,gsl_complex_mul_real(gsl_complex_conjugate(D32),-1));
	gsl_matrix_complex_set(out,0,2,gsl_complex_conjugate(D31));
	gsl_matrix_complex_set(out,1,0,gsl_complex_conjugate(D23));
	gsl_matrix_complex_set(out,1,1,gsl_complex_rect(cos_b,0.0));
	gsl_matrix_complex_set(out,1,2,gsl_complex_mul_real(D23,-1));
	gsl_matrix_complex_set(out,2,0,D31);
	gsl_matrix_complex_set(out,2,1,D32);
	gsl_matrix_complex_set(out,2,2,D33);
}
void M_csa(double s1, double s2, double s3, gsl_matrix_complex * out){
	gsl_matrix_complex_set(out,0,0,gsl_complex_rect((s1+s2)/2,0.0));
	gsl_matrix_complex_set(out,0,1,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,0,2,gsl_complex_rect((s2-s1)/2,0.0));
	gsl_matrix_complex_set(out,1,0,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,1,1,gsl_complex_rect(s3,0.0));
	gsl_matrix_complex_set(out,1,2,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,2,0,gsl_complex_rect((s2-s1)/2,0.0));
	gsl_matrix_complex_set(out,2,1,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,2,2,gsl_complex_rect((s1+s2)/2,0.0));
}
void M_dc(double Chi, gsl_matrix_complex * out){
	gsl_matrix_complex_set(out,0,0,gsl_complex_rect(-0.5*Chi,0.0));
	gsl_matrix_complex_set(out,0,1,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,0,2,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,1,0,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,1,1,gsl_complex_rect(1.0*Chi,0.0));
	gsl_matrix_complex_set(out,1,2,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,2,0,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,2,1,gsl_complex_rect(0.0,0.0));
	gsl_matrix_complex_set(out,2,2,gsl_complex_rect(-0.5*Chi,0.0));
}
void YyY(double theta, double phi, gsl_vector_complex * Y){
	gsl_complex
		holdcn,
		holdcn2;
	holdcn  = gsl_complex_exp(gsl_complex_rect(0,phi));
	holdcn2 = gsl_complex_mul_real(holdcn,(-sin(theta)/sqrt(2.0)));
	gsl_vector_complex_set(Y,0,holdcn2);
	gsl_vector_complex_set(Y,1,gsl_complex_rect(cos(theta),0));
	holdcn = gsl_complex_exp(gsl_complex_rect(0,-phi));
	holdcn2 = gsl_complex_mul_real(holdcn,(sin(theta)/sqrt(2.0)));
	gsl_vector_complex_set(Y,2,holdcn2);
}
void mat_mult(gsl_matrix_complex * x, gsl_matrix_complex * y, gsl_matrix_complex * result){
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,x,y,GSL_COMPLEX_ZERO,result);
	//gsl_blas_zsymm(CblasLeft,CblasUpper,GSL_COMPLEX_ONE,x,y,GSL_COMPLEX_ZERO,result);
}
void dagger(gsl_matrix_complex * in, gsl_matrix_complex * out){
	int m,n;
	for(m=0;m<3;m++)
		for(n=0;n<3;n++)
				gsl_matrix_complex_set(out,m,n,gsl_complex_conjugate(gsl_matrix_complex_get(in,n,m)));
}
void daggerv(gsl_vector_complex * in, gsl_vector_complex * out){
	int m;
	for(m=0;m<3;m++)
		gsl_vector_complex_set(out,m,gsl_complex_conjugate(gsl_vector_complex_get(in,m)));
}
void matinv(gsl_matrix_complex * in, gsl_matrix_complex * out){
	int m;
	gsl_matrix_complex * copy = gsl_matrix_complex_alloc(3,3);
	gsl_permutation * p = gsl_permutation_alloc(3);
	gsl_matrix_complex_memcpy(copy,in);
	gsl_linalg_complex_LU_decomp(copy, p, &m);
	gsl_linalg_complex_LU_invert(copy, p, out);
	gsl_matrix_complex_free(copy);
	gsl_permutation_free(p);
}
void csatensor(double s1, double s2, double s3, int Q, gsl_matrix_complex * out){
	int m,n;
	gsl_matrix_complex 
		* Mcs =  gsl_matrix_complex_alloc(3,3),
		* holdcm = gsl_matrix_complex_alloc(3,3);
	
	M_csa(s1, s2, s3, Mcs);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,Dcs,Mcs,GSL_COMPLEX_ZERO,holdcm);// inner1
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,holdcm,Dcsdag,GSL_COMPLEX_ZERO,out);// inner1
	
	for(m=0;m<3;m++)
		for(n=0;n<3;n++)
			gsl_matrix_complex_set(out,m,n,gsl_complex_mul_real(gsl_matrix_complex_get(out,m,n),car));
	
	if(Q==1){
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,out,Tdag,GSL_COMPLEX_ZERO,holdcm);
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,T_Q,holdcm,GSL_COMPLEX_ZERO,out);
	}	
	gsl_matrix_complex_free(Mcs);gsl_matrix_complex_free(holdcm);
}
void csdc(robj * robo, gsl_vector_complex * startY, double ** parms, double ** freqs){
	int m,n;
	double hn, hp;
	gsl_matrix_complex
		* hold1 = gsl_matrix_complex_alloc(3,3),
		* hold2 = gsl_matrix_complex_alloc(3,3),
		* hold3 = gsl_matrix_complex_alloc(3,3);
	gsl_vector_complex
		* ypre = gsl_vector_complex_alloc(3),
		* ypost = gsl_vector_complex_alloc(3),
		* ych = gsl_vector_complex_alloc(3),
		* ydag = gsl_vector_complex_alloc(3),
		* holdv = gsl_vector_complex_alloc(3);
	gsl_complex
		hold;
	
	// First resonances
	daggerv(startY,ydag);
	csatensor(robo->sigma[0][0], robo->sigma[0][1], robo->sigma[0][2], 0, inner1);
	gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,inner1,ydag,GSL_COMPLEX_ZERO,holdv);
	gsl_blas_zdotu(startY,holdv,&hold);
	freqs[0][0] = GSL_REAL(hold)/*/car*/;
	gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,inner2,ydag,GSL_COMPLEX_ZERO,holdv);
	gsl_blas_zdotu(startY,holdv,&hold);
	freqs[0][1] = fabs(GSL_REAL(hold));
	
	// Main loop
	for(n=0;n<3;n++)
			gsl_vector_complex_set(ypre,n,gsl_vector_complex_get(startY,n));
	for(m=0;m<sz-1;m++){
		if(CH==1){
			wigner(0, -pi/2, 0, hold2);
			hp = 0;
			if(robo->restyp[m]=='G'){
				wigner(angle1, parms[m][0]+sixty, (pi/2)-tetid, hold1);
				gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,hold1,hold2,GSL_COMPLEX_ZERO,hold3);
				gsl_blas_zgemv(CblasTrans,GSL_COMPLEX_ONE,hold3,ypre,GSL_COMPLEX_ZERO,ych);
				hp = chi1*fabs((((3*pow(GSL_REAL(gsl_vector_complex_get(ych,1)),2.0))-1)/2.0));
			}
			wigner(angle1, parms[m][0]-sixty, (pi/2)-tetid, hold1);
			gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,hold1,hold2,GSL_COMPLEX_ZERO,hold3);
			gsl_blas_zgemv(CblasTrans,GSL_COMPLEX_ONE,hold3,ypre,GSL_COMPLEX_ZERO,ych);
			hn = chi1*fabs((((3*pow(GSL_REAL(gsl_vector_complex_get(ych,1)),2.0))-1)/2.0));
			freqs[m][2] = sqrt( (hp*hp) + (hn*hn) );
		}
		wigner(angle1,parms[m][0],tetra,hold1);
		wigner(0,(-1*parms[m][1])-pi,angle2,hold2);
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,hold1,hold2,GSL_COMPLEX_ZERO,hold3);
		gsl_blas_zgemv(CblasTrans,GSL_COMPLEX_ONE,hold3,ypre,GSL_COMPLEX_ZERO,ypost);
		for(n=0;n<3;n++)
			gsl_vector_complex_set(ypre,n,gsl_vector_complex_get(ypost,n));

		daggerv(ypost,ydag);
		csatensor(robo->sigma[m+1][0], robo->sigma[m+1][1], robo->sigma[m+1][2], 0, inner1);
		gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,inner1,ydag,GSL_COMPLEX_ZERO,holdv);
		gsl_blas_zdotu(ypost,holdv,&hold);
		freqs[m+1][0] = GSL_REAL(hold)/*car*/;
		gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,inner2,ydag,GSL_COMPLEX_ZERO,holdv);
		gsl_blas_zdotu(ypost,holdv,&hold);
		freqs[m+1][1] = fabs(GSL_REAL(hold));
        
	}
	gsl_matrix_complex_free(hold1);gsl_matrix_complex_free(hold2);gsl_matrix_complex_free(hold3);gsl_vector_complex_free(ypre);
	gsl_vector_complex_free(ypost);gsl_vector_complex_free(ydag);gsl_vector_complex_free(ych);gsl_vector_complex_free(holdv);
}
void simplex(robj * robo, aobj * aob, int num, double (*func)(const gsl_vector * v, void * robo)){
	gsl_multimin_fminimizer * s = NULL;
	gsl_vector 
		*ss, // Updated answer
		*xx; // Initial guess
	gsl_multimin_function minex_func;
	size_t iter = 0;
	int m,status;
	double size;
	
	xx = gsl_vector_alloc(num);
	ss = gsl_vector_alloc(num);//holds the answer
	gsl_vector_set_all(ss,1.0);

	minex_func.n = num;
	minex_func.f = func;
	minex_func.params = robo;
	
	s = gsl_multimin_fminimizer_alloc(TT,num);
	gsl_multimin_fminimizer_set(s,&minex_func,aob->guess,ss);
	
	do{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s); // Take a step
	
		if(status)
			break;
	
		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size,1E-7);
	
		/*if(status==GSL_SUCCESS)
			printf("converged to minimum at\n",(gsl_vector_get(s->x,0)));
	
		printf("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",iter,gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),s->fval,size);*/
	}
	while(status==GSL_CONTINUE && iter<1000);
	
	for(m=0;m<num;m++)
		aob->ans[m] = gsl_vector_get(s->x,m);
	aob->err = s->fval;
	
	gsl_vector_free(xx);gsl_vector_free(ss);gsl_multimin_fminimizer_free(s);
}
double plane_orientation_find(const gsl_vector * v, void * rob){
	robj * robo = (robj*)rob;
	int I = robo->ind;
	double e1,e2,etot,
		PHI0 = gsl_vector_get(v,0),
		THETA0   = gsl_vector_get(v,1),
		cs, dc, cn;
	
	if( (PHI0>pi) || (PHI0<0) || (THETA0>(2*pi)) || (THETA0<0) )
		return 1E6;
	
	gsl_complex holdcn, holdcn2;
	gsl_vector_complex
		* Y0     = gsl_vector_complex_alloc(3),
		* Y0dag  = gsl_vector_complex_alloc(3),
		* holdcv = gsl_vector_complex_alloc(3);
	
	YyY(THETA0,PHI0,Y0);
	daggerv(Y0,Y0dag);
	gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,inner1,Y0dag,GSL_COMPLEX_ZERO,holdcv);
	gsl_blas_zdotu(Y0,holdcv,&holdcn);
	cs = GSL_REAL(holdcn);
	gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,inner2,Y0dag,GSL_COMPLEX_ZERO,holdcv);
	gsl_blas_zdotu(Y0,holdcv,&holdcn);
	dc = fabs(GSL_REAL(holdcn));
	
	e1 = (cs-robo->target[0][0]);
	e2 = dc-robo->target[0][1];
	etot = sqrt( (e1*e1) + (e2*e2) );
	
	gsl_vector_complex_free(Y0);gsl_vector_complex_free(Y0dag);gsl_vector_complex_free(holdcv);
	
	if( (PHI0>2*pi) || (PHI0<0) || (THETA0>pi) || (THETA0<0) )
		etot=1E6;
	return etot;
}
void measurables(robj * R, gsl_vector_complex * Qin, double PHI, double PSI, gsl_vector_complex * Qout, double * res){
	double hn, hp;
	gsl_complex
		cn0, cn1, cn2, cn3, cn4, cn5, cn6,
		Qi0, Qi1, Qi2;
	gsl_vector_complex
		* Q_i = gsl_vector_complex_alloc(3),
		* Q_i2 = gsl_vector_complex_alloc(3);
	
	int n = R->ind;
	
	cn0 = gsl_vector_complex_get(Qin,0);
	Qi0 = gsl_complex_mul(gsl_complex_exp(gsl_complex_rect(0.0,PHI)),cn0);
	Qi1 = gsl_vector_complex_get(Qin,1);
	cn0 = gsl_complex_conjugate(Qi0);// for some reason conjugate is taken of Qi1, NOT Qi0 as intended
	Qi2 = gsl_complex_mul_real(cn0,-1.0);
	gsl_vector_complex_set(Q_i,0,Qi0);gsl_vector_complex_set(Q_i,1,Qi1);gsl_vector_complex_set(Q_i,2,Qi2);
	
	if(CH==1){
		// CALC3
		hp = 0;
		if(R->restyp[n]=='G'){
			cn0 = gsl_matrix_complex_get(inner3PQ,0,0);cn1 = gsl_matrix_complex_get(inner3PQ,1,1);
			cn2 = gsl_matrix_complex_get(inner3PQ,0,1);cn3 = gsl_matrix_complex_get(inner3PQ,0,2);
			// Middle term
			cn4 = gsl_complex_add(gsl_complex_mul(gsl_complex_sub(cn1,cn0),Qi1),gsl_complex_mul_real(gsl_complex_mul(cn2,Qi0),4.0));
			cn5 = gsl_complex_mul(Qi1,cn4);
			// Last term
			cn6 = gsl_complex_mul_real(gsl_complex_mul(gsl_complex_pow_real(Qi0,2.0),cn3),2.0);
			
			hp = fabs(GSL_REAL(gsl_complex_sub(gsl_complex_add(cn0,cn5),cn6)));
		}
		cn0 = gsl_matrix_complex_get(inner3Q,0,0);cn1 = gsl_matrix_complex_get(inner3Q,1,1);
		cn2 = gsl_matrix_complex_get(inner3Q,0,1);cn3 = gsl_matrix_complex_get(inner3Q,0,2);
		// Middle term
		cn4 = gsl_complex_add(gsl_complex_mul(gsl_complex_sub(cn1,cn0),Qi1),gsl_complex_mul_real(gsl_complex_mul(cn2,Qi0),4.0));
		cn5 = gsl_complex_mul(Qi1,cn4);
		// Last term
		cn6 = gsl_complex_mul_real(gsl_complex_mul(gsl_complex_pow_real(Qi0,2.0),cn3),2.0);
		
		hn = fabs(GSL_REAL(gsl_complex_sub(gsl_complex_add(cn0,cn5),cn6)));
		
		res[2] = sqrt( (hn*hn) + (hp*hp) );
	}
	
	gsl_blas_zgemv(CblasTrans,GSL_COMPLEX_ONE,ytetra,Q_i,GSL_COMPLEX_ZERO,Q_i2);
	
	Qi0 = gsl_complex_mul( gsl_complex_exp(gsl_complex_rect(0.0,-1*(pi+PSI))) , gsl_vector_complex_get(Q_i2,0) );
	Qi1 = gsl_vector_complex_get(Q_i2,1);
	cn0 = gsl_complex_conjugate(Qi0);
	Qi2 = gsl_complex_mul_real(cn0,-1.0);
	
	// CALC1
	csatensor(R->sigma[n+1][0], R->sigma[n+1][1], R->sigma[n+1][2], 1, inner1Q);
	cn0 = gsl_matrix_complex_get(inner1Q,0,0);cn1 = gsl_matrix_complex_get(inner1Q,1,1);
	cn2 = gsl_matrix_complex_get(inner1Q,0,1);cn3 = gsl_matrix_complex_get(inner1Q,0,2);
	
	// Middle term
	cn4 = gsl_complex_add(gsl_complex_mul(gsl_complex_sub(cn1,cn0),Qi1),gsl_complex_mul_real(gsl_complex_mul(cn2,Qi0),4.0));
	cn5 = gsl_complex_mul(Qi1,cn4);
	// Last term
	cn6 = gsl_complex_mul_real(gsl_complex_mul(gsl_complex_pow_real(Qi0,2.0),cn3),2.0);
	
	res[0] = GSL_REAL(gsl_complex_sub(gsl_complex_add(cn0,cn5),cn6));
	
	// CALC2
	cn0 = gsl_matrix_complex_get(inner2Q,0,0);cn1 = gsl_matrix_complex_get(inner2Q,1,1);
	cn2 = gsl_matrix_complex_get(inner2Q,0,1);cn3 = gsl_matrix_complex_get(inner2Q,0,2);
	
	// Middle term
	cn4 = gsl_complex_add(gsl_complex_mul(gsl_complex_sub(cn1,cn0),Qi1),gsl_complex_mul_real(gsl_complex_mul(cn2,Qi0),4.0));
	cn5 = gsl_complex_mul(Qi1,cn4);
	// Last term
	cn6 = gsl_complex_mul_real(gsl_complex_mul(gsl_complex_pow_real(Qi0,2.0),cn3),2.0);
	
	res[1] = fabs(GSL_REAL(gsl_complex_sub(gsl_complex_add(cn0,cn5),cn6)));
	
	gsl_vector_complex_set(Q_i,0,Qi0);gsl_vector_complex_set(Q_i,1,Qi1);gsl_vector_complex_set(Q_i,2,Qi2);
	gsl_blas_zgemv(CblasTrans,GSL_COMPLEX_ONE,eang12,Q_i,GSL_COMPLEX_ZERO,Qout);
	
	gsl_vector_complex_free(Q_i);gsl_vector_complex_free(Q_i2);
}
double sigmaQ2(const gsl_vector * v, void * rob){
	robj * robo = (robj*)rob;
	int
		m,o,I,J,
		n = robo->ind;
	double
		PHI,PSI,
		e1,e2,e3,etot,E[3];
	
	if((RAMA==1)){
		for(m=0;m<robo->stride;m++){
			PHI = v->data[2*m];
			PSI = v->data[(2*m)+1];
			if(PHI>pi)
				PHI = -pi + fmod(PHI,pi);
			else if(PHI<-pi)
				PHI = pi + fmod(PHI,pi);
			if(PSI>pi)
				PSI = -pi + fmod(PSI,pi);
			else if(PSI<-pi)
				PSI = pi + fmod(PSI,pi);
			I = round((PHI+pi-0.01745329)/0.034906585);
			if(I==180)
				I-=1;
			J = round((PSI+pi-0.01745329)/0.034906585);
			if(J==180)
				J-=1;
			if(surf_RAM[RAM_crit[m]][I][J]<ramcut) // 5.05E-3
				return 1E5;
		}
	}
	gsl_vector_complex
		* Q_i = gsl_vector_complex_alloc(3),
		* Q_i2 = gsl_vector_complex_alloc(3);
	
	for(m=0;m<3;m++)
		gsl_vector_complex_set(Q_i, m, gsl_vector_complex_get(robo->Q0, m));
	etot=0;
	for(o=0;o<robo->stride;o++){
		measurables(robo, Q_i, v->data[o*2], v->data[(o*2)+1], Q_i2, E);
		e1 = (E[0]-robo->target[robo->ind+1][0])*wt[0];
		e2 = (E[1]-robo->target[robo->ind+1][1])*wt[1];
		etot += (e1*e1) + (e2*e2);
		if(CH==1){
			e3 = (E[2]-robo->target[robo->ind][2])*wt[2];
			etot += e3*e3;
		}
		if(o<robo->stride-1)
			for(m=0;m<3;m++)
				gsl_vector_complex_set(Q_i, m, gsl_vector_complex_get(Q_i2, m));
		robo->ind++;
	}
	robo->ind=n;
	gsl_vector_complex_free(Q_i);gsl_vector_complex_free(Q_i2);
	
	return sqrt(etot);
}
void noise(robj * robo, double radlim){
	int m;
	double radius, ang, ang2, dx, dy, dz;
	
	for(m=0;m<sz;m++){
		robo->target[m][0] += gsl_ran_flat(r, -radlim, radlim);
		robo->target[m][1] += gsl_ran_flat(r, -radlim, radlim);
		if(CH==1)
			robo->target[m][2] += gsl_ran_flat(r, -radlim, radlim);
	}
}
void printdv(int length, double * A){
	int x;
	printf("\n");	
	for(x=0;x<length;x++){
		printf("%.5lf ", A[x]);
		printf("\n");
	}
}
void printda(int rows, int cols, double ** A){
	int x,y;
	printf("\n");
	for(x=0;x<rows;x++){
		for(y=0;y<cols;y++)
			printf("%10.5lf   ",A[x][y]);
		printf("\n");
	}
}
void printcv(int length, gsl_vector_complex * A){
	int m;
	gsl_complex hold;
	for(m=0;m<length;m++){
		hold = gsl_vector_complex_get(A,m);
		printf("%lf + %lfi\n",GSL_REAL(hold),GSL_IMAG(hold));
	}
}
void printcm(int rows, int columns, gsl_matrix_complex * A){
	int m,n;
	printf("\n");
	for(m=0;m<rows;m++){
		for(n=0;n<columns;n++)
			printf("%10.5lf + %10.5lfi ",GSL_REAL(gsl_matrix_complex_get(A,m,n)),GSL_IMAG(gsl_matrix_complex_get(A,m,n)));
		printf("\n");
	}
	printf("\n");
}
void copydv(int length, double * one, double * two){
	int m;

	// two into one
	for(m=0;m<length;m++)
		one[m] = two[m];
}
int main(){
	robj
		real,
		ran,
		sys;
	
	aobj
		A;
	
	int m,n,o,p;
	double temp[3];
    
    // Look for configuration script "config.conf"
    struct dirent * de;
    DIR * dr = opendir("./");
    if(dr==NULL){
        printf("\nCould not open directory.\n");
        return 0;
    }
    int filexist = 0;
    while((de = readdir(dr)) != NULL){
        filexist = strcmp(de->d_name, "config.conf");
        if(filexist==0)
            break;
    }
    closedir(dr);
    
    // Read into configuration script to assign global variables
    if(filexist==1){
        printf("\nNo configuration file \"config.conf\" found in current directory\n\n");
        exit(0);
    }
    else
        readconf(confpath);
	
	// Constants and Angles
	i=gsl_complex_rect(0.0,1.0);
	omega = pi;
	a_nc *= d2r;NCCa *= d2r;a_n *= d2r;tetra *= d2r;HNCa *= d2r;HNCo *= d2r;gamma1 *= d2r;gamma2 *= d2r;
	tetid = tetra;//acos(-1.0/3.0);
	angle1 = ((3*pi)/2)-HNCa;
	angle2 = ((3*pi)/2)-NCCa-HNCo;
	sixty = 60*pi/180;//acos((-tan(tetra/2.0))/tan(tetid));
	
	// Random Numbers
	srand(time(NULL));
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r,rand());
	TT = gsl_multimin_fminimizer_nmsimplex2;
	
	// Allocation
	gsl_complex
		holdcn,
		holdcn2,
		holdcn3;
	
	gsl_vector_complex
		* y = gsl_vector_complex_alloc(3),
		* holdcv = gsl_vector_complex_alloc(3),
		* holdcv2 = gsl_vector_complex_alloc(3);
	
	gsl_matrix_complex
		* Mcsgly = gsl_matrix_complex_alloc(3,3),* Mdc = gsl_matrix_complex_alloc(3,3),
		* Ddc = gsl_matrix_complex_alloc(3,3),* Ddcinv = gsl_matrix_complex_alloc(3,3),
		* holdcm = gsl_matrix_complex_alloc(3,3),* holdcm2 = gsl_matrix_complex_alloc(3,3),* holdcm3 = gsl_matrix_complex_alloc(3,3),
		* holdcm4 = gsl_matrix_complex_alloc(3,3),* holdcm5 = gsl_matrix_complex_alloc(3,3),* holdcm6 = gsl_matrix_complex_alloc(3,3);
	Mcs = gsl_matrix_complex_alloc(3,3);Dcs = gsl_matrix_complex_alloc(3,3);Dcsdag = gsl_matrix_complex_alloc(3,3);T_Q = gsl_matrix_complex_alloc(3,3);Tdag = gsl_matrix_complex_alloc(3,3);
	inner1 = gsl_matrix_complex_alloc(3,3);inner1gly = gsl_matrix_complex_alloc(3,3);inner2 = gsl_matrix_complex_calloc(3,3);inner3 = gsl_matrix_complex_alloc(3,3);inner4 = gsl_matrix_complex_alloc(3,3);
	inner1Q = gsl_matrix_complex_alloc(3,3);inner1glyQ = gsl_matrix_complex_alloc(3,3);inner2Q = gsl_matrix_complex_alloc(3,3);inner3Q = gsl_matrix_complex_alloc(3,3);inner3PQ = gsl_matrix_complex_alloc(3,3);inner4Q = gsl_matrix_complex_alloc(3,3);
	ytetra = gsl_matrix_complex_alloc(3,3);real.Q0 = gsl_vector_complex_alloc(3);sys.Q0 = gsl_vector_complex_alloc(3);
	
	// Constant matrices
	wigner(gamma1, pi/2, gamma2, Dcs);// 3 angles are omN in python
	dagger(Dcs, Dcsdag);// inner1
	gsl_matrix_complex_set(inner2,0,0,gsl_complex_rect(0.25,0));// inner2
	gsl_matrix_complex_set(inner2,0,2,gsl_complex_rect(-0.75,0));// inner2
	gsl_matrix_complex_set(inner2,1,1,gsl_complex_rect(-0.5,0));// inner2
	gsl_matrix_complex_set(inner2,2,0,gsl_complex_rect(-0.75,0));// inner2
	gsl_matrix_complex_set(inner2,2,2,gsl_complex_rect(0.25,0));// inner2
	gsl_matrix_complex_memcpy(inner3,inner2);// inner3
	wigner(-HNCa,0,0,holdcm);dagger(holdcm,holdcm2);// inner4
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,inner2,holdcm2,GSL_COMPLEX_ZERO,holdcm3);// inner4
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,holdcm,holdcm3,GSL_COMPLEX_ZERO,inner4);// inner4
	for(m=0;m<3;m++)
		for(n=0;n<3;n++){
			gsl_matrix_complex_set(inner2,m,n,gsl_complex_mul_real(gsl_matrix_complex_get(inner2,m,n),chi0));
			gsl_matrix_complex_set(inner3,m,n,gsl_complex_mul_real(gsl_matrix_complex_get(inner3,m,n),chi1));
		}
	
	// Q basis
	gsl_matrix_complex
		* eIx     = gsl_matrix_complex_alloc(3,3),* eang2 = gsl_matrix_complex_alloc(3,3),* eang1 = gsl_matrix_complex_alloc(3,3),* e_Ix = gsl_matrix_complex_alloc(3,3),
		* D_omega = gsl_matrix_complex_alloc(3,3);
	eang12  = gsl_matrix_complex_alloc(3,3);
	double
		// expm(i*(pi/2)*Ix)
		re1[][3] = {
			{0.5 , 0, -0.5},
			{0   , 0,    0},
			{-0.5, 0,  0.5}},
		im1[][3] = {
			{0          , 1.0/sqrt(2),           0},
			{1.0/sqrt(2), 0          , 1.0/sqrt(2)},
			{0          , 1.0/sqrt(2),           0}},
		// expm(-i*tetra*Ly): tetra=110.5
		re3[][3] = {
			{0.324896309370267,  -0.662327256766391,   0.675103690629734},
			{0.662327256766391,  -0.350207381259467,  -0.662327256766391},
			{0.675103690629734,   0.662327256766391,   0.324896309370266}},
		// expm(-i(angle1+angle2)*Iy)
		re4[][3] = {
			{3.414675230757E-3,	8.249866936897E-2,	9.965853247692E-1},
			{-8.249866936897E-2,-9.931706495385E-1, 8.249866936897E-2},
			{9.965853247692E-1, -8.249866936897E-2, 3.414675230757E-3}};
	
	wigner(-angle2-NCCa, 0*d2r, (3*pi/2)-HNCo, D_omega); // D_omega: Middle value can vary
	// expm(-i*angle1*Iz); angle1=2.64940980452739
	gsl_matrix_complex_set(eang1,0,0,gsl_complex_rect(-0.881303452064992,-0.472550764869054));
	gsl_matrix_complex_set(eang1,1,1,GSL_COMPLEX_ONE);
	gsl_matrix_complex_set(eang1,2,2,gsl_complex_rect(-0.881303452064992,0.472550764869054));
	// expm(-i*angle2*Iz): angle2=0.609119908946021
	gsl_matrix_complex_set(eang2,0,0,gsl_complex_rect(0.820151875873772,-0.572145873445516));
	gsl_matrix_complex_set(eang2,1,1,GSL_COMPLEX_ONE);
	gsl_matrix_complex_set(eang2,2,2,gsl_complex_rect(0.820151875873772,0.572145873445516));
	for(m=0;m<3;m++){
		for(n=0;n<3;n++){
			gsl_matrix_complex_set(eIx,m,n,gsl_complex_rect(re1[m][n],im1[m][n]));
			gsl_matrix_complex_set(e_Ix,m,n,gsl_complex_rect(re1[m][n],-1*im1[m][n]));
			gsl_matrix_complex_set(ytetra,m,n,gsl_complex_rect(re3[m][n],0.0));
			gsl_matrix_complex_set(eang12,m,n,gsl_complex_rect(re4[m][n],0.0));
		}
	}
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,eIx,eang2,GSL_COMPLEX_ZERO,holdcm);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,holdcm,D_omega,GSL_COMPLEX_ZERO,T_Q);
	dagger(T_Q,Tdag);
	// inner2Q
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,inner2,Tdag,GSL_COMPLEX_ZERO,holdcm);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,T_Q,holdcm,GSL_COMPLEX_ZERO,inner2Q);
	// inner3Q
	wigner(pi/3.0,(pi/2.0)-tetra,0,holdcm);dagger(holdcm,holdcm5);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,e_Ix,holdcm5,GSL_COMPLEX_ZERO,holdcm6);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,inner3,holdcm6,GSL_COMPLEX_ZERO,holdcm5);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,eIx,holdcm5,GSL_COMPLEX_ZERO,holdcm6);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,holdcm,holdcm6,GSL_COMPLEX_ZERO,inner3Q);// inner3Q
	// inner3PQ
	wigner(-pi/3.0,(pi/2.0)-tetra,0,holdcm);dagger(holdcm,holdcm5);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,e_Ix,holdcm5,GSL_COMPLEX_ZERO,holdcm6);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,inner3,holdcm6,GSL_COMPLEX_ZERO,holdcm5);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,eIx,holdcm5,GSL_COMPLEX_ZERO,holdcm6);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,holdcm,holdcm6,GSL_COMPLEX_ZERO,inner3PQ);
	// inner4Q
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,inner4,Tdag,GSL_COMPLEX_ZERO,holdcm);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,T_Q,holdcm,GSL_COMPLEX_ZERO,inner4Q);
	
	
	// Read in ramachandran plot scores
	surf_RAM = (double***)malloc(4*sizeof(double**));
	// General: RAM_crit = 0
	FILE * f  = fopen("/home/nevlab/structcalc/prgms/RAMA/rama500general.data","r");
	surf_RAM[0] = (double**)malloc(180*sizeof(double*));
	for(m=0;m<180;m++){
		surf_RAM[0][m] = (double*)malloc(180*sizeof(double));
		for(n=0;n<180;n++)
			p = fscanf(f,"%*f %*f %lf\n",&surf_RAM[0][m][n]);
	}
	fclose(f);
	// Proline: RAM_crit = 1
	f  = fopen("/home/nevlab/structcalc/prgms/RAMA/rama500pro.data","r");
	surf_RAM[1] = (double**)malloc(180*sizeof(double*));
	for(m=0;m<180;m++){
		surf_RAM[1][m] = (double*)malloc(180*sizeof(double));
		for(n=0;n<180;n++)
			p = fscanf(f,"%*f %*f %lf\n",&surf_RAM[1][m][n]);
	}
	fclose(f);
	// Pre-Proline: RAM_crit = 2
	f  = fopen("/home/nevlab/structcalc/prgms/RAMA/rama500prepro.data","r");
	surf_RAM[2] = (double**)malloc(180*sizeof(double*));
	for(m=0;m<180;m++){
		surf_RAM[2][m] = (double*)malloc(180*sizeof(double));
		for(n=0;n<180;n++)
			p = fscanf(f,"%*f %*f %lf\n",&surf_RAM[2][m][n]);
	}
	fclose(f);
	//Post-Glycine: RAM_crit = 3
	f  = fopen("/home/nevlab/structcalc/prgms/RAMA/rama500gly_sym.data","r");
	surf_RAM[3] = (double**)malloc(180*sizeof(double*));
	for(m=0;m<180;m++){
		surf_RAM[3][m] = (double*)malloc(180*sizeof(double));
		for(n=0;n<180;n++)
			p = fscanf(f,"%*f %*f %lf\n",&surf_RAM[3][m][n]);
	}
	fclose(f);
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Setting up system
	
    // Read in targets
    f = fopen("input/realtarg.csv","r");
    p = fscanf(f,"%*[^\n]\n");
    p = fscanf(f,"%d\n",&sz);
    real.restyp = (char*)malloc((sz+1)*sizeof(char));sys.restyp = (char*)malloc((sz+1)*sizeof(char));
    real.freq = (double**)malloc(sz*sizeof(double*));
    sys.freq = (double**)malloc(sz*sizeof(double*));
    sys.target = (double**)malloc(sz*sizeof(double*));
    real.restyp[sz] = 'X';sys.restyp[sz] = 'X';
    for(m=0;m<sz;m++){
        real.freq[m] = (double*)malloc(3*sizeof(double));
        sys.freq[m] = (double*)malloc(3*sizeof(double));
        sys.target[m] = (double*)malloc(3*sizeof(double));
        p = fscanf(f,"%c,%lf,%lf,%lf\n", &real.restyp[m], &real.freq[m][0], &real.freq[m][1], &real.freq[m][2]);
        sys.restyp[m] = real.restyp[m];
    }
    fclose(f);
    
    // Read in sigmas
    sys.sigma = (double**)malloc(sz*sizeof(double*));
    f = fopen("input/sigin.csv", "r");
    p = fscanf(f,"%*[^\n]\n");
    for(m=0;m<sz;m++){
        sys.sigma[m] = (double*)malloc(3*sizeof(double));
        p = fscanf(f, "%lf,%lf,%lf\n", &sys.sigma[m][0], &sys.sigma[m][1], &sys.sigma[m][2]);
    }
    fclose(f);
	
	// Start constructing sys
	sys.params = (double***)malloc(kb_max*sizeof(double**));
	for(m=0;m<kb_max;m++){
		sys.params[m] = (double**)malloc((sz-1)*sizeof(double*));
		for(n=0;n<sz-1;n++)
			sys.params[m][n] = (double*)malloc(2*sizeof(double));
	}
	
	int **stepbkV  = (int**)malloc(2*sizeof(int*));stepbkV[0] = (int*)calloc(sz-1,sizeof(int));stepbkV[1] = (int*)calloc(sz-1,sizeof(int));
	double 
		* tolerance_all = (double*)malloc(sz*sizeof(double));
	
	A.guess = gsl_vector_alloc(2);
	int 
		S, kb=0, cnt, trials, nvar;
	double
		z_prev, theta00, phi00, z0, tolerance,
		dummy[3],
		** RMS = (double**)malloc(kb_max*sizeof(double*)),
		** rmsdf = (double**)malloc(sz*sizeof(double)),
		** freq_calc = (double**)malloc((sz-1)*sizeof(double*));
	gsl_vector_complex
		* Y00, * Y0, * Q00, * Q0;
	gsl_matrix_complex
		* Q_all = gsl_matrix_complex_alloc(sz,3);
	Y00 = gsl_vector_complex_alloc(3); Y0  = gsl_vector_complex_alloc(3);
	Q00 = gsl_vector_complex_alloc(3); Q0  = gsl_vector_complex_alloc(3);
	for(m=0;m<kb_max;m++){
        RMS[m] = (double*)calloc(4,sizeof(double));
		for(n=0;n<sz-1;n++)
			if(m==0){
				tolerance_all[n] = tol;
				rmsdf[n]         = (double*)malloc(kb_max*sizeof(double));
				freq_calc[n]     = (double*)calloc(4,sizeof(double));
			}
    }
	tolerance_all[sz-1] = tol;/*RAM_crit[sz-1]=5;*/rmsdf[sz-1] = (double*)malloc(kb_max*sizeof(double));freq_calc[sz-1] = (double*)calloc(4,sizeof(double));
	
	
	// Angle ranges
	PHI_max = (PHI0+dPHI)*d2r;PHI_min = (PHI0-dPHI)*d2r;PSI_max = (PSI0+dPSI)*d2r;PSI_min = (PSI0-dPSI)*d2r;
	
	// output.csv, freqs_RMS.csv
	FILE *fp, *fp2;
	if(file==1){
		fp = fopen("output/output.csv","w");
		fprintf(fp,"#Parameters: CH=%d; RAMA=%d; stride=%d; simits=%d; stepbk=%d; tol=%.0lf; tol2=%.0lf; noisiness=%.0lf; PHI0=%.0f PSI0=%.0f dPHI=%.0f dPSI=%.0f\n",CH,RAMA,stride,simits,stepbk,tol,tol2,noisiness,PHI0,PSI0,dPHI,dPSI);
		fprintf(fp,"%d,%d\n",sz,kb_max);
        fprintf(fp,"%s\n",sys.restyp);
		/*for(n=0;n<sz;n++)
			fprintf(fp,"%c,",sys.restyp[n]);
		fprintf(fp,"%c\n",sys.restyp[sz]);*/
		fclose(fp);
		fp2 = fopen("output/freqs_RMS.csv","w");
		fclose(fp2);
	}
	
	begin = clock();
	for(S=0;S<kb_max;S++){	
		
        // Create targets with uncertainty
        // First hand over real frequencies to system targets
        // -> reset to original targets at beginning of each run
		for(m=0;m<sz;m++)
			for(n=0;n<3;n++)
				sys.target[m][n] = real.freq[m][n];
        noise(&sys, noisiness);
		
		// Reset quantities that may have been altered during a run
		sys.stride = stride;
		nvar = 2*sys.stride;
		for(o=0;o<sz-1;o++){stepbkV[0][o] = 0;stepbkV[1][o] = 0;tolerance_all[o] = tol;}
		
		cnt=0;
		n=0;
		while(n<sz-1){
			// If not the last residue, set nvar based on original stride
			if((sz-1-n)<sys.stride){
				sys.stride = sz-1-n;
				nvar = 2*sys.stride;
				gsl_vector_free(A.guess);
				A.guess = gsl_vector_alloc(nvar);}
			
			// Set the phi/psi torsion randomization range for residue n
			for(o=0;o<sys.stride;o++){
				if(sys.restyp[n+o+1]=='P')
					RAM_crit[o] = 1;
				else if(sys.restyp[n+o+2]=='P')
					RAM_crit[o] = 2;
				else if(sys.restyp[n+o]=='G')
					RAM_crit[o] = 3;
				else
					RAM_crit[o] = 0;}
			
			printf("\nn: %d | kb: %d\n",n,S);
			sys.ind = n;
			if(n==0){
				csatensor(sys.sigma[0][0], sys.sigma[0][1], sys.sigma[0][2], 0, inner1);
				gsl_vector_free(A.guess);
				A.guess = gsl_vector_alloc(2);
				z_prev=1E7;
				for(trials=0;trials<simits;trials++){ // Few trials -> roll the dice on starting orientation
					gsl_vector_set(A.guess,0,2*pi*gsl_rng_uniform_pos(r));gsl_vector_set(A.guess,1,pi*gsl_rng_uniform_pos(r));
					simplex(&sys,&A,2,&plane_orientation_find);
					if(A.err<z_prev){
						z_prev = A.err;
						sys.PHI0 = A.ans[0];
						sys.THETA0 = A.ans[1];}
				}
				printf("Trying to initiate propagation...\n");
				printf("z_prev: %.20lf\n",z_prev);
				theta00 = sys.THETA0;
				phi00 = sys.PHI0;
				
				YyY(theta00,phi00,Y00); // Create Y vector from theta and phi
				mat_mult(eang1,e_Ix,holdcm);
				gsl_blas_zgemv(CblasTrans,GSL_COMPLEX_ONE,holdcm,Y00,GSL_COMPLEX_ZERO,Q00); // Turn Y into Q
				for(o=0;o<3;o++){
					gsl_vector_complex_set(Q0,o,gsl_vector_complex_get(Q00,o)); // Save Q00 vector for when stepbacks go to <=0
					gsl_matrix_complex_set(Q_all,0,o,gsl_vector_complex_get(Q00,o));} // Q_all important for stepbacks
				sys.stride = stride; // Set the stride variable (can change depending on where in sequence you are)
				nvar = 2*sys.stride; // 2 angles for each stride
				gsl_vector_free(A.guess); // Free memory
				A.guess = gsl_vector_alloc(nvar); // Reallocate it
			}
			
			z0 = 1E10;
			tolerance = tolerance_all[n];
			
			// Torsion search
			int in;
            sys.Q0 = Q0; // So I don't have to pass Q0 as an argument to simplex subroutine
            for(o=0;o<simits;o++){
                for(p=0;p<nvar;p++) // Generate starting guesses for nvar variables
                    if(p%2==0) // If phi
                        gsl_vector_set(A.guess,p,gsl_ran_flat(r, PHI_min, PHI_max));
                    else // If psi
                        gsl_vector_set(A.guess,p,gsl_ran_flat(r, PSI_min, PSI_max));
                simplex(&sys,&A,nvar,&sigmaQ2); // Run simplex
                if(A.err < z0){ // Save low score
                    z0 = A.err;
                    in = o;
                    for(p=0;p<nvar;p++)
                        A.low[p] = A.ans[p];
                    if(z0<tol2) // Early break to speedup algorithm
                        break;}}
			
			if(z0<tolerance){
				
				// This line returns a tolerance for indices in which it my have been erroneously removed
				//if(z0<tol)
				//	for(o=0;o<sys.stride;o++)
				//		tolerance_all[n+o] = tol;
				
				rmsdf[n][S] = z0;
				printf("z0    : %lf (%d)",z0,in);
				
				// Don't let angles go outside of -180 < x < 180
				// TODO: Put at end of program
				for(o=0;o<nvar;o++){
					if(A.low[o] > pi)
						A.low[o] = -pi + fmod(A.low[o],pi);
					else if(A.low[o] < -pi)
						A.low[o] = pi + fmod(A.low[o],pi);}
				
				for(o=0;o<3;o++)
					gsl_vector_complex_set(holdcv, o, gsl_vector_complex_get(Q0, o));
				for(o=0;o<sys.stride;o++){
					sys.params[S][n+o][0] = A.low[2*o]; // phis
					sys.params[S][n+o][1] = A.low[(2*o)+1]; // psis
					measurables(&sys, holdcv, sys.params[S][n+o][0], sys.params[S][n+o][1], holdcv2, dummy);
					for(p=0;p<3;p++){
						gsl_matrix_complex_set(Q_all, n+o+1, p, gsl_vector_complex_get(holdcv2,p));
						gsl_vector_complex_set(holdcv, p, gsl_vector_complex_get(holdcv2, p));}}
				for(o=0;o<3;o++)
					gsl_vector_complex_set(Q0, o, gsl_vector_complex_get(holdcv, o));
                
				n+=sys.stride;}
			else{
				// Residue n failed
				stepbkV[0][n]++;
				if(stepbkV[0][n]==3){ // Residue n has failed <#> times
					stepbkV[1][n]++; // Residue n has arrested the search
					if( (stepbkV[1][n] == 1) ){
						tolerance_all[n] = 1E9; // Remove tolerance for residue n
						tolerance_all[n+1] = 1E9;}
					for(o=0;o<sz-1;o++){stepbkV[0][o]=0;} // Reset tally for next search
					n=0;
					printf("!!!RESTART!!!");}
				else{ // Stepback
					n -= stepbk;
					//tolerance_all[n] += 20;
					printf("STEPBACK");}
				if(n>0){ // Set Q to stepback residue
					for(o=0;o<3;o++)
						gsl_vector_complex_set(Q0,o,gsl_matrix_complex_get(Q_all,n,o));
					sys.stride = stride;
					nvar = 2*sys.stride;
					gsl_vector_free(A.guess);
					A.guess = gsl_vector_alloc(nvar);}
				else{ // Set Q to 0
					n = 0;
					cnt = 0;
					for(o=0;o<3;o++)
						gsl_vector_complex_set(Q0,o,gsl_vector_complex_get(Q00,o));}}	
		}
		
		// Back-calculation of spectrum
		// Recover Y0
        sys.THETA0 = theta00;
		sys.PHI0 = phi00;
		YyY(sys.THETA0,sys.PHI0,Y0);
		// Calculate spectrum using Y-basis
		csdc(&sys,Y0,sys.params[S],freq_calc);
		// Calculate spectral RMS values
        double sumtot=0, sumcs=0, sumdc=0, sumch=0, c1=0,c2=0,c3=0, sumrmsdf = 0;
		for(m=0;m<sz;m++){
			sumrmsdf += rmsdf[m][S];
			c1 = pow(sys.target[m][0]-freq_calc[m][0],2.0);
			c2 = pow(sys.target[m][1]-freq_calc[m][1],2.0);
			if(CH==1)
				c3 = pow(sys.target[m][2]-freq_calc[m][2],2.0);
            sumcs+=c1;
            sumdc+=c2;
            sumch+=c3;
            sumtot+=sqrt(c1+c2+c3);
		}
        RMS[S][0] = sqrt(sumcs/(double)sz);
		RMS[S][1] = sqrt(sumdc/(double)sz);
		RMS[S][2] = sqrt(sumch/(double)sz);
        RMS[S][3] = sqrt(sumtot/(double)sz);
        
		// Append to the end of output.csv, freqs_RMS.csv
		if(file==1){
			fp = fopen("output/output.csv","a");
			fprintf(fp,"%.10lf,%.10lf\n",sys.THETA0,sys.PHI0);
			fprintf(fp,"%.10lf",sys.params[S][0][0]);
            for(m=1;m<sz-1;m++)
				fprintf(fp,",%.10lf",sys.params[S][m][0]);
			fprintf(fp,"\n");
            fprintf(fp,"%.10lf",sys.params[S][0][1]);
			for(m=1;m<sz-1;m++)
				fprintf(fp,",%.10lf",sys.params[S][m][1]);
			fprintf(fp,"\n");
			fclose(fp);
			fp2 = fopen("output/freqs_RMS.csv","a");
            fprintf(fp2,"%lf,%lf,%lf,%lf\n",RMS[S][0],RMS[S][1],RMS[S][2],RMS[S][3]);
			fclose(fp2);
		}
		
		printf("\n\none more structure done...\nkb: %d\n",S);		
        printf("Freq RMS for structure %d\n",S);
        printf("%lf\n",RMS[S][3]);
		
	}
	end = clock();
	
	////////////////////////// FILE //////////////////////////////////
	if(file==1){
		// angles.csv
        double So=0;
		FILE * g = fopen("output/angles.csv","r");
        if(g){
            m = fscanf(g, "So,%lf",&So);
            fclose(g);
        }
        else
            printf("\nWarning: No output/angles.csv found. Must manually input correct order parameter (So) into file\n");
        g = fopen("output/angles.csv", "w");
        fprintf(g, "So,%.3lf\n",So);
		fprintf(g,"HNCa,%.2lf\nNCCa,%.2lf\nHNCo,%.2lf\na_nc,%.2lf\ntetra,%.2lf\ntetid,%.2lf\n",HNCa/d2r,NCCa/d2r,HNCo/d2r,a_nc/d2r,tetra/d2r,tetid/d2r);
		fprintf(g,"omN,%.2lf,%.2lf,%.2lf\n",gamma1/d2r,90.0,gamma2/d2r);
		fprintf(g,"omH,%.2lf,%.2lf,%.2lf\n",-90.0,-90.0,90.0);
		fprintf(g,"omHN,%.2lf,%.2lf,%.2lf\n",0.0,90.0,0.0);
		fprintf(g,"chi0,%lf\n",chi0);
		fprintf(g,"chi1,%lf",chi1);
		fclose(g);
		/*// Qs.csv
		if(kb_max==1){
			g = fopen("output/Qs.csv","w");
			for(m=0;m<sz;m++)
				fprintf(g,"%lf,%lf,%lf,%lf,%lf,%lf\n",Q_all->data[6*m],Q_all->data[(6*m)+1],Q_all->data[(6*m)+2],Q_all->data[(6*m)+3],Q_all->data[(6*m)+4],Q_all->data[(6*m)+5]);
			fclose(g);
		}*/
		/*// freqs_targ.csv
		g = fopen("output/freqs_targ.csv","w");
		for(m=0;m<sz;m++)
			fprintf(g,"%lf,%lf,%lf,%lf\n",sys.target[m][0],sys.target[m][1],sys.target[m][2],sys.target[m][3]);
		fclose(g);*/
		// sigout.csv
		g = fopen("output/sigout.csv", "w");
		for(m=0;m<sz;m++)
			fprintf(g, "%.2lf,%.2lf,%.2lf\n", sys.sigma[m][0], sys.sigma[m][1], sys.sigma[m][2]);
		fclose(g);
	}
    
    // Print runtime
	elapsed = ((double)end - (double)begin)/CLOCKS_PER_SEC;
	printf("\n\nElapsed time is %.3lf seconds\n",elapsed);
}
