/* The ziggurat method for RNOR and REXP
Combine the code below with the main program in which you want
normal or exponential variates.   Then use of RNOR in any expression
will provide a standard normal variate with mean zero, variance 1,
while use of REXP in any expression will provide an exponential variate
with density exp(-x),x>0.
Before using RNOR or REXP in your main, insert a command such as
zigset(86947731);
with your own choice of seed value>0, rather than 86947731.
(If you do not invoke zigset(...) you will get all zeros for RNOR and REXP.)
For details of the method, see Marsaglia and Tsang, "The ziggurat method
for generating random variables", Journ. Statistical Software.
*/
#ifndef RND_H
#define RND_H
#define MXW_SAMPLES 1000
#define MXW_FLUX_SAMPLES 1000


#include <math.h>
#include <ctime>
#include "mymath.cpp"
#include <cstdlib>

class t_random
{
    private:
	uint32 jz,jsr;
	//static unsigned long jz,jsr=123456789;
	double sqr(double x){return x*x;};

    public:
	uint32 shr3(){return jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5), jz+jsr;};
#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3 * .2328306e-9)
	inline uint32 iuni(){return shr3();};
	//inline unsigned int iuni(){return (unsigned int)rand();};

    private:
	int32 hz;
	//static long hz;
	uint32 iz, kn[128], ke[256];
	float wn[128],fn[128], we[256],fe[256];
	double mxw_data[MXW_SAMPLES];
	double mxw_flux_data[MXW_FLUX_SAMPLES];
	static const double mxw_factor = 5.0;

    public:
#define znew   (z=36969*(z&65535)+(z>>16))
#define wnew   (w=18000*(w&65535)+(w>>16))
//#define MWC    ((znew<<16)+wnew )
#define CONG  (jcong=69069*jcong+1234567)
	inline uint32 mwc(){ return (znew<<16)+wnew;};
	inline uint32 kiss(){ return ((mwc()^CONG)+SHR3);};
	inline float uni(){return .5 + (signed) iuni() * .2328306e-9; };
/************/
/*
//#define SHR3  (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define FIB   ((b=a+b),(a=b-a))
#define LFIB4 (c++,t[c]=t[c]+t[UC(c+58)]+t[UC(c+119)]+t[UC(c+178)])
#define SWB   (c++,bro=(x<y),t[c]=(x=t[UC(c+34)])-(y=t[UC(c+19)]+bro))
#define UNI   (KISS*2.328306e-10)
#define VNI   ((long) kiss)*4.656613e-10
#define UC    (unsigned char)  //a cast operation
	typedef unsigned long UL;
*/
	/*  Global static variables: */
//	 static UL z=362436069, w=521288629, jcong=380116160;
//	 static UL a=224466889, b=7584631, t[256];
	/* Use random seeds to reset z,w,jsr,jcong,a,b, and the table t[256]*/

//	 static UL x=0,y=0,bro; static unsigned char c=0;

	/* Example procedure to set the table, using KISS: */
//	 void settable(UL i1,UL i2,UL i3,UL i4,UL i5, UL i6)
//	 { int i; z=i1;w=i2,jsr=i3; jcong=i4; a=i5; b=i6;
///	 for(i=0;i<256;i=i+1)  t[i]=KISS;
//	 }
/********************/

	uint32 z, w, jcong, a, b, t[256];
	

	inline float rnor(){return hz=iuni(), iz=hz&127, ((unsigned)abs(hz)<kn[iz]) ? hz*wn[iz] : nfix(); };
	//inline float rnor(){return normal2(); };
	inline float rexp(){return jz=iuni(), iz=jz&255, ( jz <ke[iz])? jz*we[iz] : efix(); };
	//inline float rexp(){ return -log(uni());}


	double normal1()	{
	    double res=0;
	    for(int i=0;i<12;i++) res+=uni();
	    return res-6;
	}
	double normal2(){return sqrt(-2.0*log(uni()))*cos(2.0*M_PI*uni());}


	double maxwell1(double vMax)
	{
	    int j1,j2,j3;
	    double gama;

	    gama=uni();
	    j1=0;
	    j2=MXW_SAMPLES;
	    while((j2-j1)>1)
	    {
		j3=(j1+j2)/2;
		if(gama < mxw_data[j3]) j2=j3;
		else j1=j3;
	    }
	    return (j2-uni())*vMax*mxw_factor/MXW_SAMPLES;
	}
	double maxwell_flux(double vMax)
	{
	    int j1,j2,j3;
	    double gama;

	    gama=uni();
	    j1=0;
	    j2=MXW_FLUX_SAMPLES;
	    while((j2-j1)>1)
	    {
		j3=(j1+j2)/2;
		if(gama < mxw_flux_data[j3]) j2=j3;
		else j1=j3;
	    }
	    return (j2-uni())*vMax/MXW_FLUX_SAMPLES*2.0;
	}

	double maxwell2(double vMax){ return vMax * sqrt( -log(uni()) - log(uni()) * sqr(cos(2*M_PI*uni())) ); }
	double maxwell3(double vMax){ return vMax * sqrt( sqr(rnor()) + sqr(rnor()) + sqr(rnor()) ) * M_SQRT1_2; }
	void rot(double len, double &x, double &y, double &z)
	{
	    //nahodna zmena smeru
	    double cs_theta = (1-2*uni());
	    x = len*cs_theta;
	    cs_theta = sqrt(1-sqr(cs_theta));

	    double sp,cp;
	    sincos(2*M_PI*uni(),&sp,&cp);

	    y = len*cs_theta*sp;
	    z = len*cs_theta*cp;
	};

	void rot(double &x, double &y, double &z)
	{
	    //nahodna zmena smeru
	    double len = sqrt(sqr(x)+sqr(y)+sqr(z));

	    double cs_theta = (1-2*uni());
	    x = len*cs_theta;
	    cs_theta = sqrt(1-sqr(cs_theta));

	    double sp,cp;
	    sincos(2*M_PI*uni(),&sp,&cp);

	    y = len*cs_theta*sp;
	    z = len*cs_theta*cp;
	}


	/*--------This procedure sets the seed and creates the tables------*/
	t_random(uint32 jsrseed = -1) {
	    if(jsrseed < 0) jsrseed = time(NULL);
	    jsr=123456789;
	    const double m1 = 2147483648.0, m2 = 4294967296.;
	    double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
	    double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
	    int i;
	    jsr^=(uint32)jsrseed;
	    srand(jsrseed);
	    /* Set up tables for RNOR */
	    q=vn/exp(-.5*dn*dn);
	    kn[0]=(uint32)((dn/q)*m1);	  kn[1]=0;
	    wn[0]=q/m1;		  wn[127]=dn/m1;
	    fn[0]=1.;		  fn[127]=exp(-.5*dn*dn);		
	    for(i=126;i>=1;i--)
	    {
		dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
		kn[i+1]=(uint32)((dn/tn)*m1);
		tn=dn;
		fn[i]=exp(-.5*dn*dn);
		wn[i]=dn/m1;
	    }
	    /* Set up tables for REXP */	  q = ve/exp(-de);
	    ke[0]=(uint32)((de/q) * m2);	  ke[1]=0;
	    we[0]=q/m2;		  we[255]=de/m2;
	    fe[0]=1.;		  fe[255]=exp(-de);		
	    for(i=254;i>=1;i--)
	    {
		de=-log(ve/de+exp(-de));
		ke[i+1]= (uint32)((de/te)*m2);
		te=de;
		fe[i]=exp(-de);
		we[i]=de/m2;
	    }
	    /* Set up tables for Maxwell1() */
	    double norm,v;
	    const double mxw_factor = 5.0;

	    norm=0;
	    for(i=1;i<MXW_SAMPLES;i++)
	    {
		v=(i-.5)*mxw_factor/MXW_SAMPLES;
		norm += mxw_data[i] = 4.0/sqrt(M_PI)*v*v*exp(-v*v);
	    }
	    mxw_data[0]=0;
	    //cout << norm*mxw_factor/MXW_SAMPLES << endl;
	    for(i=1;i<MXW_SAMPLES;i++)
		mxw_data[i]=mxw_data[i]/norm + mxw_data[i-1];
	    mxw_data[MXW_SAMPLES-1]=1;


	    /* Set up tables for maxwell_flux() */
	    norm=0;
	    for(i=1;i<MXW_FLUX_SAMPLES;i++)
	    {
		v=(i-.5)*1000.0/MXW_FLUX_SAMPLES*2.0;
		norm += mxw_flux_data[i] = v*exp(-2.5e-5*v*v);
		//      printf("%f\n",mxw_data[i]);
	    }
	    mxw_flux_data[0]=0;
	    for(i=1;i<MXW_FLUX_SAMPLES;i++)
		mxw_flux_data[i]=mxw_flux_data[i]/norm + mxw_flux_data[i-1];
	    mxw_flux_data[MXW_FLUX_SAMPLES-1]=1;


	    /* Set up tables for other generators */
	    z=shr3();
	    w=shr3();
	    jcong=shr3(); 
	    a=shr3();
	    b=shr3();
	    for(int i=0;i<256;i=i+1)  t[i]=kiss();



	}

    private:
	/* nfix() generates variates from the residue when rejection in RNOR occurs. */
	float nfix(void) 
	{
	    const float r = 3.442620f; 	/* The starting of the right tail */	
	    static float x, y;
	    for(;;)
	    {
		x=hz*wn[iz];
		if(iz==0)
		{	
		    /* iz==0, handle the base strip */
		    do{
			x=-log(uni())*0.2904764;   /* .2904764 is 1/r */
			y=-log(uni());			
		    } while(y+y<x*x);
		    return (hz>0)? r+x : -r-x;		
		}
		/* iz>0, handle the wedges of other strips */
		if( fn[iz]+uni()*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;
		/* start all over */
		hz=shr3();
		iz=hz&127;
		/*if(abs(hz)<kn[iz]) return (hz*wn[iz]);	*/
		if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
	    }
	}

	/* efix() generates variates from the residue when rejection in REXP occurs. */

	float efix(void)
	{	
	    float x;
	    for(;;)
	    {
		if(iz==0) return (7.69711-log(uni()));		/* iz==0 */
		x=jz*we[iz];
		if( fe[iz]+uni()*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);
		/* Start all over again */
		jz=shr3();
		iz=(jz&255);
		if(jz<ke[iz]) return (jz*we[iz]);
	    }
	}
};
#endif
