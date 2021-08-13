/* The ziggurat method for RNOR and REXP
For details of the method, see Marsaglia and Tsang, "The ziggurat method
for generating random variables", Journ. Statistical Software.
*/
#ifndef RND_H
#define RND_H
#define MXW_SAMPLES 1000
#define MXW_FLUX_SAMPLES 1000
#define RNOR_MAX 127
#define REXP_MAX 255


#include <cmath>
#include <ctime>
#include "mymath.cpp"
#include <iostream>
using namespace std;

class t_random
{
    private:
	uint32 jz,jsr;
	int32 hz;
	uint32 iz, kn[128], ke[256];
	float wn[128],fn[128], we[256],fe[256];
	double mxw_data[MXW_SAMPLES];
	double mxw_flux_data[MXW_FLUX_SAMPLES];
	static constexpr double mxw_factor = 5.0;
        inline uint32 znew(){ return z=36969*(z&65535)+(z>>16); };
        inline uint32 wnew(){ return w=18000*(w&65535)+(w>>16); };

    public:
        inline uint32 shr3()
        {
            return jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5), jz+jsr;
        };
        inline uint32 iuni(){return shr3();};
        inline uint32 cong(){ return jcong=69069*jcong+1234567; };
        inline uint32 mwc(){ return (znew()<<16)+wnew();};
        inline uint32 kiss(){ return ((mwc()^cong())+shr3());};

	inline float uni(){return .5 + (signed) iuni() * 2.3283064365386963e-10; };
	inline float uni_small(){
                float tmp = .5 + (signed) iuni() * .2328306e-9;
                return tmp > 1e-5 ? tmp : uni() * 1e-5;};
	uint32 z, w, jcong;
	

        inline float rnor()
        {
            return hz=iuni(), iz=hz&127, ((unsigned)abs(hz)<kn[iz]) ? hz*wn[iz] : nfix();
        };
        //inline float rnor(){return normal2(); };
        inline float rexp()
        {
            return jz=iuni(), iz=jz&255, ( jz <ke[iz])? jz*we[iz] : efix();
        };
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

	void deflect(double angle, double &x, double &y, double &z)
	{
            double len = tan(angle/2.0);

            //generate random vector
            double x1, y1, z1;
	    double tmp = (1-2*uni());
	    x1 = len*tmp;
	    tmp = sqrt(1-sqr(tmp));

	    double sp,cp;
	    sincos(2*M_PI*uni(),&sp,&cp);

	    y1 = len*tmp*sp;
	    z1 = len*tmp*cp;

            //generate random rotation axis by vector multiplication
            double tx, ty, tz;
            tx = y*z1 - z*y1;
            ty = z*x1 - x*z1;
            tz = x*y1 - y*x1;

            tmp = len/norm(tx, ty, tz);
            tx *= tmp;
            ty *= tmp;
            tz *= tmp;

            // rotation
            //use Boris' algorithm (Birdsall & Langdon pp. 62) for arbitrary rotation
            double xprime = x - y*tz + z*ty;
            double yprime = y - z*tx + x*tz;
            double zprime = z - x*ty + y*tx;

            tmp = 2.0/(1+len*len);
            x1 = tx*tmp;
            y1 = ty*tmp;
            z1 = tz*tmp;

            x += - yprime*z1 + zprime*y1;
            y += - zprime*x1 + xprime*z1;
            z += - xprime*y1 + yprime*x1;

	};

        double radius() { return sqrt(uni()); };


	/*--------This procedure sets the seed and creates the tables------*/
        t_random()
        {
            uint32 jsrseed = time(NULL);
            cout << "seed = "<< jsrseed <<endl;
            initialize_seed(jsrseed);
            initialize_tables();
        }
        t_random(uint32 jsrseed)
        {
            cout << "seed = "<< jsrseed <<endl;
            initialize_seed(jsrseed);
            initialize_tables();
        }
        t_random(t_random & rnd)
        {
            initialize_seed(rnd.iuni());
            copy_tables(rnd);
        }

    private:
        void initialize_seed(uint32 jsrseed)
        {
            jsr=123456789;
            jsr^=(uint32)jsrseed;
            /* Set up tables for other generators */
            z=shr3();
            w=shr3();
            jcong=shr3(); 
        }
        void initialize_tables()
        {
            mymath_test();
	    const double m1 = 2147483648.0, m2 = 4294967296.;
	    double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
	    double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
	    int i;

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

	    double norm,v;
	    /* Set up tables for Maxwell1() */
	    norm=0;
	    for(i=1;i<MXW_SAMPLES;i++)
	    {
		v=(i-.5)*mxw_factor/MXW_SAMPLES;
		norm += mxw_data[i] = 4.0/sqrt(M_PI)*v*v*exp(-v*v);
	    }
	    mxw_data[0]=0;
	    for(i=1;i<MXW_SAMPLES;i++)
		mxw_data[i]=mxw_data[i]/norm + mxw_data[i-1];
	    mxw_data[MXW_SAMPLES-1]=1;

	    /* Set up tables for maxwell_flux() */
	    norm=0;
	    for(i=1;i<MXW_FLUX_SAMPLES;i++)
	    {
		v=(i-.5)*1000.0/MXW_FLUX_SAMPLES*2.0;
		norm += mxw_flux_data[i] = v*exp(-2.5e-5*v*v);
	    }
	    mxw_flux_data[0]=0;
	    for(i=1;i<MXW_FLUX_SAMPLES;i++)
		mxw_flux_data[i]=mxw_flux_data[i]/norm + mxw_flux_data[i-1];
	    mxw_flux_data[MXW_FLUX_SAMPLES-1]=1;
	}

        void copy_tables(const t_random & rnd)
        {
            int i;
            for(i=0; i<=RNOR_MAX; i++)
            {
                kn[i] = rnd.kn[i];
                wn[i] = rnd.wn[i];
                fn[i] = rnd.fn[i];
            }
            for(i=0; i<=REXP_MAX; i++)
            {
                ke[i] = rnd.ke[i];
                we[i] = rnd.we[i];
                fe[i] = rnd.fe[i];
            }
            for(i=0; i<MXW_SAMPLES; i++)
                mxw_data[i] = rnd.mxw_data[i];
            for(i=0; i<MXW_FLUX_SAMPLES; i++)
                mxw_flux_data[i] = rnd.mxw_flux_data[i];
        }

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

/*
//
// ***** alternative FIB random number generator *****
//
#define FIB   ((b=a+b),(a=b-a))
#define LFIB4 (c++,t[c]=t[c]+t[UC(c+58)]+t[UC(c+119)]+t[UC(c+178)])
#define SWB   (c++,bro=(x<y),t[c]=(x=t[UC(c+34)])-(y=t[UC(c+19)]+bro))
#define UC    (unsigned char)  //a cast operation
	typedef unsigned long UL;

	//  Global static variables: 
	 static UL z=362436069, w=521288629, jcong=380116160;
	 static UL a=224466889, b=7584631, t[256];
	// Use random seeds to reset z,w,jsr,jcong,a,b, and the table t[256]

	 static UL x=0,y=0,bro; static unsigned char c=0;

	// Example procedure to set the table, using KISS: 
	 void settable(UL i1,UL i2,UL i3,UL i4,UL i5, UL i6)
	 { int i; z=i1;w=i2,jsr=i3; jcong=i4; a=i5; b=i6;
	 for(i=0;i<256;i=i+1)  t[i]=KISS;
	 }
	 // initialization for FIB generators
        uint32 a, b, t[256];
	    a=shr3();
	    b=shr3();
	    for(int i=0;i<256;i=i+1)  t[i]=kiss();
*/

#endif
