/*
 ============================================================================
 Name        : 2Dplot.cu
 Author      : stomo
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>

#if defined (__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

#include <thrust/complex.h>

#define EPS 0.000001  // 停止判定
#define MAXIT 16      // 最大反復回数
#define ZMAX 4.0      // 初期値の最大絶対値
#define ZOOM 200      // 拡大率
#define RMAX 1000     // 複素平面の分割数

int P;  // Nourein法の次数

// Zeros
std::vector< thrust::complex<double> > Zrs {
	thrust::complex<double> (  0.0,  1.0 ),
	thrust::complex<double> (  1.0,  2.0 ),
	thrust::complex<double> ( -1.0,  2.0 ),
	thrust::complex<double> (  3.0, -3.0 ),
	thrust::complex<double> ( -3.0, -3.0 )
};

// Coefficients
std::vector< thrust::complex<double> > Cef {
	thrust::complex<double> (  1.0,   0.0 ),  // z^5
	thrust::complex<double> (  0.0,   1.0 ),  // Z^4
	thrust::complex<double> (  3.0,   0.0 ),  // Z^3
	thrust::complex<double> (  0.0,  41.0 ),  // z^2
	thrust::complex<double> (132.0,   0.0 ),  // z^1
	thrust::complex<double> (  0.0, -90.0 )   // z^0
};

std::vector<double> Gam( Zrs.size() );  // Γ
std::vector<double> Alp( Zrs.size() );  // α

// Hornet method for polynomial
template<typename T> void Horner(
		const std::vector< thrust::complex<T> > cf,
		const thrust::complex<T> z,
		thrust::complex<T> &vf, thrust::complex<T> &df )
{
	vf = cf[0];
	df = thrust::complex<T> (0.0,0.0);
	thrust::complex<T> tmp;

    for(auto itr = cf.begin()+1; itr < cf.end(); ++itr)
    {
    	tmp = vf;
    	vf = vf*z + *itr;
    	df = df*z + tmp;
    }
}

// Nourein subfunction
template <typename T> thrust::complex<T> Chi_m(
		const std::vector< thrust::complex<T> > zr,
		const std::vector< thrust::complex<T> > cf,
		const int m,
		const thrust::complex<T> z )
{
	thrust::complex<T> tmp = thrust::complex<T> (0.0,0.0);

	for (auto itr = zr.begin(); itr < zr.end(); ++itr )
	{
		thrust::complex<T> vf, df;
		Horner( cf, *itr, vf, df );

		// tmp *= (1/f'(z_i) (-1 / (z_i -z)^{m+1})
		tmp += ( (T)(1.0) / df )*( (T)(-1.0) / pow( (*itr - z), (T)(m+1) ));
	}
	return tmp;
}

template <typename T> thrust::complex<T> Nourein(
		const std::vector< thrust::complex<T> > zr,
		const std::vector< thrust::complex<T> > cf,
		const int p,
		thrust::complex<T> z,
		int &count, T &er )
{
	assert(p>=2);

	thrust::complex<T> vf, df;
	Horner( cf, z, vf, df );
	count = 0;

	while ((count < MAXIT) && (abs(vf) > EPS))
	{
		z += Chi_m(zr,cf,p-2,z) / Chi_m(zr,cf,p-1,z);
		Horner( cf, z, vf, df );
		count++;
	}
	er = abs(vf);

	return z;
}

// 円の描画
template <typename T> void Circle2D(T r,int x,int y)
{
	for (T th1 = 0.0;  th1 <= 360.0;  th1 = th1 + 1.0)
	{
		T th2 = th1 + 10.0;
		T th1_rad = th1 / 180.0 * M_PI;
		T th2_rad = th2 / 180.0 * M_PI;

		T x1 = r * cos(th1_rad);
		T y1 = r * sin(th1_rad);
		T x2 = r * cos(th2_rad);
		T y2 = r * sin(th2_rad);

		glBegin(GL_LINES);
		glVertex2f( x1+x, y1+y );
		glVertex2f( x2+x, y2+y );
		glEnd();
	}
}

// 円の描画（塗りつぶし）
template <typename T> void Circle2DFill(T r,int x,int y)
{
	for (T th1 = 0.0;  th1 <= 360.0;  th1 = th1 + 1.0)
	{
		T th2 = th1 + 10.0;
		T th1_rad = th1 / 180.0 * M_PI;
		T th2_rad = th2 / 180.0 * M_PI;

		T x1 = r * cos(th1_rad);
		T y1 = r * sin(th1_rad);
		T x2 = r * cos(th2_rad);
		T y2 = r * sin(th2_rad);

		glBegin(GL_TRIANGLES);
		glVertex2f( x, y );
		glVertex2f( x1+x, y1+y );
		glVertex2f( x2+x, y2+y );
		glEnd();
	}
}

template <typename T> void SetGamma( std::vector<T> &Gam )
{
	for (int i=0; i<Zrs.size(); i++)
	{
		T max = (T)(0.0);
		for (int j=0; j<Zrs.size(); j++)
		{
			if (i != j)
			{
				thrust::complex<T> vf, dfi, dfj;

				Horner( Cef, Zrs[i], vf, dfi );
				Horner( Cef, Zrs[j], vf, dfj );

				if (max < abs(dfi / dfj))
				{
					max = abs(dfi / dfj);
				}
			}
		}
		Gam[i] = max;
	}
}

template <typename T> T fAlp( const int p, const T Gamma, T alp )
{
	assert(p>=2);

	return ((T)(Zrs.size()) -1.0) * Gamma * pow(alp,(T)(p-1)) - (1.0-alp)/(1.0+3.0*alp);
}

template <typename T> T dAlp( const int p, const T Gamma, T alp )
{
	assert(p>=2);

	T tmp = ((T)(Zrs.size()) -1.0) * Gamma * ((T)(p) -1.0) * pow(alp,(T)(p-2));
	return  tmp + 4.0/(1.0+6.0*alp+9.0*alp*alp);
}

template <typename T> void GetAlpha( const std::vector<T> Gam, std::vector<T> &Alp )
{
	for (int i=0; i<Zrs.size(); i++)
	{
		T alp = (T)(1.0);
		int count = 0;

		while ((count < MAXIT) && (abs(fAlp(P,Gam[i],alp)) > EPS))
		{
			alp -= fAlp(P,Gam[i],alp) / dAlp(P,Gam[i],alp);
			count++;
		}
		if (count == MAXIT)
		{
			std::cerr << "No convergence in GetAlpha\n";
			std::exit (EXIT_FAILURE);
		}
		Alp[i] = alp;
	}
}

// Apollonius円の描画
template <typename T> void DrawApollonius( const int i, const int j, const T alp )
{
	assert(alp<1.0);

	thrust::complex<T> center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
	double radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);

	Circle2D( radius, (int)(center.real()), (int)(center.imag()) );
}

// Apollonius領域の描画
template <typename T> void DrawApRegion( const T alp )
{
	const int pts = 20;    // 円周上の点数
	thrust::complex<T> center;
	T radius;
	T start, end, tic;

	int i, j;

	// p=2
	{
		// z_0, z_1
		i = 0; j=1;
		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
		start = -1.41075;
		tic = -2.0*start / pts;

		glBegin(GL_LINE_STRIP);
		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
		for (int t=1; t<=pts; t++)
		{
			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
		}
		glEnd();

		// z_0, z_2
		i = 0; j=2;
		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
		start = 1.73084;
		tic = 2.0*(M_PI - start) / pts;

		glBegin(GL_LINE_STRIP);
		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
		for (int t=1; t<=pts; t++)
		{
			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
		}
		glEnd();
	}

	// p=4
//	{
//		// z_0, z_1
//		i = 0; j=1;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = -1.22945;
//		tic = -2.0*start / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_2
//		i = 0; j=2;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = 1.91214;
//		tic = 2.0*(M_PI - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//	}

	// p=8
//	{
//		// z_0, z_1
//		i = 0; j=1;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = -1.08417;
//		tic = -2.0*start / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_2
//		i = 0; j=2;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = 2.05743;
//		tic = 2.0*(M_PI - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//	}

	// p=16
//	{
//		// z_0, z_1
//		i = 0; j=1;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = -0.0309316;
//		end   = 0.976624;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_2
//		i = 0; j=2;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start =  2.16497;
//		end   =  2*M_PI - 3.11066;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_3
//		i = 0; j=3;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = -1.0758;
//		end   = -0.928056;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_4
//		i = 0; j=4;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = -2.21354;
//		end   = -2.06579;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//	}

	// p=32
//	{
//		// z_0, z_1
//		i = 0; j=1;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = 0.325413;
//		end   = 0.902587;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_2
//		i = 0; j=2;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start =  2.23901;
//		end   =  2.81618;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_3
//		i = 0; j=3;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = -0.903981;
//		end   = -1.01722;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//
//		// z_0, z_4
//		i = 0; j=4;
//		center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
//		radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
//		start = -2.12437;
//		end   = -2.23761;
//		tic = (end - start) / pts;
//
//		glBegin(GL_LINE_STRIP);
//		glVertex2d( center.real() + radius*cos(start), center.imag() + radius*sin(start) );
//		for (int t=1; t<=pts; t++)
//		{
//			glVertex2d( center.real() + radius*cos( start+tic*t )  , center.imag() + radius*sin( start+tic*t ) );
//		}
//		glEnd();
//	}
}

template <typename T> int FixPoint( thrust::complex<T> z )
{
	int i = 0;
	int col = 0;
	double min = (double)(MAXIT);

	for (auto itr = Zrs.begin(); itr < Zrs.end(); ++itr )
	{
		if (abs( z - *itr) < min)
		{
			min = abs( z - *itr);
			col = i;
		}
		i++;
	}
	return col;
}

void display(void)
{
	//////////////////////////////////////////////
	// 背景を白に
	glClearColor(1.0, 1.0, 1.0, 1.0); // 塗りつぶしの色を指定
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	//////////////////////////////////////////////
	// 点の描画
	glBegin(GL_POINTS);

	double min = (double)(MAXIT);
	double max = 0.0;

	double x = (double)(-ZMAX);
	for (int i=0; i<RMAX; i++)
	{
		double y = (double)(-ZMAX);
		for (int j=0; j<RMAX; j++)
		{
			int count;
			double er;
			thrust::complex<double> z0 = thrust::complex<double>( x, y );
			thrust::complex<double> z = Nourein(Zrs,Cef,P,z0,count,er);

			double bright;
			if (count > MAXIT)
				bright = 0.0;
			else
			{
				// 反復回数1回が最も明るく（bright=1）となるように修正（count-1）
				//bright = double(MAXIT - (count-1)) / double(MAXIT);
				bright = double(MAXIT - (count)) / double(MAXIT);
			}
//			std::cout << bright << std::endl;
			if (bright > max)
				max = bright;
			if (bright < min)
				min = bright;

			switch( FixPoint(z) )  // 塗りつぶし色の設定
			{
			case 0:
				glColor3d(bright,0.0,0.0);
				break;
			case 1:
				glColor3d(0.0,bright,0.0);
				break;
			case 2:
				glColor3d(0.0,0.0,bright);
				break;
			case 3:
				glColor3d(bright,0.0,bright);
				break;
			case 4:
				glColor3d(0.0,bright,bright);
				break;
			default:
				glColor3d(0.0,0.0,0.0);
				break;
			}

			glVertex2d( z0.real(), z0.imag() );  // 点の描画
			y += (double)(2*ZMAX / RMAX);
		}
		x += (double)(2*ZMAX / RMAX);
	}
	glEnd();
//	std::cout << "min = " << min << ", max = " << max << std::endl;
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	// 零点の描画
	glColor3d(1.0,1.0,1.0);   // 白の点を描画
	glLineWidth(1.0);
	for (auto itr = Zrs.begin(); itr < Zrs.end(); ++itr)
		Circle2DFill( (double)(0.05), (*itr).real(), (*itr).imag() );
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	// Apolloniusの描画
	SetGamma( Gam );       // Γ
	GetAlpha( Gam, Alp );  // α

	glColor3d(1.0,1.0,1.0);   // 白の円を描画
	glLineWidth(1.0);         // 線の太さ（ディフォルトは1.0）
	for (int i=0; i<Zrs.size(); i++)
	{
		for (int j=0; j<Zrs.size(); j++)
		{
			if (i!=j)
				DrawApollonius(i,j,Alp[i]);
		}
	}
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	// Apollonius領域の描画
	SetGamma( Gam );       // Γ
	GetAlpha( Gam, Alp );  // α

	glColor3d(1.0,1.0,1.0);   // 白の円を描画
	glLineWidth(2.0);         // 線の太さ（縮小時は4.0にする）
	DrawApRegion( Alp[0] );
	//////////////////////////////////////////////

	glFlush();
}

void resize(int w, int h)
{
	// Window全体をView portにする
	glViewport(0,0,w,h);

	// 変換行列の初期化
	glLoadIdentity();

	// Screen上の表示領域をView portの大きさに比例させる
	glOrtho( -(double)w/ZOOM, (double)w/ZOOM, -(double)h/ZOOM, (double)h/ZOOM, -1.0, 1.0);
}

int main(int argc, char *argv[])
{
	if (argc<2)
	{
		std::cerr << "Usage: a.out [Order]\n";
		exit (EXIT_FAILURE);
	}
	P = atoi(argv[1]);  // Nourein法の次数

	glutInit(&argc, argv);          // OpenGL初期化
	glutInitWindowSize(1000,1000);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
