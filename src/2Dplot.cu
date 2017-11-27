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
#include <GLUT/glut.h>
#include <thrust/complex.h>

#define EPS 0.000001  // 停止判定
#define MAXIT 30      // 最大反復回数
#define ZMAX 4.0      // 初期値の最大絶対値
#define ZOOM 200      // 拡大率
//#define RMAX 1000     // 複素平面の分割数（ iMacでは 2000 とする）
#define RMAX 500     // 複素平面の分割数（ iMacでは 2000 とする）

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
template<typename T> void Horner( std::vector< thrust::complex<T> > cf, thrust::complex<T> z,
					thrust::complex<T> &vf, thrust::complex<T> &df )
{
	vf = Cef[0];
	df = thrust::complex<T> (0.0,0.0);
	thrust::complex<T> tmp;

    for(auto itr = Cef.begin()+1; itr < Cef.end(); ++itr)
    {
    	tmp = vf;
    	vf = vf*z + *itr;
    	df = df*z + tmp;
    }
}

// Nourein subfunction
template <typename T> thrust::complex<T> vc( const int K, thrust::complex<T> z )
{
	thrust::complex<T> tmp = thrust::complex<T> (0.0,0.0);;

	for (auto itr = Zrs.begin(); itr < Zrs.end(); ++itr )
	{
		thrust::complex<T> vf, df;
		Horner( Cef, *itr, vf, df );

		// tmp *= (1/f'(z_i) (-1 / (z_i -z)^{K+1})
		tmp += ( (T)(1.0) / df )*( (T)(-1.0) / pow( (*itr - z), (T)(K+1) ));
	}
	return tmp;
}

template <typename T> thrust::complex<T> Nourein( const int p, thrust::complex<T> z, int &count, T &er )
{
	assert(p>=2);

	thrust::complex<T> vf, df;
	Horner( Cef, z, vf, df );
	count = 0;

	while ((count < MAXIT) && (abs(vf) > EPS))
	{
		z += vc(p-2,z) / vc(p-1,z);
		Horner( Cef, z, vf, df );
		count++;
	}
	er = abs(vf);

	return z;
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

template <typename T> void DrawApollonius( const int i, const int j, const T alp )
{
	const int pts = 180;    // 円周上の点数

	thrust::complex<T> center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
	double radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
	double tic = (double)(2.0*M_PI / pts);

	// Apollonius円の描画
	glColor3d(1.0,1.0,1.0);   // 白の円を描画
	glLineWidth(1.0);         // 線の太さ（ディフォルトは1.0）

	glBegin(GL_LINE_LOOP);
	for (int i=1; i<pts; i++)
	{
		glVertex2d( center.real() + radius*cos( tic*i )  , center.imag() + radius*sin( tic*i ) );
	}
	glEnd();
	glFlush();
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
	glClearColor(0.0, 0.0, 0.0, 0.0); // 塗りつぶしの色を指定
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	//////////////////////////////////////////////
	// 点の描画
	glBegin(GL_POINTS);

	double x = (double)(-ZMAX);
	for (int i=0; i<RMAX; i++)
	{
		double y = (double)(-ZMAX);
		for (int j=0; j<RMAX; j++)
		{
			int count;
			double er;
			thrust::complex<double> z0 = thrust::complex<double>( x, y );
			thrust::complex<double> z = Nourein(P,z0,count,er);

			int grad = 16;  // 明るさの階調
			double bright;
			if (count > grad)
				bright = 0.0;
			else
			{
				// 反復回数1回が最も明るく（bright=1）となるように修正（count-1）
				bright = double(grad - (count-1)) / double(grad);
			}

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
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	// 零点の描画
	for (auto itr = Zrs.begin(); itr < Zrs.end(); ++itr)
	{
		glColor3d(1.0,1.0,1.0);   // 白の点を描画
		glPointSize(8.0);      // 点の大きさ（ディフォルトは1.0)
		glBegin(GL_POINTS);
		glVertex2d( (*itr).real(), (*itr).imag() );
		glEnd();
	}
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	// Apolloniusの描画
	SetGamma( Gam );       // Γ
	GetAlpha( Gam, Alp );  // α
//	for (int i=0; i<Zrs.size(); i++)
	int i=0;
	{
		for (int j=0; j<Zrs.size(); j++)
		{
			if (i!=j)
				DrawApollonius(i,j,Alp[i]);
		}
	}
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
	glutInitWindowSize(1100,1100);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
