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
#include <vector>
#include <thrust/complex.h>

#if defined (__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

#define EPS 0.000001  // 停止判定
#define MAXIT 30      // 最大反復回数
#define ZMAX 0.8      // 初期値の最大絶対値
#define ZOOM 600      // 拡大率
#define RMAX 1500     // 複素平面の分割数
#define XOFS 1.0      // x軸のオフセット

// Zeros
std::vector< thrust::complex<double> > Zrs {
	thrust::complex<double> (  1.0,  0.0 ),
	thrust::complex<double> (  0.0,  0.5 ),
	thrust::complex<double> (  0.0, -0.5 ),
	thrust::complex<double> ( -1.0,  0.0 )
};

// Coefficients
std::vector< thrust::complex<double> > Cef {
	thrust::complex<double> (      1.0, 0.0 ),  // Z^4
	thrust::complex<double> (      0.0, 0.0 ),  // Z^3
	thrust::complex<double> ( -3.0/4.0, 0.0 ),  // z^2
	thrust::complex<double> (      0.0, 0.0 ),  // z^1
	thrust::complex<double> ( -1.0/4.0, 0.0 )   // z^0
};

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

template<typename T> thrust::complex<T> Newton( thrust::complex<T> z, int &count, double &er )
{
	thrust::complex<T> vf, df;
	Horner( Cef, z, vf, df );

	count = 0;
	while ((count < MAXIT) && (abs(vf) > EPS))
	{
		z -= vf / df;
		Horner( Cef, z, vf, df );
		count++;
	}
	er = abs(vf);

	return z;
}

template<typename T> int FixPoint( thrust::complex<T> z )
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

template<typename T> void DrawApollonius( int i, int j, T alp )
{
	const int pts = 180;    // 円周上の点数

	thrust::complex<T> center = (Zrs[i] - alp*alp*Zrs[j]) / (1.0 - alp*alp);
	T radius = alp*abs(Zrs[i] - Zrs[j]) / (1.0 - alp*alp);
	T tic = (T)(2.0*M_PI / pts);

	// Z_i の描画
	glColor3d(1.0,1.0,1.0);   // 白の点を描画
	glPointSize(8.0);      // 点の大きさ（ディフォルトは1.0)
	glBegin(GL_POINTS);
	glVertex2d( Zrs[i].real() -XOFS, Zrs[i].imag() );
	glEnd();

	// Apollonius円の描画
	glColor3d(1.0,1.0,1.0);   // 白の円を描画
	glLineWidth(2.0);         // 線の太さ（ディフォルトは1.0）

	glBegin(GL_LINE_LOOP);
	for (int i=1; i<pts; i++)
	{
		glVertex2d( center.real() + radius*cos( tic*i ) -XOFS  , center.imag() + radius*sin( tic*i ) );
	}
	glEnd();
	glFlush();
}

void display(void)
{
	//////////////////////////////////////////////
	// 背景を白に
	glClearColor(1.0, 1.0, 1.0, 1.0); // 塗りつぶしの色を指定（黒）
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	//////////////////////////////////////////////
	// 点の描画
	glBegin(GL_POINTS);

	double x = (double)(-ZMAX + XOFS);
	for (int i=0; i<RMAX; i++)
	{
		double y = (double)(-ZMAX);
		for (int j=0; j<RMAX; j++)
		{
			int count;
			double er;
			//double p;
			thrust::complex<double> z0 = thrust::complex<double>( x, y );
			thrust::complex<double> z = Newton(z0,count,er);

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
			default:
				glColor3d(0.0,0.0,0.0);
				break;
			}

			glVertex2d( z0.real()-XOFS, z0.imag() );  // 点の描画（原点補正あり）

			y += (double)(2.0*ZMAX / RMAX);
		}
		x += (double)(2.0*ZMAX / RMAX);
	}

	glEnd();
	glFlush();
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	// Apollonius円の描画
	double alp = (double)1.0 / (2 * Zrs.size() - 3.0);

	int j=0;
	for (auto itr = Zrs.begin(); itr < Zrs.end(); ++itr )
	{
		DrawApollonius( 0, j, alp );
		j++;
	}
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
	glutInit(&argc, argv);          // OpenGL初期化
	glutInitWindowSize(500,500);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
