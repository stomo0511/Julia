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
#include <GLUT/glut.h>
#include <thrust/complex.h>
#include <thrust/host_vector.h>

#define EPS 0.000001  // 停止判定
#define MAXIT 30      // 最大反復回数
#define ZMAX 4.0      // 初期値の最大絶対値
#define ZOOM 200      // 拡大率
#define RMAX 2000     // 複素平面の分割数
#define ORD  32        // Nourein法の次数

#define NFP 5 // 零点の数
thrust::complex<double> fps[NFP];

void setZero( thrust::complex<double> *fps )
{
	fps[0] = thrust::complex<double> (  0.0,  1.0 );
	fps[1] = thrust::complex<double> (  1.0,  2.0 );
	fps[2] = thrust::complex<double> ( -1.0,  2.0 );
	fps[3] = thrust::complex<double> (  3.0, -3.0 );
	fps[4] = thrust::complex<double> ( -3.0, -3.0 );
}

// Polynomial
thrust::complex<double> vf( thrust::complex<double> z )
{
	thrust::complex<double> iu = thrust::complex<double> ( 0.0, 1.0 );
	return z*z*z*z*z + iu*z*z*z*z + + 3.0*z*z*z + 41.0*iu*z*z + 132.0*z -90.0*iu;
}

// derived function of the polynomial
thrust::complex<double> df( thrust::complex<double> z )
{
	thrust::complex<double> iu = thrust::complex<double> ( 0.0, 1.0 );
	return 5.0*z*z*z*z + 4.0*iu*z*z*z + 9.0*z*z + 82.0*iu*z + 132.0;
}

// Nourein subfunction
thrust::complex<double> vc( const int K, thrust::complex<double> z )
{
	thrust::complex<double> f = thrust::complex<double> (0.0,0.0);;

	for (int i=0; i<NFP; i++)
	{
		thrust::complex<double> tmp = thrust::complex<double> (1.0,0.0);

		// tmp = (z_i -z)^{k+1}
		for (int k=0; k<=K; k++)
		{
			tmp = tmp * (fps[i] - z);
		}
		// tmp = -1.0 /  (z_i -z)^{k+1}
		tmp = -1.0 / tmp;

		f += ( 1.0 / df(fps[i]) )*tmp;
	}
	return f;
}

thrust::complex<double> Nourein( const int p, thrust::complex<double> z, int &count, double &er )
{
	assert(p>=2);

	count = 0;

	while ((count < MAXIT) && (abs(vf(z)) > EPS))
	{
		z += vc(p-2,z) / vc(p-1,z);
		count++;
	}
	er = abs(vf(z));

	return z;
}

void DrawApollonius( int i, int j, double alp )
{
	const int pts = 180;    // 円周上の点数

	thrust::complex<double> center = (fps[i] - alp*alp*fps[j]) / (1.0 - alp*alp);
	double radius = alp*abs(fps[i] - fps[j]) / (1.0 - alp*alp);
	double tic = (double)(2.0*M_PI / pts);

	// Z_i の描画
	glColor3d(1.0,1.0,1.0);   // 白の点を描画
	glPointSize(8.0);      // 点の大きさ（ディフォルトは1.0)
	glBegin(GL_POINTS);
	glVertex2d( fps[i].real(), fps[i].imag() );
	glEnd();

	// Apollonius円の描画
	glColor3d(1.0,1.0,1.0);   // 白の円を描画
	glLineWidth(2.0);         // 線の太さ（ディフォルトは1.0）

	glBegin(GL_LINE_LOOP);
	for (int i=1; i<pts; i++)
	{
		glVertex2d( center.real() + radius*cos( tic*i )  , center.imag() + radius*sin( tic*i ) );
	}
	glEnd();
	glFlush();
}

int FixPoint( thrust::complex<double> z )
{
	int col = 0;
	double min = (double)MAXIT;

	for (int i=0; i<NFP; i++)
	{
		if (abs(z - fps[i]) < min)
		{
			min = abs(z - fps[i]);
			col = i;
		}
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

	double x = (double)(-ZMAX);
	for (int i=0; i<RMAX; i++)
	{
		double y = (double)(-ZMAX);
		for (int j=0; j<RMAX; j++)
		{
			int count;
			double er;
			thrust::complex<double> z0 = thrust::complex<double>( x, y );
			thrust::complex<double> z = Nourein(ORD,z0,count,er);

			int grad = 16;  // 明るさの階調
			double bright;
			if (count > grad)
				bright = 0.0;
			else
			{
				// 反復回数1回が最も明るく（bright=1）となるように修正（count-1）
				bright = double(grad - (count-1)) / double(grad);
			}
			// 明るさの補正
//			bright += 0.2;

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
	for (int i=0; i<NFP; i++)
	{
		glColor3d(1.0,1.0,1.0);   // 白の点を描画
		glPointSize(8.0);      // 点の大きさ（ディフォルトは1.0)
		glBegin(GL_POINTS);
		glVertex2d( fps[i].real(), fps[i].imag() );
		glEnd();
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
	setZero(fps);     // 零点のセット

	glutInit(&argc, argv);          // OpenGL初期化
	glutInitWindowSize(1100,1100);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
