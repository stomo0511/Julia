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

#include <thrust/complex.h>

#define EPS 0.000001  // 停止判定
#define MAXIT 40      // 最大反復回数
#define ZMAX 4.0      // 初期値の最大絶対値
#define ZOOM 200      // 拡大率
#define RMAX 2000     // 複素平面の分割数

#if defined (__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

template<typename T> thrust::complex<T> vf( thrust::complex<T> z )
{
	thrust::complex<T> iu = thrust::complex<T> ( 0.0, 1.0 );
	return z*z*z*z*z + iu*z*z*z*z + 3.0*z*z*z + 41.0*iu*z*z + 132.0*z -90.0*iu;
}

template<typename T> thrust::complex<T> df( thrust::complex<T> z )
{
	thrust::complex<T> iu = thrust::complex<T> ( 0.0, 1.0 );
	return 5.0*z*z*z*z + 4.0*iu*z*z*z + 9.0*z*z + 82.0*iu*z + 132.0;
}

template<typename T> thrust::complex<T> Newton( thrust::complex<T> z, int &count )
{
	count = 0;

	while ((count < MAXIT) && (abs(vf(z)) > EPS))
	{
		z -= vf(z) / df(z);
		count++;
	}

	return z;
}

template<typename T> int FixPoint( thrust::complex<T> z )
{
	const int nfp = 5;  // 不動点の数
	thrust::complex<T> *fps = new thrust::complex<T> [nfp];

	fps[0] = thrust::complex<T> (  0.0,  1.0 );
	fps[1] = thrust::complex<T> (  1.0,  2.0 );
	fps[2] = thrust::complex<T> ( -1.0,  2.0 );
	fps[3] = thrust::complex<T> (  3.0, -3.0 );
	fps[4] = thrust::complex<T> ( -3.0, -3.0 );

	int col = 0;
	double min = 0.001;

	for (int i=0; i<nfp; i++)
	{
		if (abs(z - fps[i]) < min)
		{
			min = abs(z - fps[i]);
			col = i;
		}
	}
	delete[] fps;

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
		int count;
		double y = (double)(-ZMAX);
		for (int j=0; j<RMAX; j++)
		{
			thrust::complex<double> z0 = thrust::complex<double>( x, y );
			thrust::complex<double> z = Newton(z0,count);
//			double brit = (double)(1.0/MAXIT)*(MAXIT - count);

			double brit;
			if (count > 13)
				brit = 0.0;
			else
				brit = (13.0 - double(count)) / 13.0;
			// 明るさの補正
			brit += 0.1;

			switch( FixPoint(z) )  // 塗りつぶし色の設定
			{
			case 0:
				glColor3d(brit,0.0,0.0);
				break;
			case 1:
				glColor3d(0.0,brit,0.0);
				break;
			case 2:
				glColor3d(0.0,0.0,brit);
				break;
			case 3:
				glColor3d(brit,0.0,brit);
				break;
			case 4:
				glColor3d(0.0,brit,brit);
				break;
			default:
				glColor3d(0.0,0.0,0.0);
				break;
			}

			glVertex2d( z0.real(), z0.imag() );  // 点の描画
//			if (count >= MAXIT-1)
//			{
//				std::cout << "z0 = " << z0 << ", count = " << count;
//				std::cout << ", z = " << z << ", bright = " << brit;
//				std::cout << ", color = " << FixPoint(z) << std::endl;
//			}
			y += (double)(2*ZMAX / RMAX);
		}
		x += (double)(2*ZMAX / RMAX);
	}
	glEnd();
	glFlush();
	//////////////////////////////////////////////
}

void resize(int w, int h)
{
	// Window全体をView portにする
	glViewport(0,0,w,h);

	// 変換行列の初期化
	glLoadIdentity();

	// Screen上の表示領域をView portの大きさに比例させる
	glOrtho( -w/ZOOM, w/ZOOM, -h/ZOOM, h/ZOOM, -1.0, 1.0);
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);          // OpenGL初期化
	glutInitWindowSize(1100,1100);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
