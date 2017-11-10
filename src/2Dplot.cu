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
#include <numeric>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <vector>

#include <thrust/complex.h>

#define EPS 0.0000001  // 停止判定
#define MAXIT 40    // 最大反復回数
#define ZMAX 1.5     // 初期値の最大絶対値
#define ZOOM 500     // 拡大率
#define RMAX 4000    // 複素平面の分割数
#define OFFS 0.3     // 明度のオフセット

#if defined (__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

template<typename T> thrust::complex<T> Update( thrust::complex<T> z )
{
	thrust::complex<T> vf = z*z*z -1.0;
	thrust::complex<T> df = 3.0*z*z;

	return z - vf / df;
}

template<typename T> thrust::complex<T> Newton( thrust::complex<T> z, int &count )
{
	double diff = (double)(MAXIT);

	count = 0;
	while ((count < MAXIT) && (diff > EPS))
	{
		thrust::complex<T> d = Update(z);
		diff = abs( d - z );
		count++;
		z = d;
//		std::cout << z << std::endl;
	}

	return z;
}

template<typename T> int FixPoint( thrust::complex<T> z )
{
	const int nfp = 3;  // 不動点の数
	thrust::complex<T> *fps = new thrust::complex<T> [nfp];

	fps[0] = thrust::complex<T> ( 1.0, 0.0 );
	fps[1] = thrust::complex<T> ( -0.5, 0.866025 );
	fps[2] = thrust::complex<T> (  0.5, 0.866025 );

	int col = 0;
	double min = (double)(MAXIT);

	for (int i=0; i<nfp; i++)
	{
		if (abs(z - fps[i]) < min)
		{
			min = abs(z - fps[0]);
			col = i;
		}
	}
	delete[] fps;

	return col;
}

void display(void)
{
	glClearColor(0.0, 0.0, 0.0, 1.0); // 塗りつぶしの色を指定（黒）
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	///////////////////////////////////
	// 座標軸の描画
	glBegin(GL_LINES);
	glColor3d(1.0, 0.0, 0.0);         // 赤 (1,0,0) で描画
	glVertex2d(-ZMAX,  0.0);
	glVertex2d( ZMAX,  0.0);

	glVertex2d( 0.0, -ZMAX);
	glVertex2d( 0.0,  ZMAX);

	glEnd();
	glFlush();                        // OpenGL命令のフラッシュ
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	// 点の描画
	glBegin(GL_POINTS);

	double x = (double)(-ZMAX);

//	#pragma omp parallel for
	for (int i=0; i<RMAX; i++)
	{
		int count;
		double y = (double)(-ZMAX);
		for (int j=0; j<RMAX; j++)
		{
			thrust::complex<double> z0 = thrust::complex<double>( x, y );
			thrust::complex<double> z = Newton(z0,count);
			double brit = (double)(1.0/MAXIT)*(MAXIT - count);
			//double brit = (1.0/MAXIT)*count + OFFS;

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
			default:
				glColor3d(0.0,0.0,0.0);
				break;
			}

			glVertex2d( z0.real(), z0.imag() );  // 点の描画
			if (count >= MAXIT-10)
			{
				std::cout << "z0 = " << z0 << ", count = " << count;
				std::cout << ", z = " << z << ", bright = " << brit;
				std::cout << ", color = " << FixPoint(z) << std::endl;
			}
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
	glutInitDisplayMode(GLUT_RGBA); // RGBモードに設定
	glutInitWindowSize(1000,1000);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
