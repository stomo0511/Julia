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

#define EPS 0.0000001  // 停止判定
#define MAXIT 40      // 最大反復回数
#define ZMAX 1.5      // 初期値の最大絶対値
#define ZOOM 500      // 拡大率
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
	return z*z*z -1.0;
}

template<typename T> thrust::complex<T> df( thrust::complex<T> z )
{
	return 3.0*z*z;
}

template<typename T> thrust::complex<T> Newton( thrust::complex<T> z, int &count, double &er )
{
	count = 0;

	while ((count < MAXIT) && (abs(vf(z)) > EPS))
	{
		z -= vf(z) / df(z);
		count++;
	}
	er = abs(vf(z));

	return z;
}

template<typename T> int FixPoint( thrust::complex<T> z )
{
	const int nfp = 3;  // 不動点の数
	thrust::complex<T> *fps = new thrust::complex<T> [nfp];

	fps[0] = thrust::complex<T> ( 1.0, 0.0 );
	fps[1] = thrust::complex<T> ( -0.5,  0.866025 );
	fps[2] = thrust::complex<T> ( -0.5, -0.866025 );

	int col = 0;
	double min = (double)(MAXIT);

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
		double y = (double)(-ZMAX);
		for (int j=0; j<RMAX; j++)
		{
			int count;
			double er, p;
			thrust::complex<double> z0 = thrust::complex<double>( x, y );
			thrust::complex<double> z = Newton(z0,count,er);

			//			double brit = (double)(1.0/MAXIT)*(MAXIT - count);

//			p = 0.0;
//			if (count < MAXIT)
//			{
//				p = log2( -12.0 / log10(er));
//			}
//			double brit = p;

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
			default:
				glColor3d(0.0,0.0,0.0);
				break;
			}

			glVertex2d( z0.real(), z0.imag() );  // 点の描画
			//if (count >= MAXIT-4)
//			if (count == 4)
//			{
//				std::cout << "z0 = " << z0 << ", count = " << count;
//				std::cout << ", z = " << z << ", bright = " << brit;
//				std::cout << ", color = " << FixPoint(z) << ", p = " << p;
//				std::cout << ", er = " << er << std::endl;
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
	glutInitWindowSize(1000,1000);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
