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

#include <thrust/complex.h>

#define EPS 0.000001  // 停止判定
#define MAXIT 30      // 最大反復回数
#define ZMAX 4.0      // 初期値の最大絶対値
#define ZOOM 200      // 拡大率
#define RMAX 2000     // 複素平面の分割数
#define ORD  32       // 次数

#if defined (__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

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

void Reg0()
{
	gl
	glBegin(GL_POLYGON);
	glEnd();
}

void display(void)
{
	//////////////////////////////////////////////
	// 背景を白に
	glClearColor(1.0, 1.0, 1.0, 1.0); // 塗りつぶしの色を指定
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	setZero(fps);     // 零点のセット

	// 領域 0 の描画（ \pm4 + \pm 4 i）の枠
	Reg0();

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
	glOrtho( -(double)w/ZOOM, (double)w/ZOOM, -(double)h/ZOOM, (double)h/ZOOM, -1.0, 1.0);
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
