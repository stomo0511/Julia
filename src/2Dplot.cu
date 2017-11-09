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

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	glColor3d(1.0, 0.0, 0.0);         // 赤 (1,0,0) で描画

	///////////////////////////////////
	// 2点間に線を描画
//	glBegin(GL_LINES);
//
//	glVertex2d(-3.0,  0.0);
//	glVertex2d( 3.0,  0.0);
//
//	glVertex2d( 0.0, -3.0);
//	glVertex2d( 0.0,  3.0);
//
//	glEnd();

	thrust::complex<double> p = thrust::complex<double>(1.0, 0.0);

	int points = 4;
	glBegin(GL_LINE_LOOP);
	for (int i=0; i<points; i++)
	{
		glVertex2d(p.real(),p.imag());
		p += thrust::complex<double>(0.2, 0.0);
	}
	glEnd();
	///////////////////////////////////

	glFlush();                        // OpenGL命令のフラッシュ
}

void resize(int w, int h)
{
	// Window全体をView portにする
	glViewport(0,0,w,h);

	// 変換行列の初期化
	glLoadIdentity();

	// Screen上の表示領域をView portの大きさに比例させる
	glOrtho( -w/250.0, w/250.0, -h/250.0, h/250.0, -1.0, 1.0);
}

void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 1.0); // 塗りつぶしの色を指定
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);          // OpenGL初期化
	glutInitDisplayMode(GLUT_RGBA); // RGBモードに設定
	glutInitWindowSize(1000,1000);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	init();
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
