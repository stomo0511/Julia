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

#define ZMAX 4.0      // 初期値の最大絶対値
#define ZOOM 200      // 拡大率
#define RMAX 2000     // 複素平面の分割数

// Zeros
std::vector< thrust::complex<double> > Zrs {
	thrust::complex<double> (  0.0,  1.0 ),
	thrust::complex<double> (  1.0,  2.0 ),
	thrust::complex<double> ( -1.0,  2.0 ),
	thrust::complex<double> (  3.0, -3.0 ),
	thrust::complex<double> ( -3.0, -3.0 )
};

void PBisec( const int i, const int j, double &a, double &b )
{
	const double a1 = Zrs[i].real();
	const double a2 = Zrs[j].real();
	const double b1 = Zrs[i].imag();
	const double b2 = Zrs[j].imag();

//	std::cout << "a1 = " << a1 << ", b1 = " << b1 << ", a2 = " << a2 << ", b2 = " << b2 << std::endl;

	a = (b2 - b1) / (a2 - a1);
	b = (a2*b1 - a1*b2) / (a2 - a1);

//	std::cout << "a = " << a << ", b = " << b << std::endl;

	a = -1.0 / a;  // 傾き = 2点を通る直線に垂直
	// 2点の中点座標
	double mx = (a1+a2)/2.0;
	double my = (b1+b2)/2.0;
//	std::cout << "mx = " << mx << ", my = " << my << std::endl;

	b = my - a*mx;
}

// 領域0 枠（黒）
void Reg0()
{
	glColor3d(0.0,0.0,0.0);
	glBegin(GL_POLYGON);
	glVertex2d( 4.0, 4.0);
	glVertex2d( 4.0,-4.0);
	glVertex2d(-4.0,-4.0);
	glVertex2d(-4.0, 4.0);
	glEnd();
}

//零点 0 の Voronoi Cell
void Cell0()
{
	glColor3d(1.0,0.0,0.0);
	glBegin(GL_POLYGON);
	glVertex2d( 0.0, 2.0 );
	glVertex2d(  2.3571, -0.3571 );
	glVertex2d( 0.0, -2.125 );
	glVertex2d( -2.3571, -0.3571 );
	glEnd();
}

// 零点 1 の Voronoi Cell
void Cell1()
{
	glColor3d(0.0,1.0,0.0);
	glBegin(GL_POLYGON);
	glVertex2d( 0.0, 2.0 );
	glVertex2d( 2.3571, -0.3571 );
	glVertex2d( 4.0, 0.3 );
	glVertex2d( 4.0, 4.0 );
	glVertex2d( 0.0, 4.0 );
	glEnd();
}

// 零点 2 の Voronoi Cell
void Cell2()
{
	glColor3d(0.0,0.0,1.0);
	glBegin(GL_POLYGON);
	glVertex2d( 0.0, 2.0 );
	glVertex2d( -2.3571, -0.3571 );
	glVertex2d( -4.0, 0.3 );
	glVertex2d( -4.0, 4.0 );
	glVertex2d( 0.0, 4.0 );
	glEnd();
}

// 零点 3 の Voronoi Cell
void Cell3()
{
	glColor3d(1.0,0.0,1.0);
	glBegin(GL_POLYGON);
	glVertex2d( 2.3571, -0.3571 );
	glVertex2d( 0.0, -2.125 );
	glVertex2d( 0.0, -4.0 );
	glVertex2d( 4.0, -4.0 );
	glVertex2d( 4.0,  0.3 );
	glEnd();
}

// 零点 4 の Voronoi Cell
void Cell4()
{
	glColor3d(0.0,1.0,1.0);
	glBegin(GL_POLYGON);
	glVertex2d( -2.3571, -0.3571 );
	glVertex2d( 0.0, -2.125 );
	glVertex2d( 0.0, -4.0 );
	glVertex2d( -4.0, -4.0 );
	glVertex2d( -4.0,  0.3 );
	glEnd();
}

void display(void)
{
	//////////////////////////////////////////////
	// 背景を白に
	glClearColor(1.0, 1.0, 1.0, 1.0); // 塗りつぶしの色を指定
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	// 領域 0 の描画（ \pm4 + \pm 4 i）の枠
//	Reg0();

	//零点 0 の Voronoi Cell
	Cell0();

	//零点 1 の Voronoi Cell
	Cell1();

	//零点 2 の Voronoi Cell
	Cell2();

	//零点 3 の Voronoi Cell
	Cell3();

	//零点 4 の Voronoi Cell
	Cell4();

	// Z_i の描画
	for (auto itr = Zrs.begin(); itr < Zrs.end(); itr++)
	{
		glColor3d(1.0,1.0,1.0);   // 白の点を描画
		glPointSize(8.0);      // 点の大きさ（ディフォルトは1.0)
		glBegin(GL_POINTS);
		glVertex2d( (*itr).real(), (*itr).imag() );
		glEnd();
	}

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
	// 垂直二等分線
	double a, b;
	PBisec( 0, 1, a, b );
	std::cout << "(0,1) の垂直二等分線の傾き : " << a << ", 切片 : " << b << std::endl;

	PBisec( 0, 2, a, b );
	std::cout << "(0,2) の垂直二等分線の傾き : " << a << ", 切片 : " << b << std::endl;

	PBisec( 0, 3, a, b );
	std::cout << "(0,3) の垂直二等分線の傾き : " << a << ", 切片 : " << b << std::endl;

	PBisec( 0, 4, a, b );
	std::cout << "(0,4) の垂直二等分線の傾き : " << a << ", 切片 : " << b << std::endl;

	PBisec( 1, 3, a, b );
	std::cout << "(1,3) の垂直二等分線の傾き : " << a << ", 切片 : " << b << std::endl;

	PBisec( 2, 4, a, b );
	std::cout << "(2,4) の垂直二等分線の傾き : " << a << ", 切片 : " << b << std::endl;

	glutInit(&argc, argv);          // OpenGL初期化
	glutInitWindowSize(1000,1000);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
