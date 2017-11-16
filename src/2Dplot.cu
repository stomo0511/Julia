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
#define MAXIT 30      // 最大反復回数
#define ZMAX 0.8      // 初期値の最大絶対値
#define ZOOM 600      // 拡大率
#define RMAX 2000     // 複素平面の分割数
#define XOFS 1.0      // x軸のオフセット

#if defined (__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

#define NFP 4  // 不動点の数
thrust::complex<double> fps[NFP];

void setZero( thrust::complex<double> *fps )
{
	fps[0] = thrust::complex<double> (  1.0,  0.0 );
	fps[1] = thrust::complex<double> (  0.0,  0.5 );
	fps[2] = thrust::complex<double> (  0.0, -0.5 );
	fps[3] = thrust::complex<double> ( -1.0,  0.0 );
}

template<typename T> thrust::complex<T> vf( thrust::complex<T> z )
{
	return z*z*z*z -(3.0*z*z)/(4.0) -(1.0)/(4.0);
}

template<typename T> thrust::complex<T> df( thrust::complex<T> z )
{
	return 4.0*z*z*z - (3.0*z)/(2.0);
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
	int col = 0;
	double min = (double)(MAXIT);

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

void DrawApollonius( int i, int j, double alp )
{
	const int pts = 180;    // 円周上の点数

	thrust::complex<double> center = (fps[i] - alp*alp*fps[j]) / (1.0 - alp*alp);
	double radius = alp*abs(fps[i] - fps[j]) / (1.0 - alp*alp);
	double tic = (double)(2.0*M_PI / pts);

	// Z_i の描画
	glColor3d(1.0,1.0,1.0);   // 白の点を描画
	glPointSize(10.0);      // 点の大きさ（ディフォルトは1.0)
	glBegin(GL_POINTS);
	glVertex2d( fps[i].real() -XOFS, fps[i].imag() );
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
	// 背景を黒に
	glClearColor(0.0, 0.0, 0.0, 1.0); // 塗りつぶしの色を指定（黒）
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	setZero(fps);      // 零点のセット

	//////////////////////////////////////////////
	// 点の描画
	glBegin(GL_POINTS);

//	double x = (double)(-ZMAX);
//	double x = (double)0.0;
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

			double brit;
			if (count > 13)
				brit = 0.0;
			else
				brit = (13.0 - double(count)) / 13.0;
			// 明るさの補正
			brit += 0.15;

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
	double alp = (double)1.0 / (2*NFP - 3.0);

	for (int j=1; j<NFP; j++)
	{
		DrawApollonius( 0, j, alp );
	}
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

void saveImage( const unsigned int imageWidth, const unsigned int imageHeight )
{
	const unsigned int channelNum = 3; // RGBなら3, RGBAなら4
	void* dataBuffer = ( GLubyte* )malloc( imageWidth * imageHeight * channelNum );

	// 読み取るOpneGLのバッファを指定 GL_FRONT:フロントバッファ　GL_BACK:バックバッファ
	glReadBuffer( GL_FRONT );

	// OpenGLで画面に描画されている内容をバッファに格納
	glReadPixels(
			0,                 //読み取る領域の左下隅のx座標
			0,                 //読み取る領域の左下隅のy座標 //0 or getCurrentWidth() - 1
			imageWidth,             //読み取る領域の幅
			imageHeight,            //読み取る領域の高さ
			GL_RGB, //it means GL_BGR,           //取得したい色情報の形式
			GL_UNSIGNED_BYTE,  //読み取ったデータを保存する配列の型
			dataBuffer      //ビットマップのピクセルデータ（実際にはバイト配列）へのポインタ
	);
	glFlush();

//	GLubyte* p = static_cast<GLubyte*>( dataBuffer );
//	std::string fname = "outputImage.jpg";
//	IplImage* outImage = cvCreateImage( cvSize( imageWidth, imageHeight ), IPL_DEPTH_8U, 3 );
//
//	for ( unsigned int j = 0; j < imageHeight; ++ j )
//	{
//		for ( unsigned int i = 0; i < imageWidth; ++i )
//		{
//			outImage->imageData[ ( imageHeight - j - 1 ) * outImage->widthStep + i * 3 + 0 ] = *p;
//			outImage->imageData[ ( imageHeight - j - 1 ) * outImage->widthStep + i * 3 + 1 ] = *( p + 1 );
//			outImage->imageData[ ( imageHeight - j - 1 ) * outImage->widthStep + i * 3 + 2 ] = *( p + 2 );
//			p += 3;
//		}
//	}
//
//	cvSaveImage( fname.c_str(), outImage );
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
