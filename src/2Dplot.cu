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
#include <GLUT/glut.h>
#include <thrust/complex.h>

#define EPS 0.0000001  // 停止判定
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

template<typename T> thrust::complex<T> vf( thrust::complex<T> cf, thrust::complex<T> z )
{
	thrust::omplex<T> tmp = cf[0];
    for(auto itr = cf.begin()+1; itr != cf.end(); ++itr)
    {
    	tmp = tmp*z + *itr;
    }
	return tmp;
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
	glPointSize(8.0);      // 点の大きさ（ディフォルトは1.0)
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
	// 背景を白に
	glClearColor(1.0, 1.0, 1.0, 1.0); // 塗りつぶしの色を指定（黒）
	glClear(GL_COLOR_BUFFER_BIT);     // 塗りつぶし

	setZero(fps);      // 零点のセット

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
	glOrtho( -(double)w/ZOOM, (double)w/ZOOM, -(double)h/ZOOM, (double)h/ZOOM, -1.0, 1.0);
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
	glutInitWindowSize(500,500);  // 初期Windowサイズ指定
	glutCreateWindow(argv[0]);      // Windowを開く
	glutDisplayFunc(display);       // Windowに描画
	glutReshapeFunc(resize);
	glutMainLoop();                 // イベント待ち

	return EXIT_SUCCESS;
}
