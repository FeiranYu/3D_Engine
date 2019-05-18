#define _CRT_SECURE_NO_WARNINGS
#include<fstream>
#include<iostream>
#include<math.h>
#include<Windows.h>
#include<time.h>
#include<iostream>
#include<tchar.h>

using namespace std;

const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

double angle = 0;

double CamerPosX = 0;
double CamerPosY = 0;
double CamerPosZ = 6;

class color
{
public:
	color(int R, int G, int B) :r(R), g(G), b(B) {}
	color() :r(0), g(0), b(0) {}
	int r, g, b;
};

class vector
{
public:
	vector(double X, double Y, double Z,double W=1) :x(X), y(Y), z(Z), w(W) {}
	vector() :x(0), y(0), z(0), w(1) {}
	double x, y, z, w;

};

class point
{
public:
	point(color COLOR, vector VECTOR) :Color(COLOR), Vector(VECTOR) {}
	point() {};
	color Color;
	vector Vector;
	double u, v;
	bool hide = 0;
};

class matrix
{
public:
	double m[4][4];
};

struct TriangleIndex
{
	int a, b, c;
	double ua, va, ub, vb, uc, vc;
};


double interp(double y, double y0, double y1, double x0, double x1)
{
	if (y < y0)return x0;
	if (y > y1)return x1;
	if (y0 == y1)return x0;
	return (y-y0) / (y1 - y0)*(x1 - x0) + x0;
}

double vecLength(const vector &v)
{
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

void vecNormal(vector &v)
{
	double length = vecLength(v);
	if (length != 0.0f) {
		double inv = 1.0f / length;
		v.x *= inv;
		v.y *= inv;
		v.z *= inv;
	}
}

void vecMul(const vector& a, const vector& b, vector& c)
{
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	c.w = 1;
}

void vecSub(const vector &a, const vector &b, vector& c)
{
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	c.w = a.w - b.w;
}

double  vecDotMul(const vector& a, const vector&b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w;
}


void matMul(const matrix& a, const matrix& b, matrix& c)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			c.m[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				c.m[i][j] += a.m[i][k] * b.m[k][j];
			}
		}
	}
}

void showMat(const matrix& mat)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << mat.m[i][j] << "\t";
		}
		cout << endl;
	}
}

void Swap(double &a, double &b)
{
	double temp = a;
	a = b;
	b = temp;
}

matrix worldMat1;
matrix worldMat2;
matrix camerMat;
matrix projectMat;
matrix finalMat;
matrix tm1, tm2;	//矩阵相乘时的临时矩阵

//位移矩阵
void initWorldMat1(int x, int y, int z)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			worldMat1.m[i][j] = 0;

	for (int i = 0; i < 4; i++)
	{
		worldMat1.m[i][i] = 1;
	}
	worldMat1.m[3][0] = x;
	worldMat1.m[3][1] = y;
	worldMat1.m[3][2] = z;
}


//旋转矩阵
void initWorldMat2(int axis, double angle)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			worldMat2.m[i][j] = 0;
	if (axis == 0)	// x轴
	{
		worldMat2.m[0][0] = 1;
		worldMat2.m[1][1] = cos(angle);
		worldMat2.m[1][2] = sin(angle);
		worldMat2.m[2][1] = -sin(angle);
		worldMat2.m[2][2] = cos(angle);
		worldMat2.m[3][3] = 1;
	}
	if (axis == 1)		// y轴
	{
		worldMat2.m[0][0] = cos(angle);
		worldMat2.m[0][2] = -sin(angle);
		worldMat2.m[1][1] = 1;
		worldMat2.m[2][0] = sin(angle);
		worldMat2.m[2][2] = cos(angle);
		worldMat2.m[3][3] = 1;
	}
	if (axis == 2)			// z轴
	{
		worldMat2.m[0][0] = cos(angle);
		worldMat2.m[0][1] = sin(angle);
		worldMat2.m[1][0] = -sin(angle);
		worldMat2.m[1][1] = cos(angle);
		worldMat2.m[3][3] = 1;
	}
}

void initCamerMat(vector const up, vector const at, vector const eye)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			camerMat.m[i][j] = 0;

	vector xaxis, yaxis, zaxis;

	vecSub(at, eye, zaxis);
	vecNormal(zaxis);
	vecMul(up, zaxis, xaxis);
	vecNormal(xaxis);
	vecMul(zaxis, xaxis, yaxis);

	camerMat.m[0][0] = xaxis.x;
	camerMat.m[1][0] = xaxis.y;
	camerMat.m[2][0] = xaxis.z;
	camerMat.m[3][0] = -vecDotMul(xaxis, eye);

	camerMat.m[0][1] = yaxis.x;
	camerMat.m[1][1] = yaxis.y;
	camerMat.m[2][1] = yaxis.z;
	camerMat.m[3][1] = -vecDotMul(yaxis, eye);

	camerMat.m[0][2] = zaxis.x;
	camerMat.m[1][2] = zaxis.y;
	camerMat.m[2][2] = zaxis.z;
	camerMat.m[3][2] = -vecDotMul(zaxis, eye);

	camerMat.m[3][3] = 1;
	
}

void initProjectMat(float fov, float aspect, float zn, float zf)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			projectMat.m[i][j] = 0;
	projectMat.m[0][0] = 1 / (tan(fov*0.5f)*aspect);
	projectMat.m[1][1] = 1 / tan(fov*0.5f);
	projectMat.m[2][2] = zf / (zf - zn);
	projectMat.m[2][3] = 1.0f;
	projectMat.m[3][2] = (zn*zf) / (zn - zf);
	
}

void CalFinalMat()
{
	initWorldMat1(0, 2, 0);
	initWorldMat2(0, angle);
	vector up;
	vector direction;
	vector pos;
	up.x = 0; up.y = 1; up.z = 0; up.w = 1;
	direction.x = 1; direction.y = 0; direction.z = -1; direction.w = 1;
	pos.x = CamerPosX; pos.y = CamerPosY; pos.z = CamerPosZ; pos.w = 1;
	initCamerMat(up, direction, pos);
	double aspect = (double)SCREEN_WIDTH / (double)SCREEN_HEIGHT;
	initProjectMat(3.1415926/2,aspect, 1, 500);


	matMul(worldMat2, worldMat1, tm1);
	matMul(tm1, camerMat, tm2);
	matMul(tm2, projectMat, finalMat);
}

void Transform(const vector & rawVec, const matrix & mat, vector & finalVec)
{
	finalVec.x = 0;
	finalVec.y = 0;
	finalVec.z = 0;
	finalVec.w = 0;

	finalVec.x = rawVec.x * mat.m[0][0] + rawVec.y * mat.m[1][0] + rawVec.z * mat.m[2][0] + rawVec.w * mat.m[3][0];
	finalVec.y = rawVec.x * mat.m[0][1] + rawVec.y * mat.m[1][1] + rawVec.z * mat.m[2][1] + rawVec.w * mat.m[3][1];
	finalVec.z = rawVec.x * mat.m[0][2] + rawVec.y * mat.m[1][2] + rawVec.z * mat.m[2][2] + rawVec.w * mat.m[3][2];
	finalVec.w = rawVec.x * mat.m[0][3] + rawVec.y * mat.m[1][3] + rawVec.z * mat.m[2][3] + rawVec.w * mat.m[3][3];
}


// 世界坐标到相机坐标
void transform_homogenize(vector & y, vector & x) {
	float rhw = 1.0f / x.w;
	y.x = (x.x * rhw + 1.0f) * SCREEN_WIDTH * 0.5f;
	y.y = (1.0f - x.y * rhw) * SCREEN_HEIGHT * 0.5f;
	y.z = x.z * rhw;
	y.w = x.w;
}

const int point_sum = 8;
const int index_sum = 12;

vector vec_mesh[8];
color c_mesh[8];
point FinalMesh[8];
vector finalV[point_sum];


TriangleIndex index[12] = { 
{0,1,2,0,0,0,1,1,1},
{2,3,0,1,1,1,0,0,0},
{0,1,4,0,0,0,1,1,0},
{5,1,4,1,1,0,1,1,0},
{1,2,5,0,0,0,1,1,0},
{6,2,5,1,1,0,1,1,0},
{3,2,7,0,0,0,1,1,0},
{6,2,7,1,1,0,1,1,0},
{3,0,7,0,0,0,1,1,0},
{4,0,7,1,1,0,1,1,0},
{4,5,7,1,1,1,0,0,1},
{6,5,7,0,0,1,0,0,1}
};


void initMesh()
{
	vec_mesh[0] = vector(1, 1, 1);
	vec_mesh[1] = vector(1, -1, 1);
	vec_mesh[2] = vector(-1, -1, 1);
	vec_mesh[3] = vector(-1, 1, 1);
	vec_mesh[4] = vector(1, 1, -1);
	vec_mesh[5] = vector(1, -1, -1);
	vec_mesh[6] = vector(-1, -1, -1);
	vec_mesh[7] = vector(-1, 1, -1);
	for (int i = 0; i < 8; i++)
	{
		c_mesh[i] = color(255, 255, 255);
		FinalMesh[i].Color = c_mesh[i];
	}
	FinalMesh[0].u = 0;
	FinalMesh[0].v = 0;
	FinalMesh[1].u = 0;
	FinalMesh[1].v = 1;
	FinalMesh[2].u = 1;
	FinalMesh[2].v = 1;
	FinalMesh[3].u = 1;
	FinalMesh[3].v = 0;
}

bool CVVCheck(const vector &vec)
{
	float w = vec.w;
	if (vec.x<-w || vec.x>w)
		return true;
	if (vec.y<-w || vec.y>w)
		return true;
	if (vec.z<0.0f || vec.z>w)
		return true;
	return false;

}


void CalTransform()
{
	
	for (int i = 0; i < point_sum; i++)
	{
		vector finalVV;
		vector tempV;
		Transform(vec_mesh[i], tm2, tempV);
		Transform(vec_mesh[i], finalMat, finalV[i]);
		transform_homogenize(finalVV, finalV[i]);
		//cout << "v:" << vec_mesh[i].x << " " << vec_mesh[i].y << " " << vec_mesh[i].z<<" "<<vec_mesh[i].w << endl;
		//cout << "tempV: " << tempV.x << " " << tempV.y << " " << tempV.z<<" "<<tempV.w << endl;
		//cout << "finalV: " << finalV[i].x << " " << finalV[i].y << " " << finalV[i].z << " " << finalV[i].w << endl;
		//cout << "finalVV: " << finalVV.x << " " << finalVV.y << " " << finalVV.z << " " << finalVV.w << endl;
		FinalMesh[i].hide = CVVCheck(finalV[i]);
		FinalMesh[i].Vector = finalVV;
	}

}

const int bits = 24;
const int MAX_DEPTH = -9999;

BYTE buffer[SCREEN_WIDTH*SCREEN_HEIGHT*bits / 8];
double depth[SCREEN_WIDTH*SCREEN_HEIGHT];

HDC screen_hdc;
HDC hCompatibleDC;
HBITMAP hCompatibleBitmap;
HBITMAP hOldBitmap;
BITMAPINFO binfo;

HINSTANCE hInstance;
WNDCLASS Draw;
HWND hwnd;
MSG msg;

//=========================================================
//	读取图片
unsigned char *pBmpBuf;//读入图像数据的指针
int bmpWidth;//图像的宽
int bmpHeight;//图像的高
RGBQUAD *pColorTable;//颜色表指针
int biBitCount;//图像类型，每像素位数


color loadTexture(double u, double v)
{
	int y = v * bmpHeight;
	int x = u * bmpWidth;
	return color(pBmpBuf[y*bmpWidth * 3 + x * 3 + 2], pBmpBuf[y*bmpWidth * 3 + x * 3 + 1], pBmpBuf[y*bmpWidth * 3 + x * 3]);
}

//-----------------------------------------------------------------------------------------
//给定一个图像位图数据、宽、高、颜色表指针及每像素所占的位数等信息,将其写到指定文件中
bool readBmp(char *bmpName)
{
	FILE *fp = fopen(bmpName, "rb");//二进制读方式打开指定的图像文件
	if (fp == 0)
		return 0;

	//跳过位图文件头结构BITMAPFILEHEADER
	fseek(fp, sizeof(BITMAPFILEHEADER), 0);
	/*
	BITMAPFILEHEADER filehead;
	fread(&filehead, 1, sizeof(BITMAPFILEHEADER), fp);
	showBmpHead(filehead);//显示文件头
*/

//定义位图信息头结构变量，读取位图信息头进内存，存放在变量head中
	BITMAPINFOHEADER infohead;
	fread(&infohead, sizeof(BITMAPINFOHEADER), 1, fp); //获取图像宽、高、每像素所占位数等信息
	bmpWidth = infohead.biWidth;
	bmpHeight = infohead.biHeight;
	biBitCount = infohead.biBitCount;//定义变量，计算图像每行像素所占的字节数（必须是4的倍数）
	//showBmpInforHead(infohead);//显示信息头 


	int lineByte = (bmpWidth * biBitCount / 8 + 3) / 4 * 4;//灰度图像有颜色表，且颜色表表项为256
	if (biBitCount == 8)
	{
		//申请颜色表所需要的空间，读颜色表进内存
		pColorTable = new RGBQUAD[256];
		fread(pColorTable, sizeof(RGBQUAD), 256, fp);
	}

	//申请位图数据所需要的空间，读位图数据进内存
	pBmpBuf = new unsigned char[lineByte * bmpHeight];
	fread(pBmpBuf, 1, lineByte * bmpHeight, fp);
	fclose(fp);//关闭文件
	return 1;//读取文件成功
}


//=========================================================


void CleanBuffer()
{
	for (int y = 0; y < SCREEN_HEIGHT; y++)
	{
		for (int x = 0; x < SCREEN_WIDTH; x++)
		{
			buffer[y*SCREEN_WIDTH * 3 + (x + 1) * 3 - 1] = 255-255* (double)y/ (double)SCREEN_HEIGHT;
			buffer[y*SCREEN_WIDTH * 3 + (x + 1) * 3 - 2] = 255-255 * (double)y / (double)SCREEN_HEIGHT;
			buffer[y*SCREEN_WIDTH * 3 + (x + 1) * 3 - 3] = 255-255 * (double)y / (double)SCREEN_HEIGHT;
			depth[y*SCREEN_WIDTH + x] = 20;
		}
	}
}

void DrawPoint(const point &p)
{
	//cout << "Point z " << p.Vector.z << " depth " << depth[(int)p.Vector.y*SCREEN_WIDTH + (int)p.Vector.x] << endl;
	if (p.Vector.x > 0 && p.Vector.x < SCREEN_WIDTH&&p.Vector.y>0 && p.Vector.y < SCREEN_HEIGHT&&p.Vector.z<depth[(int)p.Vector.y*SCREEN_WIDTH+(int)p.Vector.x])
	{
		depth[(int)p.Vector.y*SCREEN_WIDTH + (int)p.Vector.x] = p.Vector.z;
		buffer[(int)p.Vector.y*SCREEN_WIDTH * 3 + ((int)p.Vector.x + 1) * 3 - 1] = p.Color.r;
		buffer[(int)p.Vector.y*SCREEN_WIDTH * 3 + ((int)p.Vector.x + 1) * 3 - 2] = p.Color.g;
		buffer[(int)p.Vector.y*SCREEN_WIDTH * 3 + ((int)p.Vector.x + 1) * 3 - 3] = p.Color.b;
	}
}


void DrawLine(const point &leftp,const point &rightp)
{
	for (int i = leftp.Vector.x; i < rightp.Vector.x; i++)
	{
		point p;
		p.Vector = vector(i, leftp.Vector.y, leftp.Vector.z);
		p.Vector.z =1/ interp(i, leftp.Vector.x, rightp.Vector.x, 1/leftp.Vector.z,1/ rightp.Vector.z);
		//cout << "leftz " << leftp.Vector.z << " rightz " << rightp.Vector.z << " p z " << p.Vector.z << endl;
		p.u = interp(i, leftp.Vector.x, rightp.Vector.x, leftp.u/leftp.Vector.z, rightp.u/rightp.Vector.z)*p.Vector.z;
		p.v = interp(i, leftp.Vector.x, rightp.Vector.x, leftp.v/leftp.Vector.z, rightp.v/rightp.Vector.z)*p.Vector.z;
		//cout << p.Vector.z << endl;

		//if(i==leftp.Vector.x+20)
			//cout << "u " << p.u << " v " << p.v << endl;
		p.Color = loadTexture(p.u, p.v);
		DrawPoint (p);
	}
}

void DrawTriangle(const point &p1, const point &p2, const point &p3, color Color=color(255,255,255))
{
	point P1, P2, P3;
	P1 = p1; P2 = p2; P3 = p3;

	//cout << "P1 u " << P1.u << " v " << P1.v << " P2 u " << P2.u << " v " << P2.v << " P3 u " << P3.u << " v " << P3.v << endl;


	point temp;
	if (P1.hide == 1 || P2.hide == 1 || P3.hide == 1)
		return;

	if (P2.Vector.y < P1.Vector.y)
	{
		temp = P2;
		P2 = P1;
		P1 = temp;
	}
	if (P3.Vector.y < P1.Vector.y)
	{
		temp = P3;
		P3 = P1;
		P1 = temp;
	}
	if (P2.Vector.y > P3.Vector.y)
	{
		temp = P2;
		P2 = P3;
		P3 = temp;
	}
	// y: P1<P2<P3

	if (P1.Vector.y == P2.Vector.y)
	{
		for (int y = P1.Vector.y; y < P3.Vector.y; y++)
		{
			int left, right;
			left=interp(y, P1.Vector.y, P3.Vector.y, P1.Vector.x, P3.Vector.x);
			right = interp(y, P2.Vector.y, P3.Vector.y, P2.Vector.x, P3.Vector.x);

			double zleft_t = interp(y, P1.Vector.y, P3.Vector.y, 1 / P1.Vector.w, 1 / P3.Vector.w);
			double zright_t = interp(y, P2.Vector.y, P3.Vector.y, 1 / P2.Vector.w, 1 / P3.Vector.w);

			double uleft = interp(y, P1.Vector.y, P3.Vector.y, P1.u / P1.Vector.w, P3.u / P3.Vector.w)*(1/zleft_t);
			double vleft = interp(y, P1.Vector.y, P3.Vector.y, P1.v / P1.Vector.w, P3.v / P3.Vector.w)*(1/zleft_t);
			double uright = interp(y, P2.Vector.y, P3.Vector.y, P2.u/ P2.Vector.w, P3.u / P3.Vector.w)*(1/zright_t);
			double vright = interp(y, P2.Vector.y, P3.Vector.y, P2.v/ P2.Vector.w, P3.v / P3.Vector.w)*(1/zright_t);



			if (left > right)
			{
				int temp = left;
				left = right;
				right = temp;
				double ztemp = zleft_t;
				zleft_t = zright_t;
				zright_t = ztemp;
				double utemp = uleft;
				uleft = uright;
				uright = utemp;
				double vtemp = vleft;
				vleft = vright;
				vright = vtemp;
			}
			point pleft, pright;
			pleft.Vector = vector(left, y, 1/zleft_t);
			pleft.Color = Color;
			pright.Vector = vector(right, y, 1/zright_t);
			pright.Color = Color;
			//cout << "P1 z" << P1.Vector.z << "P2 z " << P2.Vector.z << " P3 z " << P3.Vector.z << " zleft " << pleft.Vector.z << endl;
			pleft.u = uleft;
			pleft.v = vleft;
			pright.u = uright;
			pright.v = vright;

			DrawLine(pleft, pright);

		}
	}
	else if (P2.Vector.y == P3.Vector.y)
	{
		for (int y = P1.Vector.y; y < P2.Vector.y; y++)
		{
			int left, right;
			left = interp(y, P1.Vector.y, P2.Vector.y, P1.Vector.x, P2.Vector.x);
			right = interp(y, P1.Vector.y, P3.Vector.y, P1.Vector.x, P3.Vector.x);
			double zleft_t = interp(y, P1.Vector.y, P2.Vector.y, 1/P1.Vector.w, 1/P2.Vector.w);
			double zright_t = interp(y, P1.Vector.y, P3.Vector.y, 1/P1.Vector.w, 1/P3.Vector.w);

			double uleft = interp(y, P1.Vector.y, P2.Vector.y, P1.u/P1.Vector.w, P2.u/P2.Vector.w)*(1/zleft_t);
			double vleft = interp(y, P1.Vector.y, P2.Vector.y, P1.v/P1.Vector.w, P2.v/P2.Vector.w)*(1/zleft_t);
			double uright = interp(y, P1.Vector.y, P3.Vector.y, P1.u/P1.Vector.w, P3.u/P3.Vector.w)*(1/zright_t);
			double vright = interp(y, P1.Vector.y, P3.Vector.y, P1.v/P1.Vector.w, P3.v/P3.Vector.w)*(1/zright_t);
			if (left > right)
			{
				int temp = left;
				left = right;
				right = temp;
				double ztemp = zleft_t;
				zleft_t = zright_t;
				zright_t = ztemp;
				double utemp = uleft;
				uleft = uright;
				uright = utemp;
				double vtemp = vleft;
				vleft = vright;
				vright = vtemp;
			}
			point pleft, pright;
			pleft.Vector = vector(left, y, 1/zleft_t);
			pleft.Color = Color;
			pright.Vector = vector(right, y, 1/zright_t);
			pright.Color = Color;
			pleft.u = uleft;
			pleft.v = vleft;
			pright.u = uright;
			pright.v = vright;
			//cout << "P1 z" << P1.Vector.z << "P2 z " << P2.Vector.z << " P3 z " << P3.Vector.z << " zleft " << pleft.Vector.z << endl;

			DrawLine(pleft, pright);
		}
	}
	else
	{
		for (int y = P1.Vector.y+1; y < P2.Vector.y; y++)
		{
			int left, right;
			left = interp(y, P1.Vector.y, P2.Vector.y, P1.Vector.x, P2.Vector.x);
			right = interp(y, P1.Vector.y, P3.Vector.y, P1.Vector.x, P3.Vector.x);
			double zleft_t = interp(y, P1.Vector.y, P2.Vector.y, 1/P1.Vector.w, 1/P2.Vector.w);
			double zright_t = interp(y, P1.Vector.y, P3.Vector.y,1/P1.Vector.w, 1/P3.Vector.w);


			double uleft = interp(y, P1.Vector.y, P2.Vector.y, P1.u / P1.Vector.w, P2.u / P2.Vector.w)*(1 / zleft_t);
			double vleft = interp(y, P1.Vector.y, P2.Vector.y, P1.v / P1.Vector.w, P2.v / P2.Vector.w)*(1 / zleft_t);
			double uright = interp(y, P1.Vector.y, P3.Vector.y, P1.u / P1.Vector.w, P3.u / P3.Vector.w)*(1 / zright_t);
			double vright = interp(y, P1.Vector.y, P3.Vector.y, P1.v / P1.Vector.w, P3.v / P3.Vector.w)*(1 / zright_t);
			if (left > right)
			{
				int temp = left;
				left = right;
				right = temp;
				double ztemp = zleft_t;
				zleft_t = zright_t;
				zright_t = ztemp;
				double utemp = uleft;
				uleft = uright;
				uright = utemp;
				double vtemp = vleft;
				vleft = vright;
				vright = vtemp;
			}
			point pleft, pright;
			pleft.Vector = vector(left, y, 1/zleft_t);
			pleft.Color = Color;
			pright.Vector = vector(right, y, 1/zright_t);
			pright.Color = Color;

			pleft.u = uleft;
			pleft.v = vleft;
			pright.u = uright;
			pright.v = vright;
			//cout << "P1 z" << P1.Vector.z << "P2 z " << P2.Vector.z << " P3 z " << P3.Vector.z << " zleft " << pleft.Vector.z << endl;

			DrawLine(pleft,pright);
		}
		for (int y = P2.Vector.y+1; y < P3.Vector.y; y++)
		{
			int left, right;
			left = interp(y, P2.Vector.y, P3.Vector.y, P2.Vector.x, P3.Vector.x);
			right = interp(y, P1.Vector.y, P3.Vector.y, P1.Vector.x, P3.Vector.x);
			double zleft_t = interp(y, P2.Vector.y, P3.Vector.y, 1/P2.Vector.w, 1/P3.Vector.w);
			double zright_t = interp(y, P1.Vector.y, P3.Vector.y, 1/P1.Vector.w, 1/P3.Vector.w);

			double uleft = interp(y, P2.Vector.y, P3.Vector.y, P2.u / P2.Vector.w, P3.u / P3.Vector.w)*(1 / zleft_t);
			double vleft = interp(y, P2.Vector.y, P3.Vector.y, P2.v / P2.Vector.w, P3.v / P3.Vector.w)*(1 / zleft_t);
			double uright = interp(y, P1.Vector.y, P3.Vector.y, P1.u / P1.Vector.w, P3.u / P3.Vector.w)*(1 / zright_t);
			double vright = interp(y, P1.Vector.y, P3.Vector.y, P1.v / P1.Vector.w, P3.v / P3.Vector.w)*(1 / zright_t);
			if (left > right)
			{
				int temp = left;
				left = right;
				right = temp;
				double ztemp = zleft_t;
				zleft_t = zright_t;
				zright_t = ztemp;
				double utemp = uleft;
				uleft = uright;
				uright = utemp;
				double vtemp = vleft;
				vleft = vright;
				vright = vtemp;
			}
			point pleft, pright;
			pleft.Vector = vector(left, y, 1/zleft_t);
			pleft.Color = Color;
			pright.Vector = vector(right, y, 1/zright_t);
			pright.Color = Color;

			pleft.u = uleft;
			pleft.v = vleft;
			pright.u = uright;
			pright.v = vright;
			//cout<<"P1 z"<<P1.Vector.z << "P2 z " << P2.Vector.z << " P3 z " << P3.Vector.z << " zleft " << pleft.Vector.z << endl;

			DrawLine(pleft, pright);
		}
	}
}

LRESULT CALLBACK WindowProc(
	_In_	HWND hwnd,
	_In_	UINT uMsg,
	_In_	WPARAM wParam,
	_In_	LPARAM lParam
)
{
	switch (uMsg)
	{
		case WM_DESTROY:
		{
			PostQuitMessage(0);
			return 0;
		}
		case WM_KEYDOWN:
		{
			if (wParam == VK_ESCAPE)
				exit(0);
			if (wParam == VK_LEFT)
				angle += 0.05;
			if (wParam == VK_RIGHT)
				angle -= 0.05;
			if (wParam == VK_UP);
			
			if (wParam == VK_DOWN);
			
			if (wParam == 'W')
				CamerPosZ -= 0.05;
			if (wParam == 'S')
				CamerPosZ += 0.05;
			if (wParam == 'A')
				CamerPosX += 0.05;
			if (wParam == 'D')
				CamerPosX -= 0.05;
			if (wParam == 'Q')
				CamerPosY -= 0.05;
			if (wParam == 'E')
				CamerPosY += 0.05;
		}
	}
	return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

void PutBufferToScreen()
{
	SetDIBits(screen_hdc, hCompatibleBitmap, 0, SCREEN_HEIGHT, buffer, (BITMAPINFO*)&binfo, DIB_RGB_COLORS);
	BitBlt(screen_hdc, -1, -1, SCREEN_WIDTH, SCREEN_HEIGHT, hCompatibleDC, 0, 0, SRCCOPY);
}

void GameLoop()
{
	CleanBuffer();
	CalFinalMat();
	CalTransform();

	/*
	for (int y = (SCREEN_HEIGHT-bmpHeight)/2; y < bmpHeight+ (SCREEN_HEIGHT - bmpHeight) / 2; y++)
	{
		for (int x = (SCREEN_WIDTH-bmpWidth)/2; x < bmpWidth+ (SCREEN_WIDTH - bmpWidth) / 2; x++)
		{
			point p;
			p.Vector = vector(x, bmpHeight- y, 10);
			p.Color = color(pBmpBuf[(y - (SCREEN_HEIGHT - bmpHeight) / 2) * 3 * bmpWidth + (x-(SCREEN_WIDTH-bmpWidth)/2) * 3 + 2],
							pBmpBuf[(y - (SCREEN_HEIGHT - bmpHeight) / 2) * 3 * bmpWidth + (x-(SCREEN_WIDTH - bmpWidth)/2) * 3 + 1], 
							pBmpBuf[(y - (SCREEN_HEIGHT - bmpHeight) / 2) * 3 * bmpWidth + (x- (SCREEN_WIDTH - bmpWidth)/2) * 3]);
			DrawPoint(p);
		}
	}
	*/
	for (int i = 0; i < index_sum; i++)
	{
		int a = index[i].a;
		int b = index[i].b;
		int c = index[i].c;
		FinalMesh[a].u = index[i].ua;
		FinalMesh[a].v = index[i].va;
		FinalMesh[b].u = index[i].ub;
		FinalMesh[b].v = index[i].vb;
		FinalMesh[c].u = index[i].uc;
		FinalMesh[c].v = index[i].vc;
	
		DrawTriangle(FinalMesh[a], FinalMesh[b], FinalMesh[c]);
	}
	PutBufferToScreen();
}

void initWindow()
{
	hInstance = GetModuleHandle(NULL);

	Draw.cbClsExtra = 0;
	Draw.cbWndExtra = 0;
	Draw.hCursor = LoadCursor(hInstance, IDC_ARROW);
	Draw.hIcon = LoadIcon(hInstance, IDI_APPLICATION);
	Draw.lpszMenuName = NULL;
	Draw.style = WS_MINIMIZEBOX | CS_HREDRAW | CS_VREDRAW;
	Draw.hbrBackground = (HBRUSH)COLOR_WINDOW;
	Draw.lpfnWndProc = WindowProc;
	Draw.lpszClassName = _T("DDraw");
	Draw.hInstance = hInstance;

	RegisterClass(&Draw);

	hwnd = CreateWindow(
		_T("DDraw"),
		L"Draw",
		WS_OVERLAPPEDWINDOW,
		38,
		20,
		SCREEN_WIDTH + 15,
		SCREEN_HEIGHT + 38,
		NULL,
		NULL,
		hInstance,
		NULL
	);

	ShowWindow(hwnd, SW_SHOW);
	UpdateWindow(hwnd);

	//init bitbuffer
	ZeroMemory(&binfo, sizeof(BITMAPINFO));
	binfo.bmiHeader.biBitCount = bits;
	binfo.bmiHeader.biCompression = BI_RGB;
	binfo.bmiHeader.biHeight = -SCREEN_HEIGHT;
	binfo.bmiHeader.biPlanes = 1;
	binfo.bmiHeader.biSizeImage = 0;
	binfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	binfo.bmiHeader.biWidth = SCREEN_WIDTH;

	screen_hdc = GetDC(hwnd);
	hCompatibleDC = CreateCompatibleDC(screen_hdc);
	hCompatibleBitmap = CreateCompatibleBitmap(screen_hdc, SCREEN_WIDTH, SCREEN_HEIGHT);
	hOldBitmap = (HBITMAP)SelectObject(hCompatibleDC, hCompatibleBitmap);

	while (1)
	{	
		if (PeekMessage(&msg, NULL, 0, 0,PM_REMOVE))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		GameLoop();
	}
}

int main()
{
	srand(time(0));
	char readPath[] = "1.bmp";
	if (readBmp(readPath))
	{
		cout << "readOK!" << endl;
		cout << "\nwidth=" << bmpWidth << "\nheight=" << bmpHeight << endl;
	}
	initMesh();
	initWindow();
	return 0;
}