#include<iostream>
#include<math.h>

using namespace std;

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 600

class vector
{
public:
	vector(int X, int Y, int Z) :x(X), y(Y), z(Z), w(1) {}
	vector() :x(0), y(0), z(0), w(1) {}
	double x, y, z, w;
};

class matrix
{
public:
	double m[4][4];
};

void vecMul(const vector& a, const vector& b, vector& c)
{
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	c.w = 1;
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

vector mesh[8];
matrix worldMat1;
matrix worldMat2;
matrix camerMat;
matrix projectMat;

matrix finalMat;


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
void initWorldMat2(int axis, int angle)
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

void initCamerMat(vector const upVec, vector const directionVec, vector const posVec)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			camerMat.m[i][j] = 0;
	vector rightVec;
	vecMul(directionVec, upVec, rightVec);
	camerMat.m[0][0] = rightVec.x;
	camerMat.m[1][0] = rightVec.y;
	camerMat.m[2][0] = rightVec.z;
	camerMat.m[0][1] = directionVec.x;
	camerMat.m[1][1] = directionVec.y;
	camerMat.m[2][1] = directionVec.z;
	camerMat.m[0][2] = upVec.x;
	camerMat.m[1][2] = upVec.y;
	camerMat.m[2][2] = upVec.z;
	camerMat.m[3][0] = -posVec.x;
	camerMat.m[3][1] = -posVec.y;
	camerMat.m[3][2] = -posVec.z;
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

	
	/*
	
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			projectMat.m[i][j] = 0;
	projectMat.m[0][0] = 2 * zn * SCREEN_WIDTH / (SCREEN_HEIGHT * SCREEN_HEIGHT);
	projectMat.m[1][1] = 2 * zn / SCREEN_HEIGHT;
	projectMat.m[2][2] = zf / (zf - zn);
	projectMat.m[2][3] = 1;
	projectMat.m[3][2] = (zn * zf) / (zn - zf);
	*/
}


vector raw0;

void CalFinalMat()
{
	initWorldMat1(0, 0, 0);
	initWorldMat2(0, 1.5);
	vector up;
	vector direction;
	vector pos;
	up.x = 0; up.y = 0; up.z = 1; up.w = 0;
	direction.x = 1; direction.y = 0; direction.z = 0; direction.w = 0;
	pos.x = 0; pos.y = 0; pos.z = 5; pos.w = 0;
	initCamerMat(up, direction, pos);
	initProjectMat(3.1, (double)SCREEN_HEIGHT/(double)SCREEN_WIDTH, 1.0f, 1000);
	matrix tm1, tm2;
	matMul(worldMat1, worldMat2, tm1);
	matMul(tm1, camerMat, tm2);
	matMul(tm2, projectMat, finalMat);

	cout << "worldMat1" << endl;
	showMat(worldMat1);
	cout << "worldMat2" << endl;
	showMat(worldMat2);
	cout << "cameraMat" << endl;
	showMat(camerMat);
	cout << "projectMat" << endl;
	showMat(projectMat);
	cout << "finalMat" << endl;
	showMat(finalMat);
}

void Transform(const vector & rawVec, const matrix & mat, vector & finalVec)
{
	finalVec.x = 0;
	finalVec.y = 0;
	finalVec.z = 0;
	finalVec.w = 0;

	finalVec.x = rawVec.x * mat.m[0][0] + rawVec.y * mat.m[0][1] + rawVec.z * mat.m[0][2] + rawVec.w * mat.m[0][3];
	finalVec.y = rawVec.x * mat.m[1][0] + rawVec.y * mat.m[1][1] + rawVec.z * mat.m[1][2] + rawVec.w * mat.m[1][3];
	finalVec.z = rawVec.x * mat.m[2][0] + rawVec.y * mat.m[2][1] + rawVec.z * mat.m[2][2] + rawVec.w * mat.m[2][3];
	finalVec.w = rawVec.x * mat.m[3][0] + rawVec.y * mat.m[3][1] + rawVec.z * mat.m[3][2] + rawVec.w * mat.m[3][3];
}


// 世界坐标到相机坐标
void transform_homogenize(vector & y, vector & x) {
	float rhw = 1.0f / x.w;
	y.x = (x.x * rhw + 1.0f) * SCREEN_WIDTH * 0.5f;
	y.y = (1.0f - x.y * rhw) * SCREEN_HEIGHT * 0.5f;
	y.z = x.z * rhw;
	y.w = 1.0f;
}

int main()
{
	CalFinalMat();
	vector v = vector(-1, 1, -1);
	vector finalV;
	vector finalVV;
	Transform(v, finalMat, finalV);
	transform_homogenize(finalVV, finalV);
	cout << "v:" << v.x << " " << v.y << " " << v.z << endl;
	cout << "finalV:" << finalV.x << " " << finalV.y << " " << finalV.z << " " << finalV.w << endl;
	cout << "finalVV:" << finalVV.x << " " << finalVV.y << " " << finalVV.z << " " << finalVV.w << endl;

	cin.get();
	return 0;
}