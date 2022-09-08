#include "olcConsoleGameEngine.h"
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

// 2d structure holding texture coordinates
struct vec2d {
	float u = 0;
	float v = 0;
	float w = 1;
};

struct vec3d {
	float x = 0;
	float y = 0;
	float z = 0;
	float w = 1.0f; // required for various matrix vector multiplication
};

struct triangle {
	vec3d p[3];
	vec2d t[3];
	// represent the colours of the triangle
	wchar_t sym;
	short col;
};

struct mesh {
	vector<triangle> tris;

	bool LoadFromObjectFile(string sFilename, bool bHasTexture = false) { // change to true if object file has textures
		ifstream f(sFilename);
		if (!f.is_open()) {
			// if we cannot open the file 
			return false;
		}

		// local cache of vertices 
		vector<vec3d> vertices;
		vector<vec2d> texs;

		while (!f.eof()) {
			char line[128]; // assume that each line is a maximum of 128 characters long
			f.getline(line, 128);

			stringstream s;
			s << line;

			char junk;

			if (line[0] == 'v') {
				if (line[1] == 't') {
					// reading texture
					vec2d v;
					s >> junk >> junk >> v.u >> v.v;
					texs.push_back(v);
				}
				else {
					// reading vertex
					vec3d v;
					s >> junk >> v.x >> v.y >> v.z;
					vertices.push_back(v);
				}
			}

			if (!bHasTexture) {
				if (line[0] == 'f') {
					// reading triangle face
					int f[3];
					s >> junk >> f[0] >> f[1] >> f[2];
					tris.push_back({ vertices[f[0] - 1], vertices[f[1] - 1], vertices[f[2] - 1] });
				}
			}
			else {
				if (line[0] == 'f') {
					s >> junk;
					string tokens[6];
					int nTokenCount = -1;

					while (!s.eof()) {
						char c = s.get();
						if (c == ' ' || c == '/') {
							nTokenCount++;
						}
						else {
							tokens[nTokenCount].append(1, c);
						}
					}

					tokens[nTokenCount].pop_back();
					tris.push_back({ vertices[stoi(tokens[0]) - 1], vertices[stoi(tokens[2]) - 1], vertices[stoi(tokens[4]) - 1],
						texs[stoi(tokens[1]) - 1], texs[stoi(tokens[3]) - 1], texs[stoi(tokens[5]) - 1] });
				}
			}
		}

		return true;
	}
};

struct mat4x4 {
	float m[4][4] = { 0 };
};

class engine3D : public olcConsoleGameEngine {
public:
	engine3D() {
		m_sAppName = L"3D Demo";
	}

private:
	mesh meshCube;
	mat4x4 matProj;
	vec3d vCamera;
	vec3d vLookDir;
	float fYaw; // rotation in the xz plane required to implement first-person camera
	float fTheta;

	olcSprite* sprTex1;

	vec3d Matrix_MultiplyVector(mat4x4& m, vec3d& i)
	{
		vec3d v;
		v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
		v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
		v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
		v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
		return v;
	}

	mat4x4 Matrix_MakeIdentity()
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4x4 Matrix_MakeRotationX(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[1][2] = sinf(fAngleRad);
		matrix.m[2][1] = -sinf(fAngleRad);
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4x4 Matrix_MakeRotationY(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][2] = sinf(fAngleRad);
		matrix.m[2][0] = -sinf(fAngleRad);
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4x4 Matrix_MakeRotationZ(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][1] = sinf(fAngleRad);
		matrix.m[1][0] = -sinf(fAngleRad);
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4x4 Matrix_MakeTranslation(float x, float y, float z)
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		matrix.m[3][0] = x;
		matrix.m[3][1] = y;
		matrix.m[3][2] = z;
		return matrix;
	}

	mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
	{
		float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
		mat4x4 matrix;
		matrix.m[0][0] = fAspectRatio * fFovRad;
		matrix.m[1][1] = fFovRad;
		matrix.m[2][2] = fFar / (fFar - fNear);
		matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
		matrix.m[2][3] = 1.0f;
		matrix.m[3][3] = 0.0f;
		return matrix;
	}

	mat4x4 Matrix_MultiplyMatrix(mat4x4& m1, mat4x4& m2)
	{
		mat4x4 matrix;
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
		return matrix;
	}

	mat4x4 Matrix_PointAt(vec3d& pos, vec3d& target, vec3d& up)
	{
		// pos - where the vector should be 
		// target is akin to forward vector 
		
		// Calculate new forward direction
		vec3d newForward = Vector_Sub(target, pos);
		newForward = Vector_Normalize(newForward);

		// Calculate new Up direction
		vec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
		vec3d newUp = Vector_Sub(up, a);
		newUp = Vector_Normalize(newUp);

		// New Right direction is easy, its just cross product
		vec3d newRight = Vector_CrossProduct(newUp, newForward);

		// Construct Dimensioning and Translation Matrix	
		mat4x4 matrix;
		matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
		return matrix;

	}

	// Only for Rotation/Translation Matrices
	mat4x4 Matrix_QuickInverse(mat4x4& m) {
		mat4x4 matrix;
		matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
		matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
		matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	vec3d Vector_Add(vec3d& v1, vec3d& v2) {
		return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
	}

	vec3d Vector_Sub(vec3d& v1, vec3d& v2)
	{
		return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
	}

	vec3d Vector_Mul(vec3d& v1, float k)
	{
		return { v1.x * k, v1.y * k, v1.z * k };
	}

	vec3d Vector_Div(vec3d& v1, float k)
	{
		return { v1.x / k, v1.y / k, v1.z / k };
	}

	float Vector_DotProduct(vec3d& v1, vec3d& v2)
	{
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	float Vector_Length(vec3d& v)
	{
		return sqrtf(Vector_DotProduct(v, v));
	}

	vec3d Vector_Normalize(vec3d& v)
	{
		float l = Vector_Length(v);
		return { v.x / l, v.y / l, v.z / l };
	}

	vec3d Vector_CrossProduct(vec3d& v1, vec3d& v2)
	{
		vec3d v;
		v.x = v1.y * v2.z - v1.z * v2.y;
		v.y = v1.z * v2.x - v1.x * v2.z;
		v.z = v1.x * v2.y - v1.y * v2.x;
		return v;
	}

	vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd, float &t) {
		// function will return a vector if the line crosses the plane
		plane_n = Vector_Normalize(plane_n);
		float plane_d = -Vector_DotProduct(plane_n, plane_p);
		float ad = Vector_DotProduct(lineStart, plane_n);
		float bd = Vector_DotProduct(lineEnd, plane_n);
		t = (-plane_d - ad) / (bd - ad); // normalized distance along line between two points where intersection has happened
		vec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
		vec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
		return Vector_Add(lineStart, lineToIntersect);
	}

	int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2) {
		// returns an integer representing number of triangles returned by the function

		// make sure plane normal is indeed normal 
		plane_n = Vector_Normalize(plane_n);

		// return signed shortest distance from point to plane, plane normal must be normalized 
		auto dist = [&](vec3d& p) {
			vec3d n = Vector_Normalize(p);
			return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
		};

		// create two temp storage arrays to classify points on either side of plane 
		// if distance sign is positive, point lies on inside of plane 
		vec3d* inside_points[3];
		int nInsidePointCount = 0;
		vec3d* outside_points[3];
		int nOutsidePointCount = 0;
		vec2d* inside_tex[3];
		int nInsideTexCount = 0;
		vec2d* outside_tex[3];
		int nOutsideTexCount = 0;

		// get signed distance of each point in triangle to plane 
		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		// determine if inside or outside point based on sign of distance 
		if (d0 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[0];
			inside_tex[nInsideTexCount++] = &in_tri.t[0];
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[0];
			outside_tex[nOutsideTexCount++] = &in_tri.t[0];
		}
		if (d1 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[1];
			inside_tex[nInsideTexCount++] = &in_tri.t[1];
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[1];
			outside_tex[nOutsideTexCount++] = &in_tri.t[0];
		}
		if (d2 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[2];
			inside_tex[nInsideTexCount++] = &in_tri.t[2];
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[2];
			outside_tex[nOutsideTexCount++] = &in_tri.t[0];
		}

		// classify triangle points and break triangle into smaller triangles if needed
		// four possible outcomes, inside points ranging from 0 to 3
		if (nInsidePointCount == 0) {
			// clip whole triangle since all points lie outside of plane 
			return 0;
		}
		if (nInsidePointCount == 1 && nOutsidePointCount == 2) {
			// triangle should be clipped, triangle turns into smaller triangle 

			// copy appearance info to new triangle 
			out_tri1.col = in_tri.col;
			out_tri1.sym = in_tri.sym;

			// keep valid inside point 
			out_tri1.p[0] = *inside_points[0];
			out_tri1.t[0] = *inside_tex[0];

			// two new points are at the locations where original sides of the triangle intersect with the plane 
			float t;
			out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[1].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[1].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[1].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);
			out_tri1.t[2].u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[2].w = t * (outside_tex[1]->w - inside_tex[0]->w) + inside_tex[0]->w;

			return 1; // the newly formed smaller triangle 
		}
		if (nInsidePointCount == 2 && nOutsidePointCount == 1) {
			// triangle should be clipped 
			// since two points lie inside the plane we need to create a quad with two new triangles 

			// copy appearance info to new triangle 
			out_tri1.col = in_tri.col;
			out_tri1.sym = in_tri.sym;
			out_tri2.col = in_tri.col;
			out_tri2.sym = in_tri.sym;

			// the first triangle consists of the two inside points and a new point determined by location where 
			//  one side of the triangle intersects with the plane 
			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.t[0] = *inside_tex[0];
			out_tri2.t[1] = *inside_tex[1];
			
			float t;
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[2].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[2].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			// second triangle consists of one inside point, a new point determined by the intersection of the other side 
			//   of the triangle and the plane, and the newly created point above
			out_tri2.p[0] = *inside_points[1];
			out_tri2.t[0] = *inside_tex[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.t[1] = out_tri1.t[2];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0], t);
			out_tri2.t[2].u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
			out_tri2.t[2].v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
			out_tri2.t[2].w = t * (outside_tex[0]->w - inside_tex[1]->w) + inside_tex[1]->w;

			return 2;
		}
		if (nInsidePointCount == 3) {
			// all points lie on inside of plane, allow triangle to pass through entirely 
			out_tri1 = in_tri;
			
			return 1; // only returning original triangle 
		}
	}

	// code taken from "Command Line Webcam?" on Youtube by javid9x
	// takes floating point of luminance value between 0 and 1 and returns the symbol and console colour combinations that represent it in white, gray, and black 
	CHAR_INFO GetColour(float lum)
	{
		short bg_col, fg_col;
		wchar_t sym;
		int pixel_bw = (int)(13.0f * lum);
		switch (pixel_bw)
		{
		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
		default:
			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
		}

		CHAR_INFO c;
		c.Attributes = bg_col | fg_col;
		c.Char.UnicodeChar = sym;
		return c;
	}

	float* pDepthBuffer = nullptr;

public:
	// two functions we must override

	bool OnUserCreate() override {
		// in this case we have manually created a cube, to display other shapes we can import object files
		//meshCube.LoadFromObjectFile("mountains.obj");
		
		
		meshCube.tris = {
			// creating a cube composed of triangles, every 3 floating point values represents (x, y, z)
			// SOUTH FACE
			{0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f},
			{0.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  1.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f},
			// EAST FACE
			{1.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f},
			{1.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f},
			// NORTH FACE
			{1.0f, 0.0f, 1.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f},
			{1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f},
			// WEST FACE
			{0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f},
			{0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f},
			// TOP FACE
			{0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f},
			{0.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f},
			// BOTTOM FACE
			{1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f},
			{1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f}
		};

		sprTex1 = new olcSprite(L"../ConsoleGame/Jario.spr");

		// Projection matrix 
		matProj = Matrix_MakeProjection(90.0f, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);

		return true;
	}
	bool OnUserUpdate(float fElapsedTime) override {
		// commands to move camera 
		if (GetKey(VK_UP).bHeld) {
			vCamera.y += 8.0f * fElapsedTime;
		} 
		if (GetKey(VK_DOWN).bHeld) {
			vCamera.y -= 8.0f * fElapsedTime;
		}
		// do not use in FPS mode 
		if (GetKey(VK_LEFT).bHeld) {
			vCamera.x -= 8.0f * fElapsedTime;
		}
		if (GetKey(VK_RIGHT).bHeld) {
			vCamera.x += 8.0f * fElapsedTime;
		}

		vec3d vForward = Vector_Mul(vLookDir, 8.0f * fElapsedTime);

		if (GetKey(L'W').bHeld) {
			vCamera = Vector_Add(vCamera, vForward);
		}
		if (GetKey(L'S').bHeld) {
			vCamera = Vector_Sub(vCamera, vForward);
		}

		if (GetKey(L'A').bHeld) {
			fYaw -= 2.0f * fElapsedTime;
		}
		if (GetKey(L'D').bHeld) {
			fYaw += 2.0f * fElapsedTime;
		}

		// perform rotations about the x and z axis, perform at different speed else we have gimbal lock
		// set up rotation matrices 
		mat4x4 matRotZ, matRotX;
		// fTheta += 1.0f * fElapsedTime; // uncomment to implement world spin

		// z-axis rotation 
		matRotZ = Matrix_MakeRotationZ(fTheta * 0.5f);

		// x-axis rotation
		matRotX = Matrix_MakeRotationX(fTheta);

		mat4x4 matTrans;
		matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 5.0f);

		mat4x4 matWorld;
		matWorld = Matrix_MakeIdentity();
		matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
		matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);

		vec3d vUp = { 0, 1, 0 };
		vec3d vTarget = {0, 0, 1};
		mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw); // rotate by "yaw" radians 
		vLookDir = Matrix_MultiplyVector(matCameraRot, vTarget); // unit vector rotated around the origin
		vTarget = Vector_Add(vCamera, vLookDir); // offset to current camera's location

		mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);

		// make view matrix from camera 
		mat4x4 matView = Matrix_QuickInverse(matCamera);

		// store triangles so they can be rasterized later
		vector<triangle> vecTrianglesToRaster;

		// draw triangles 
		for (auto tri : meshCube.tris) {
			// the objects exist in 3d space but we are in a 2d screen, need to use projection
			triangle triProjected, triTransformed, triViewed;

			// multiply each point of the triangle by the world matrix 
			triTransformed.p[0] = Matrix_MultiplyVector(matWorld, tri.p[0]);
			triTransformed.p[1] = Matrix_MultiplyVector(matWorld, tri.p[1]);
			triTransformed.p[2] = Matrix_MultiplyVector(matWorld, tri.p[2]);
			triTransformed.t[0] = tri.t[0];
			triTransformed.t[1] = tri.t[1];
			triTransformed.t[2] = tri.t[2];

			// calculate triangle's 
			vec3d normal, line1, line2;
			line1 = Vector_Sub(triTransformed.p[1], triTransformed.p[0]);
			line2 = Vector_Sub(triTransformed.p[2], triTransformed.p[0]);

			// take cross product of lines to get normal to triangle surface 
			normal = Vector_CrossProduct(line1, line2);

			// need to normalize normal 
			normal = Vector_Normalize(normal);

			// get ray from triangle to camera 
			vec3d vCameraRay = Vector_Sub(triTransformed.p[0], vCamera);


			// we only want to project the triangles if we can "see" them, that is the ray is aligned with normal 
			// calculate dot product between camera's field of view line and the normal and determine whether it is less than zero
			if (Vector_DotProduct(normal, vCameraRay) < 0.0f) {
				// illuminate the triangle
				vec3d light_direction = { 0.0f, 1.0f, -1.0f }; // all the light is coming from a single direction 
				light_direction = Vector_Normalize(light_direction);

				// determine alignment of light direction and triangle surface normal 
				float dotProduct = max(0.1f, Vector_DotProduct(light_direction, normal));
				
				// choose console colours as required 
				CHAR_INFO c = GetColour(dotProduct);
				triTransformed.col = c.Attributes;
				triTransformed.sym = c.Char.UnicodeChar;

				// convert world space to view space 
				triViewed.p[0] = Matrix_MultiplyVector(matView, triTransformed.p[0]);
				triViewed.p[1] = Matrix_MultiplyVector(matView, triTransformed.p[1]);
				triViewed.p[2] = Matrix_MultiplyVector(matView, triTransformed.p[2]);
				triViewed.col = triTransformed.col;
				triViewed.sym = triTransformed.sym;
				triViewed.t[0] = triTransformed.t[0];
				triViewed.t[1] = triTransformed.t[1];
				triViewed.t[2] = triTransformed.t[2];

				// clip viewed triangle against near plane to form up to two additional triangles 
				int nClippedTriangles = 0;
				triangle clipped[2];
				// first two parameters refer to a point on the near-plane (just in-front of the camera on z-axis) 
				//   and the normal straight along the z-axis
				nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

				for (int n = 0; n < nClippedTriangles; n++) {
					// project triangles from 3d into 2d
					triProjected.p[0] = Matrix_MultiplyVector(matProj, clipped[n].p[0]);
					triProjected.p[1] = Matrix_MultiplyVector(matProj, clipped[n].p[1]);
					triProjected.p[2] = Matrix_MultiplyVector(matProj, clipped[n].p[2]);
					triProjected.col = clipped[n].col;
					triProjected.sym = clipped[n].sym;
					triProjected.t[0] = clipped[n].t[0];
					triProjected.t[1] = clipped[n].t[1];
					triProjected.t[2] = clipped[n].t[2];

					triProjected.t[0].u = triProjected.t[0].u / triProjected.p[0].w;
					triProjected.t[1].u = triProjected.t[1].u / triProjected.p[1].w;
					triProjected.t[2].u = triProjected.t[2].u / triProjected.p[2].w;

					triProjected.t[0].v = triProjected.t[0].v / triProjected.p[0].w;
					triProjected.t[1].v = triProjected.t[1].v / triProjected.p[1].w;
					triProjected.t[2].v = triProjected.t[2].v / triProjected.p[2].w;

					triProjected.t[0].w = 1.0f / triProjected.p[0].w;
					triProjected.t[1].w = 1.0f / triProjected.p[1].w;
					triProjected.t[2].w = 1.0f / triProjected.p[2].w;

					// scale into view 
					triProjected.p[0] = Vector_Div(triProjected.p[0], triProjected.p[0].w);
					triProjected.p[1] = Vector_Div(triProjected.p[1], triProjected.p[1].w);
					triProjected.p[2] = Vector_Div(triProjected.p[2], triProjected.p[2].w);

					// X/Y are inverted so put them back
					triProjected.p[0].x *= -1.0f;
					triProjected.p[1].x *= -1.0f;
					triProjected.p[2].x *= -1.0f;
					triProjected.p[0].y *= -1.0f;
					triProjected.p[1].y *= -1.0f;
					triProjected.p[2].y *= -1.0f;

					// offset vertices into visible normalized space 
					vec3d vOffsetView = { 1, 1, 0 };
					triProjected.p[0] = Vector_Add(triProjected.p[0], vOffsetView);
					triProjected.p[1] = Vector_Add(triProjected.p[1], vOffsetView);
					triProjected.p[2] = Vector_Add(triProjected.p[2], vOffsetView);

					triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

					// store triangles for sorting 
					vecTrianglesToRaster.push_back(triProjected);


					//// rasterize triangle
					//FillTriangle(triProjected.p[0].x, triProjected.p[0].y,
					//	triProjected.p[1].x, triProjected.p[1].y,
					//	triProjected.p[2].x, triProjected.p[2].y,
					//	triProjected.sym, triProjected.col);

					// useful for debugging, determine wireframe outline 
					/*DrawTriangle(triProjected.p[0].x, triProjected.p[0].y,
						triProjected.p[1].x, triProjected.p[1].y,
						triProjected.p[2].x, triProjected.p[2].y,
						PIXEL_SOLID, FG_WHITE); */
				}
			}
		}

		// sort triangles from back to front since we are drawing triangles closer to camera on top of triangles further away
		//sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(), [](triangle& t1, triangle& t2) {
		//	// get the midpoint z-values of both triangles
		//	float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
		//	float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
		//	return z1 > z2;
		//});

		// clear screen 
		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		// clear depth buffer
		for (int i = 0; i < ScreenWidth() * ScreenHeight(); i++) {
			pDepthBuffer[i] = 0.0f;
		}	

		// rasterize triangles 
		for (auto& triToRaster : vecTrianglesToRaster) {
			// clip triangles against all four screen edges to yield many triangles
			triangle clipped[2];
			list<triangle> listTriangles; // create a list to ensure we only test new triangles generated against planes

			// Add initial triangle
			listTriangles.push_back(triToRaster);
			int nNewTriangles = 1;

			for (int p = 0; p < 4; p++)
			{
				int nTrisToAdd = 0;
				while (nNewTriangles > 0)
				{
					// Take triangle from front of queue
					triangle test = listTriangles.front();
					listTriangles.pop_front();
					nNewTriangles--;

					// Clip it against a plane. We only need to test each 
					// subsequent plane, against subsequent new triangles
					// as all triangles after a plane clip are guaranteed
					// to lie on the inside of the plane. 
					switch (p)
					{
						// all four screen edges
						case 0:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
						case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
						case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
						case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					}

					// Clipping may yield a variable number of triangles, so
					// add these new ones to the back of the queue for subsequent
					// clipping against next planes
					for (int w = 0; w < nTrisToAdd; w++)
						listTriangles.push_back(clipped[w]);
				}
				nNewTriangles = listTriangles.size();
			}

			// Draw the transformed, viewed, clipped, projected, and sorted triangles
			for (auto& t : listTriangles)
			{
				TexturedTriangle(t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
					t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
					t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, sprTex1);
				// FillTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, t.sym, t.col);
				DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, PIXEL_SOLID, FG_WHITE); // if you would like to see the wireframe
			}
		}

		return true;
	}

	void TexturedTriangle(int x1, int y1, float u1, float v1, float w1,
		int x2, int y2, float u2, float v2, float w2, 
		int x3, int y3, float u3, float v3, float w3,
		olcSprite* tex) {
		// sort x-y coordinates from high to low y values
		if (y2 < y1) {
			swap(y1, y2);
			swap(x1, x2);
			swap(u1, u2);
			swap(v1, v2);
			swap(w1, w2);
		}
		if (y3 < y1) {
			swap(y1, y3);
			swap(x1, x3);
			swap(u1, u3);
			swap(v1, v3);
			swap(w1, w3);
		}
		if (y3 < y2) {
			swap(y2, y3);
			swap(x2, x3);
			swap(u2, u3);
			swap(v2, v3);
			swap(w2, w3);
		}

		// gradient and utility information
		int dy1 = y2 - y1;
		int dx1 = x2 - x1; 
		float dv1 = v2 - v1;
		float du1 = u2 - u1;
		float dw1 = w2 - w1;

		int dy2 = y3 - y1;
		int dx2 = x3 - x1;
		float dv2 = v3 - v1;
		float du2 = u3 - u1;
		float dw2 = w3 - w1;

		float tex_u, tex_v, tex_w;

		float dax_step = 0, dbx_step = 0, du1_step = 0, dv1_step = 0, du2_step = 0, dv2_step = 0, dw1_step = 0, dw2_step = 0;

		if (dy1 != 0) {
			dax_step = dx1 / (float)abs(dy1);
			du1_step = du1 / (float)abs(dy1);
			dv1_step = dv1 / (float)abs(dy1);
			dw1_step = dw1 / (float)abs(dy1);
		}
		if (dy2 != 0) {
			dbx_step = dx2 / (float)abs(dy2);
			du2_step = du2 / (float)abs(dy2);
			dv2_step = dv2 / (float)abs(dy2);
			dw2_step = dw2 / (float)abs(dy2);
		}
		
		if (dy1 != 0) {
			for (int i = y2; i <= y3; i++) {
				int ax = x2 + (float)(i - y2) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				// sending texture coordinates
				// start u, v
				float tex_su = u2 + (float)(i - y2) * du1_step;
				float tex_sv = v2 + (float)(i - y2) * dv1_step;
				float tex_sw = w2 + (float)(i - y2) * dw1_step;
				 
				// end u, v
				float tex_eu = u1 + (float)(i - y1) * du2_step;
				float tex_ev = v1 + (float)(i - y1) * dv2_step;
				float tex_ew = w1 + (float)(i - y1) * dw2_step;

				if (ax > bx) {
					swap(ax, bx);
					swap(tex_su, tex_eu);
					swap(tex_sv, tex_ev);
					swap(tex_sw, tex_ew);
				}

				tex_u = tex_su;
				tex_v = tex_sv;
				tex_w = tex_sw;

				float tstep = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				// sampling the texture: each one of these locations represents a pixel on the screen
				for (int j = ax; j < bx; j++)
				{
					tex_u = (1.0f - t) * tex_su + t * tex_eu;
					tex_v = (1.0f - t) * tex_sv + t * tex_ev;
					tex_w = (1.0f - t) * tex_sw + t * tex_ew;

					if (tex_w > pDepthBuffer[i * ScreenWidth() + j])
					{
						// draw a single pixel 
						Draw(j, i, tex->SampleGlyph(tex_u / tex_w, tex_v / tex_w), tex->SampleColour(tex_u / tex_w, tex_v / tex_w));
						pDepthBuffer[i * ScreenWidth() + j] = tex_w;
					}
					t += tstep;
				}
			}
		}
	}
};



int main()
{
	engine3D demo;
	if (demo.ConstructConsole(300, 200, 2, 2)) { // width, height, fontw, fonth
		demo.Start();
	}


	return 0;
}