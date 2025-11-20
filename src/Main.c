#if defined __linux__
#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#elif defined _WINE
#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#elif defined _WIN32
#include "F:/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#elif defined(__APPLE__)
#error "Apple not supported!"
#else
#error "Platform not supported!"
#endif

typedef struct Vec3D{
	float x;
	float y;
	float z;
	float w;
} Vec3D;

Vec3D Vec3D_new(float x,float y,float z){
	return (Vec3D){ x,y,z,1.0f };
}

typedef struct Tri3D
{
	Vec3D p[3];
	Pixel c;
} Tri3D;

typedef struct mesh{
	Vector tris;
} mesh;

typedef struct M4x4D{
	float m[4][4];
} M4x4D;

mesh meshCube;
M4x4D matProj;
Vec3D vCamera;
Vec3D vLookDir;
float fYaw;	
float fTheta;

M4x4D Matrix_Null()
{
	return (M4x4D){{
		{ 0.0f,0.0f,0.0f,0.0f },
		{ 0.0f,0.0f,0.0f,0.0f },
		{ 0.0f,0.0f,0.0f,0.0f },
		{ 0.0f,0.0f,0.0f,0.0f }
	}};
}
Vec3D Matrix_MultiplyVector(M4x4D m, Vec3D i)
{
	Vec3D v;
	v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
	v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
	v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
	v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
	return v;
}
M4x4D Matrix_MakeIdentity()
{
	M4x4D matrix = Matrix_Null();
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}
M4x4D Matrix_MakeRotationX(float fAngleRad)
{
	M4x4D matrix = Matrix_Null();
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[1][2] = sinf(fAngleRad);
	matrix.m[2][1] = -sinf(fAngleRad);
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0f;
	return matrix;
}
M4x4D Matrix_MakeRotationY(float fAngleRad)
{
	M4x4D matrix = Matrix_Null();
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][2] = sinf(fAngleRad);
	matrix.m[2][0] = -sinf(fAngleRad);
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0f;
	return matrix;
}
M4x4D Matrix_MakeRotationZ(float fAngleRad)
{
	M4x4D matrix = Matrix_Null();
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][1] = sinf(fAngleRad);
	matrix.m[1][0] = -sinf(fAngleRad);
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}
M4x4D Matrix_MakeTranslation(float x, float y, float z)
{
	M4x4D matrix = Matrix_Null();
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	matrix.m[3][0] = x;
	matrix.m[3][1] = y;
	matrix.m[3][2] = z;
	return matrix;
}
M4x4D Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
{
	float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
	M4x4D matrix = Matrix_Null();
	matrix.m[0][0] = fAspectRatio * fFovRad;
	matrix.m[1][1] = fFovRad;
	matrix.m[2][2] = fFar / (fFar - fNear);
	matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matrix.m[2][3] = 1.0f;
	matrix.m[3][3] = 0.0f;
	return matrix;
}
M4x4D Matrix_MultiplyMatrix(M4x4D m1, M4x4D m2)
{
	M4x4D matrix = Matrix_Null();
	for (int c = 0; c < 4; c++)
		for (int r = 0; r < 4; r++)
			matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
	return matrix;
}

Vec3D Vec3D_Add(Vec3D v1, Vec3D v2)
{
	return Vec3D_new( v1.x + v2.x, v1.y + v2.y, v1.z + v2.z );
}
Vec3D Vec3D_Sub(Vec3D v1, Vec3D v2)
{
	return Vec3D_new( v1.x - v2.x, v1.y - v2.y, v1.z - v2.z );
}
Vec3D Vec3D_Mul(Vec3D v1, float k)
{
	return Vec3D_new( v1.x * k, v1.y * k, v1.z * k );
}
Vec3D Vec3D_Div(Vec3D v1, float k)
{
	if(k==0.0f) return Vec3D_new( 0.0f,0.0f,0.0f );
	return Vec3D_new( v1.x / k, v1.y / k, v1.z / k );
}
float Vec3D_DotProduct(Vec3D v1, Vec3D v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z * v2.z;
}
float Vec3D_Length(Vec3D v)
{
	return sqrtf(Vec3D_DotProduct(v, v));
}
Vec3D Vec3D_Normalise(Vec3D v)
{
	float l = Vec3D_Length(v);
	return Vec3D_new( v.x / l, v.y / l, v.z / l );
}
Vec3D Vec3D_Perp(Vec3D v)
{
	return Vec3D_new( -v.z, v.y, v.x );
}

Vec3D Vec3D_CrossProduct(Vec3D v1, Vec3D v2)
{
	Vec3D v = Vec3D_new(0.0f,0.0f,0.0f);
	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;
	return v;
}
Vec3D Vec3D_IntersectPlane(Vec3D plane_p, Vec3D* plane_n, Vec3D lineStart, Vec3D lineEnd)
{
	*plane_n = Vec3D_Normalise(*plane_n);
	float plane_d = -Vec3D_DotProduct(*plane_n, plane_p);
	float ad = Vec3D_DotProduct(lineStart, *plane_n);
	float bd = Vec3D_DotProduct(lineEnd, *plane_n);
	float t = (-plane_d - ad) / (bd - ad);
	Vec3D lineStartToEnd = Vec3D_Sub(lineEnd, lineStart);
	Vec3D lineToIntersect = Vec3D_Mul(lineStartToEnd, t);
	return Vec3D_Add(lineStart, lineToIntersect);
}

float dist(Vec3D plane_p, Vec3D plane_n,Vec3D p){
	Vec3D n = Vec3D_Normalise(p);
	return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vec3D_DotProduct(plane_n, plane_p));
}

int Triangle_ClipAgainstPlane(Vec3D plane_p, Vec3D plane_n, Tri3D in_tri, Tri3D* out_tri1, Tri3D* out_tri2){
	// Make sure plane normal is indeed normal
	plane_n = Vec3D_Normalise(plane_n);
	
	// Create two temporary storage arrays to classify points either side of plane
	// If distance sign is positive, point lies on "inside" of plane
	Vec3D* inside_points[3];  int nInsidePointCount = 0;
	Vec3D* outside_points[3]; int nOutsidePointCount = 0;
	// Get signed distance of each point in Tri3D to plane
	float d0 = dist(plane_p,plane_n,in_tri.p[0]);
	float d1 = dist(plane_p,plane_n,in_tri.p[1]);
	float d2 = dist(plane_p,plane_n,in_tri.p[2]);
	if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
	if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
	if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }
	// Now classify Tri3D points, and break the input Tri3D into 
	// smaller output Tri3Ds if required. There are four possible
	// outcomes...
	if (nInsidePointCount == 0)
	{
		// All points lie on the outside of plane, so clip whole Tri3D
		// It ceases to exist
		return 0; // No returned Tri3Ds are valid
	}
	if (nInsidePointCount == 3)
	{
		// All points lie on the inside of plane, so do nothing
		// and allow the Tri3D to simply pass through
		*out_tri1 = in_tri;
		return 1; // Just the one returned original Tri3D is valid
	}
	if (nInsidePointCount == 1 && nOutsidePointCount == 2)
	{
		// Triangle should be clipped. As two points lie outside
		// the plane, the Tri3D simply becomes a smaller Tri3D
		// Copy appearance info to new Tri3D
		out_tri1->c =  in_tri.c;
		// The inside point is valid, so keep that...
		out_tri1->p[0] = *inside_points[0];
		// but the two new points are at the locations where the 
		// original sides of the Tri3D (lines) intersect with the plane
		out_tri1->p[1] = Vec3D_IntersectPlane(plane_p, &plane_n, *inside_points[0], *outside_points[0]);
		out_tri1->p[2] = Vec3D_IntersectPlane(plane_p, &plane_n, *inside_points[0], *outside_points[1]);
		return 1; // Return the newly formed single Tri3D
	}
	if (nInsidePointCount == 2 && nOutsidePointCount == 1)
	{
		// Triangle should be clipped. As two points lie inside the plane,
		// the clipped Tri3D becomes a "quad". Fortunately, we can
		// represent a quad with two new Tri3Ds
		// Copy appearance info to new Tri3Ds
		out_tri1->c =  in_tri.c;
		out_tri2->c =  in_tri.c;
		// The first Tri3D consists of the two inside points and a new
		// point determined by the location where one side of the Tri3D
		// intersects with the plane
		out_tri1->p[0] = *inside_points[0];
		out_tri1->p[1] = *inside_points[1];
		out_tri1->p[2] = Vec3D_IntersectPlane(plane_p, &plane_n, *inside_points[0], *outside_points[0]);
		// The second Tri3D is composed of one of he inside points, a
		// new point determined by the intersection of the other side of the 
		// Tri3D and the plane, and the newly created point above
		out_tri2->p[0] = *inside_points[1];
		out_tri2->p[1] = out_tri1->p[2];
		out_tri2->p[2] = Vec3D_IntersectPlane(plane_p, &plane_n, *inside_points[1], *outside_points[0]);
		return 2; // Return two newly formed Tri3Ds which form a quad
	}
}

M4x4D Matrix_PointAt(Vec3D pos, Vec3D target, Vec3D up)
{
	// Calculate new forward direction
	Vec3D newForward = Vec3D_Sub(target, pos);
	newForward = Vec3D_Normalise(newForward);
	// Calculate new Up direction
	Vec3D a = Vec3D_Mul(newForward, Vec3D_DotProduct(up,newForward));
	Vec3D newUp = Vec3D_Sub(up, a);
	newUp = Vec3D_Normalise(newUp);
	// New Right direction is easy, its just cross product
	Vec3D newRight = Vec3D_CrossProduct(newUp, newForward);
	// Construct Dimensioning and Translation Matrix	
	M4x4D matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;
}
M4x4D Matrix_QuickInverse(M4x4D m) // Only for Rotation/Translation Matrices
{
	M4x4D matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
	matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
	matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
	matrix.m[3][3] = 1.0f;
	return matrix;
}

typedef struct Rect3 {
	Vec3D p;
	Vec3D d;
} Rect3;

Bool Rect3_Overlap(Rect3 r1,Rect3 r2){
	return !(r1.p.x<r2.p.x-r1.d.x || r1.p.y<r2.p.y-r1.d.y || r1.p.z<r2.p.z-r1.d.z || r1.p.x>r2.p.x+r2.d.x || r1.p.y>r2.p.y+r2.d.y || r1.p.z>r2.p.z+r2.d.z);
}
void Rect3_Static(Rect3* r1,Rect3 r2){
	if(Rect3_Overlap(*r1,r2)){
        Vec3D m1 = Vec3D_Add(r1->p,Vec3D_Mul(r1->d,0.5f));
        Vec3D m2 = Vec3D_Add(r2.p, Vec3D_Mul(r2.d, 0.5f));

        Vec3D l = Vec3D_Add(r1->d,r2.d);
        Vec3D d = Vec3D_Sub(m2,m1);
        d = Vec3D_new(d.x / l.x,d.y / l.y,d.z / l.z);

        if(F32_Abs(d.x)>F32_Abs(d.y)){
            if(F32_Abs(d.x)>F32_Abs(d.z)){
				if(d.x>0.0f){
                	m1.x = m2.x - l.x * 0.5f;
                	r1->p.x = m1.x - r1->d.x * 0.5f;
            	}else{
            	    m1.x = m2.x + l.x * 0.5f;
            	    r1->p.x = m1.x - r1->d.x * 0.5f;
            	}
			}else{
				if(d.z>0.0f){
                	m1.z = m2.z - l.z * 0.5f;
                	r1->p.z = m1.z - r1->d.z * 0.5f;
            	}else{
            	    m1.z = m2.z + l.z * 0.5f;
            	    r1->p.z = m1.z - r1->d.z * 0.5f;
            	}
			}
        }else{
            if(F32_Abs(d.y)>F32_Abs(d.z)){
				if(d.y>0.0f){
                	m1.y = m2.y - l.y * 0.5f;
                	r1->p.y = m1.y - r1->d.y * 0.5f;
            	}else{
            	    m1.y = m2.y + l.y * 0.5f;
            	    r1->p.y = m1.y - r1->d.y * 0.5f;
            	}
			}else{
				if(d.z>0.0f){
                	m1.z = m2.z - l.z * 0.5f;
                	r1->p.z = m1.z - r1->d.z * 0.5f;
            	}else{
            	    m1.z = m2.z + l.z * 0.5f;
            	    r1->p.z = m1.z - r1->d.z * 0.5f;
            	}
			}
        }
    }
}

#define PLANE_FRONT 	0
#define PLANE_LEFT 		1
#define PLANE_BACK 		2
#define PLANE_RIGHT 	3
#define PLANE_TOP 		4
#define PLANE_BOTTOM 	5

mesh meshCube;
M4x4D matProj;
Vec3D vCamera = { 0.0f,0.0f,-1.0f,1.0f };
Vec3D vLookDir;
float fYaw;	
float fPitch;	
float fTheta;

Vec3D vLength = { 0.5f,1.8f,0.5f,1.0f };
Vector Cubes;
Vec2 MouseBefore = { 0.0f,0.0f };

void MakeCube(Vec3D p,Vec3D d,Pixel c){
	Tri3D tris[12] = {
	// SOUTH
	{ 0.0f, 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,	c },
	{ 0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,	c },
	// EAST                                                      
	{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,	c },
	{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f, 1.0f,	c },
	// NORTH                                                     
	{ 1.0f, 0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,	c },
	{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,	c },
	// WEST                                                      
	{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,	c },
	{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,	c },
	// TOP                                                       
	{ 0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,	c },
	{ 0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,	c },
	// BOTTOM                                                    
	{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,	c },
	{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,	c },
	};

	for(int i = 0;i<12;i++){
		for(int j = 0;j<3;j++){
			tris[i].p[j] = Vec3D_Add(p,Vec3D_new(tris[i].p[j].x * d.x,tris[i].p[j].y * d.y,tris[i].p[j].z * d.z));
		}
	}
	Vector_PushCount(&meshCube.tris,tris,12);
}

void MakePlane(Vec3D p,Vec3D d,int Plane,Pixel c){
	Tri3D tris[12] = {
	// SOUTH
	{ 0.0f, 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,	c },
	{ 0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,	c },
	// EAST                                                      
	{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,	c },
	{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f, 1.0f,	c },
	// NORTH                                                     
	{ 1.0f, 0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,	c },
	{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,	c },
	// WEST                                                      
	{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,	c },
	{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,	c },
	// TOP                                                       
	{ 0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,	c },
	{ 0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,	c },
	// BOTTOM                                                    
	{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,	c },
	{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,	c },
	};

	for(int i = Plane*2;i<(Plane+1)*2;i++){
		for(int j = 0;j<3;j++){
			tris[i].p[j] = Vec3D_Add(p,Vec3D_new(tris[i].p[j].x * d.x,tris[i].p[j].y * d.y,tris[i].p[j].z * d.z));
		}
		Vector_Push(&meshCube.tris,&tris[i]);
	}
}

void BuildCube(Vec3D p,Vec3D d,Pixel c){
	Vector_Push(&Cubes,(Rect3[]){ { p,d } });
	MakeCube(p,d,c);
}

void Setup(AlxWindow* w){
	meshCube = (mesh){ Vector_New(sizeof(Tri3D)) };
	Cubes = Vector_New(sizeof(Rect3));

	//meshCube.LoadFromObjectFile("mountains.obj");

	//MakeCube((Vec3D){ 1.0f,0.0f,0.0f },(Vec3D){ 2.0f,2.0f,2.0f },WHITE);
	//for(int i = 0;i<100;i++){
	//	for(int j = 0;j<10;j++){
	//		MakePlane((Vec3D){ i,0.0f,j },(Vec3D){ 1.0f,0.0f,1.0f },PLANE_TOP,WHITE);
	//	}
	//}

	BuildCube((Vec3D){ 0.0f,0.0f,0.0f },(Vec3D){ 1.0f,1.0f,1.0f },WHITE);
	BuildCube((Vec3D){ 2.0f,0.0f,0.0f },(Vec3D){ 1.0f,1.0f,1.0f },WHITE);
	BuildCube((Vec3D){ 0.0f,0.0f,2.0f },(Vec3D){ 1.0f,1.0f,1.0f },WHITE);
	BuildCube((Vec3D){ 2.0f,0.0f,2.0f },(Vec3D){ 1.0f,1.0f,1.0f },WHITE);

	matProj = Matrix_MakeProjection(90.0f, (float)GetHeight() / (float)GetWidth(), 0.1f, 1000.0f);
}

int compare(const void* e1,const void* e2) {
	Tri3D t1 = *(Tri3D*)e1;
	Tri3D t2 = *(Tri3D*)e2;
	float z1 = (t1.p[0].z+t1.p[1].z+t1.p[2].z)/3;
    float z2 = (t2.p[0].z+t2.p[1].z+t2.p[2].z)/3;
    return z1 == z2 ? 0 : (z1 < z2 ? 1 : -1);
}

void Update(AlxWindow* w){
	for(int i = 0;i<Cubes.size;i++){
		//Rect3 r1 = *(Rect3*)Vector_Get(&Cubes,i);
		//Rect3 r2 = (Rect3){ Vec3D_Sub(vCamera,Vec3D_Mul(vLength,0.5f)),vLength };
		//Rect3_Static(&r2,r1);
		//vCamera = r2.p;
	}
	
	w->ElapsedTime = 0.005;

	if(Stroke(ALX_MOUSE_L).PRESSED){
		MouseBefore = GetMouse();
	}
	if(Stroke(ALX_MOUSE_L).DOWN){
		if(GetMouse().x!=MouseBefore.x || GetMouse().y!=MouseBefore.y){
			Vec2 d = Vec2_Sub(GetMouse(),MouseBefore);
			Vec2 a = Vec2_Mulf(Vec2_Div(d,(Vec2){ window.Width,window.Height }),2 * F32_PI);
	
			fYaw += a.x;
			fPitch += a.y;
	
			MouseBefore = GetMouse();
		}
	}
	
	if (Stroke(ALX_KEY_R).DOWN)
		vCamera.y += 8.0f * w->ElapsedTime;

	if (Stroke(ALX_KEY_F).DOWN)
		vCamera.y -= 8.0f * w->ElapsedTime;

	if (Stroke(ALX_KEY_UP).DOWN)
		fPitch -= 2.0f * w->ElapsedTime;

	if (Stroke(ALX_KEY_DOWN).DOWN)
		fPitch += 2.0f * w->ElapsedTime;

	if (Stroke(ALX_KEY_LEFT).DOWN)
		fYaw -= 2.0f * w->ElapsedTime;

	if (Stroke(ALX_KEY_RIGHT).DOWN)
		fYaw += 2.0f * w->ElapsedTime;
	

	M4x4D matCameraRot = Matrix_MakeRotationY(fYaw);
	Vec3D vForward = Vec3D_Mul(Matrix_MultiplyVector(matCameraRot, Vec3D_new( 0.0f,0.0f,1.0f )), 8.0f * w->ElapsedTime);
	
	if (Stroke(ALX_KEY_W).DOWN)
		vCamera = Vec3D_Add(vCamera, vForward);

	if (Stroke(ALX_KEY_S).DOWN)
		vCamera = Vec3D_Sub(vCamera, vForward);

	if (Stroke(ALX_KEY_A).DOWN)
		vCamera = Vec3D_Sub(vCamera, Vec3D_Perp(vForward));

	if (Stroke(ALX_KEY_D).DOWN)
		vCamera = Vec3D_Add(vCamera, Vec3D_Perp(vForward));

	M4x4D matRotZ, matRotX;
	
	fTheta += 1.0f * (float)w->ElapsedTime;
	matRotZ = Matrix_MakeRotationZ(fTheta * 0.5f);
	matRotX = Matrix_MakeRotationX(fTheta);

	M4x4D matTrans;
	matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 0.0f);

	M4x4D matWorld;
	matWorld = Matrix_MakeIdentity();
	//matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
	matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);

	Vec3D vUp = Vec3D_new( 0.0f,1.0f,0.0f );
	Vec3D vTarget = Vec3D_new( 0.0f,0.0f,1.0f );
	M4x4D matCameraRotX = Matrix_MakeRotationX(fPitch);
	vLookDir = Matrix_MultiplyVector(matCameraRotX, vTarget);
	vLookDir = Matrix_MultiplyVector(matCameraRot, vLookDir);
	
	vTarget = Vec3D_Add(vCamera, vLookDir);
	M4x4D matCamera = Matrix_PointAt(vCamera, vTarget, vUp);
	M4x4D matView = Matrix_QuickInverse(matCamera);

	Vector vecTrianglesToRaster = Vector_New(sizeof(Tri3D));

	for (int i = 0;i<meshCube.tris.size;i++){
		Tri3D tri = *(Tri3D*)Vector_Get(&meshCube.tris,i);
		Tri3D triProjected, triTransformed, triViewed;

		triTransformed.p[0] = Matrix_MultiplyVector(matWorld, tri.p[0]);
		triTransformed.p[1] = Matrix_MultiplyVector(matWorld, tri.p[1]);
		triTransformed.p[2] = Matrix_MultiplyVector(matWorld, tri.p[2]);

		Vec3D normal, line1, line2;

		line1 = Vec3D_Sub(triTransformed.p[1], triTransformed.p[0]);
		line2 = Vec3D_Sub(triTransformed.p[2], triTransformed.p[0]);

		normal = Vec3D_CrossProduct(line1, line2);

		normal = Vec3D_Normalise(normal);
		
		Vec3D vCameraRay = Vec3D_Sub(triTransformed.p[0], vCamera);

		if (Vec3D_DotProduct(normal, vCameraRay) < 0.0f){
			Vec3D light_direction = Vec3D_new( 0.0f, 1.0f, -1.0f );
			light_direction = Vec3D_Normalise(light_direction);

			float dp = fmaxf(0.1f, Vec3D_DotProduct(light_direction, normal));

			//triTransformed.c = Pixel_Mul(triTransformed.c,Pixel_toRGBA(dp,dp,dp,1.0f));
			triTransformed.c = Pixel_toRGBA(dp,dp,dp,1.0f);

			triViewed.p[0] = Matrix_MultiplyVector(matView, triTransformed.p[0]);
			triViewed.p[1] = Matrix_MultiplyVector(matView, triTransformed.p[1]);
			triViewed.p[2] = Matrix_MultiplyVector(matView, triTransformed.p[2]);
			triViewed.c = triTransformed.c;

			int nClippedTriangles = 0;
			Tri3D clipped[2];
			nClippedTriangles = Triangle_ClipAgainstPlane(Vec3D_new( 0.0f, 0.0f, 0.1f ), Vec3D_new( 0.0f, 0.0f, 1.0f ), triViewed, &clipped[0], &clipped[1]);

			for (int n = 0; n < nClippedTriangles; n++){
				triProjected.p[0] = Matrix_MultiplyVector(matProj, clipped[n].p[0]);
				triProjected.p[1] = Matrix_MultiplyVector(matProj, clipped[n].p[1]);
				triProjected.p[2] = Matrix_MultiplyVector(matProj, clipped[n].p[2]);
				triProjected.c = clipped[n].c;

				triProjected.p[0] = Vec3D_Div(triProjected.p[0], triProjected.p[0].w);
				triProjected.p[1] = Vec3D_Div(triProjected.p[1], triProjected.p[1].w);
				triProjected.p[2] = Vec3D_Div(triProjected.p[2], triProjected.p[2].w);

				triProjected.p[0].x *= -1.0f;
				triProjected.p[1].x *= -1.0f;
				triProjected.p[2].x *= -1.0f;
				triProjected.p[0].y *= -1.0f;
				triProjected.p[1].y *= -1.0f;
				triProjected.p[2].y *= -1.0f;

				Vec3D vOffsetView = Vec3D_new( 1,1,0 );
				triProjected.p[0] = Vec3D_Add(triProjected.p[0], vOffsetView);
				triProjected.p[1] = Vec3D_Add(triProjected.p[1], vOffsetView);
				triProjected.p[2] = Vec3D_Add(triProjected.p[2], vOffsetView);
				triProjected.p[0].x *= 0.5f * (float)GetWidth();
				triProjected.p[0].y *= 0.5f * (float)GetHeight();
				triProjected.p[1].x *= 0.5f * (float)GetWidth();
				triProjected.p[1].y *= 0.5f * (float)GetHeight();
				triProjected.p[2].x *= 0.5f * (float)GetWidth();
				triProjected.p[2].y *= 0.5f * (float)GetHeight();

				Vector_Push(&vecTrianglesToRaster,&triProjected);
			}			
		}
	}

	qsort(vecTrianglesToRaster.Memory,vecTrianglesToRaster.size,vecTrianglesToRaster.ELEMENT_SIZE,compare);

	Clear(BLACK);

	for (int i = 0;i<vecTrianglesToRaster.size;i++)
	{
		Tri3D triToRaster = *(Tri3D*)Vector_Get(&vecTrianglesToRaster,i);

		Tri3D clipped[2];
		Vector listTriangles = Vector_New(sizeof(Tri3D));

		Vector_Push(&listTriangles,&triToRaster);
		int nNewTriangles = 1;

		for (int p = 0; p < 4; p++)
		{
			int nTrisToAdd = 0;
			while (nNewTriangles > 0)
			{
				Tri3D test = *(Tri3D*)Vector_Get(&listTriangles,0);
				Vector_Remove(&listTriangles,0);
				nNewTriangles--;

				switch (p)
				{
				case 0:	nTrisToAdd = Triangle_ClipAgainstPlane(Vec3D_new( 0.0f, 0.0f, 0.0f ), 					Vec3D_new( 0.0f, 1.0f, 0.0f ), 	test, &clipped[0], &clipped[1]); break;
				case 1:	nTrisToAdd = Triangle_ClipAgainstPlane(Vec3D_new( 0.0f, (float)GetHeight() - 1, 0.0f ), Vec3D_new( 0.0f, -1.0f, 0.0f ), test, &clipped[0], &clipped[1]); break;
				case 2:	nTrisToAdd = Triangle_ClipAgainstPlane(Vec3D_new( 0.0f, 0.0f, 0.0f ), 					Vec3D_new( 1.0f, 0.0f, 0.0f ), 	test, &clipped[0], &clipped[1]); break;
				case 3:	nTrisToAdd = Triangle_ClipAgainstPlane(Vec3D_new( (float)GetWidth() - 1, 0.0f, 0.0f ), 	Vec3D_new( -1.0f, 0.0f, 0.0f ), test, &clipped[0], &clipped[1]); break;
				}

				for (int w = 0; w < nTrisToAdd; w++)
					Vector_Push(&listTriangles,&clipped[w]);
			}
			nNewTriangles = listTriangles.size;
		}

		for (int j = 0;j<listTriangles.size;j++){
			Tri3D t = *(Tri3D*)Vector_Get(&listTriangles,j);
			RenderTriangle(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),t.c);
			//RenderTriangleWire(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),WHITE,1.0f);
		}

		Vector_Free(&listTriangles);
	}
	Vector_Free(&vecTrianglesToRaster);
}

void Delete(AlxWindow* w){
	Vector_Free(&meshCube.tris);
	Vector_Free(&Cubes);
}

int main(){
    if(Create("3D Test no Tex",2500,1200,1,1,Setup,Update,Delete))
        Start();
    return 0;
}