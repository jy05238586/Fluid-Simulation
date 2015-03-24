#include "stdafx.h"
#include "stdio.h"
#include "math.h"
//This program requires the OpenGL and GLUT libraries
// You can obtain them for free from http://www.opengl.org
#include "GL/glut.h"

struct GLvector
{
	GLfloat fX;
	GLfloat fY;
	GLfloat fZ;     
};

//These tables are used so that everything can be done in little loops that you can look at all at once
// rather than in pages and pages of unrolled code.

//a2fVertexOffset lists the positions, relative to vertex0, of each of the 8 vertices of a cube
static const GLfloat a2fVertexOffset[8][3] =
{
	{0.0, 0.0, 0.0},{1.0, 0.0, 0.0},{1.0, 1.0, 0.0},{0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0},{1.0, 0.0, 1.0},{1.0, 1.0, 1.0},{0.0, 1.0, 1.0}
};

//a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
static const GLint a2iEdgeConnection[12][2] = 
{
	{0,1}, {1,2}, {2,3}, {3,0},
	{4,5}, {5,6}, {6,7}, {7,4},
	{0,4}, {1,5}, {2,6}, {3,7}
};

//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
static const GLfloat a2fEdgeDirection[12][3] =
{
	{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
	{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
	{0.0, 0.0, 1.0},{0.0, 0.0, 1.0},{ 0.0, 0.0, 1.0},{0.0,  0.0, 1.0}
};

//a2iTetrahedronEdgeConnection lists the index of the endpoint vertices for each of the 6 edges of the tetrahedron
static const GLint a2iTetrahedronEdgeConnection[6][2] =
{
	{0,1},  {1,2},  {2,0},  {0,3},  {1,3},  {2,3}
};

//a2iTetrahedronEdgeConnection lists the index of verticies from a cube 
// that made up each of the six tetrahedrons within the cube
static const GLint a2iTetrahedronsInACube[6][4] =
{
	{0,5,1,6},
	{0,1,2,6},
	{0,2,3,6},
	{0,3,7,6},
	{0,7,4,6},
	{0,4,5,6},
};

static const GLfloat afAmbientWhite [] = {0.25, 0.25, 0.25, 1.00}; 
static const GLfloat afAmbientRed   [] = {0.25, 0.00, 0.00, 1.00}; 
static const GLfloat afAmbientGreen [] = {0.00, 0.25, 0.00, 1.00}; 
static const GLfloat afAmbientBlue  [] = {0.00, 0.00, 0.25, 1.00}; 
static const GLfloat afDiffuseWhite [] = {0.75, 0.75, 0.75, 1.00}; 
static const GLfloat afDiffuseRed   [] = {0.75, 0.00, 0.00, 1.00}; 
static const GLfloat afDiffuseGreen [] = {0.00, 0.75, 0.00, 1.00}; 
static const GLfloat afDiffuseBlue  [] = {0.00, 0.00, 0.75, 1.00}; 
static const GLfloat afSpecularWhite[] = {1.00, 1.00, 1.00, 1.00}; 
static const GLfloat afSpecularRed  [] = {1.00, 0.25, 0.25, 1.00}; 
static const GLfloat afSpecularGreen[] = {0.25, 1.00, 0.25, 1.00}; 
static const GLfloat afSpecularBlue [] = {0.25, 0.25, 1.00, 1.00}; 


GLenum    ePolygonMode = GL_FILL;
GLint     iDataSetSize = 16;
GLfloat   fStepSize = 1.0/iDataSetSize;
GLfloat   fTargetValue = 48.0;
GLfloat   fTime = 0.0;
GLvector  sSourcePoint[3];
GLboolean bSpin = true;
GLboolean bMove = true;
GLboolean bLight = true;


GLfloat fGetOffset(GLfloat fValue1, GLfloat fValue2, GLfloat fValueDesired)
{
	GLdouble fDelta = fValue2 - fValue1;

	if(fDelta == 0.0)
	{
		return 0.5;
	}
	return (fValueDesired - fValue1)/fDelta;
}

GLfloat fSample(GLfloat fX, GLfloat fY, GLfloat fZ)
{
	GLdouble fResult = 0.0;
	GLdouble fDx, fDy, fDz;
	fDx = fX - sSourcePoint[0].fX;
	fDy = fY - sSourcePoint[0].fY;
	fDz = fZ - sSourcePoint[0].fZ;
	fResult += 1/(fDx*fDx + fDy*fDy + fDz*fDz);

	fDx = fX - sSourcePoint[1].fX;
	fDy = fY - sSourcePoint[1].fY;
	fDz = fZ - sSourcePoint[1].fZ;
	//fResult += 1.0/(fDx*fDx + fDy*fDy + fDz*fDz);

	fDx = fX - sSourcePoint[2].fX;
	fDy = fY - sSourcePoint[2].fY;
	fDz = fZ - sSourcePoint[2].fZ;
	//fResult += 1.5/(fDx*fDx + fDy*fDy + fDz*fDz);

	return fResult;
}

GLvoid vNormalizeVector(GLvector &rfVectorResult, GLvector &rfVectorSource)
{
	GLfloat fOldLength;
	GLfloat fScale;

	fOldLength = sqrtf( (rfVectorSource.fX * rfVectorSource.fX) +
		(rfVectorSource.fY * rfVectorSource.fY) +
		(rfVectorSource.fZ * rfVectorSource.fZ) );

	if(fOldLength == 0.0)
	{
		rfVectorResult.fX = rfVectorSource.fX;
		rfVectorResult.fY = rfVectorSource.fY;
		rfVectorResult.fZ = rfVectorSource.fZ;
	}
	else
	{
		fScale = 1.0/fOldLength;
		rfVectorResult.fX = rfVectorSource.fX*fScale;
		rfVectorResult.fY = rfVectorSource.fY*fScale;
		rfVectorResult.fZ = rfVectorSource.fZ*fScale;
	}
}

GLvoid vGetNormal(GLvector &rfNormal, GLfloat fX, GLfloat fY, GLfloat fZ)
{
	rfNormal.fX = fSample(fX-0.01, fY, fZ) - fSample(fX+0.01, fY, fZ);
	rfNormal.fY = fSample(fX, fY-0.01, fZ) - fSample(fX, fY+0.01, fZ);
	rfNormal.fZ = fSample(fX, fY, fZ-0.01) - fSample(fX, fY, fZ+0.01);
	vNormalizeVector(rfNormal, rfNormal);
}

GLvoid vMarchCube1(GLfloat fX, GLfloat fY, GLfloat fZ, GLfloat fScale)
{
	extern GLint aiCubeEdgeFlags[256];
	extern GLint a2iTriangleConnectionTable[256][16];

	GLint iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
	GLfloat fOffset;
	GLvector sColor;
	GLfloat afCubeValue[8];
	GLvector asEdgeVertex[12];
	GLvector asEdgeNorm[12];

	//Make a local copy of the values at the cube's corners
	for(iVertex = 0; iVertex < 8; iVertex++)
	{
		afCubeValue[iVertex] = fSample(fX + a2fVertexOffset[iVertex][0]*fScale,
			fY + a2fVertexOffset[iVertex][1]*fScale,
			fZ + a2fVertexOffset[iVertex][2]*fScale);
	}

	//Find which vertices are inside of the surface and which are outside
	iFlagIndex = 0;
	for(iVertexTest = 0; iVertexTest < 8; iVertexTest++)
	{
		if(afCubeValue[iVertexTest] <= fTargetValue) 
			iFlagIndex |= 1<<iVertexTest;
	}

	//Find which edges are intersected by the surface
	iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

	//If the cube is entirely inside or outside of the surface, then there will be no intersections
	if(iEdgeFlags == 0) 
	{
		return;
	}

	//Find the point of intersection of the surface with each edge
	//Then find the normal to the surface at those points
	for(iEdge = 0; iEdge < 12; iEdge++)
	{
		//if there is an intersection on this edge
		if(iEdgeFlags & (1<<iEdge))
		{
			fOffset = fGetOffset(afCubeValue[ a2iEdgeConnection[iEdge][0] ], 
				afCubeValue[ a2iEdgeConnection[iEdge][1] ], fTargetValue);

			asEdgeVertex[iEdge].fX = fX + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][0]  +  fOffset * a2fEdgeDirection[iEdge][0]) * fScale;
			asEdgeVertex[iEdge].fY = fY + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][1]  +  fOffset * a2fEdgeDirection[iEdge][1]) * fScale;
			asEdgeVertex[iEdge].fZ = fZ + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][2]  +  fOffset * a2fEdgeDirection[iEdge][2]) * fScale;

			vGetNormal(asEdgeNorm[iEdge], asEdgeVertex[iEdge].fX, asEdgeVertex[iEdge].fY, asEdgeVertex[iEdge].fZ);
		}
	}


	//Draw the triangles that were found.  There can be up to five per cube
	for(iTriangle = 0; iTriangle < 5; iTriangle++)
	{
		if(a2iTriangleConnectionTable[iFlagIndex][3*iTriangle] < 0)
			break;

		for(iCorner = 0; iCorner < 3; iCorner++)
		{
			iVertex = a2iTriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];

			glColor3f(sColor.fX, sColor.fY, sColor.fZ);
			glNormal3f(asEdgeNorm[iVertex].fX,   asEdgeNorm[iVertex].fY,   asEdgeNorm[iVertex].fZ);
			glVertex3f(asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
		}
	}
}

GLvoid vMarchingCubes()
{
	GLint iX, iY, iZ;
	for(iX = 0; iX < iDataSetSize; iX++)
		for(iY = 0; iY < iDataSetSize; iY++)
			for(iZ = 0; iZ < iDataSetSize; iZ++)
			{
				vMarchCube1(iX*fStepSize, iY*fStepSize, iZ*fStepSize, fStepSize);
			}
}