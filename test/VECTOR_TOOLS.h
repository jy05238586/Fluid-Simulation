//**************************************************************************************
// Vector and matrix library.
//**************************************************************************************

#ifndef __VECTOR_TOOLS_H__
#define __VECTOR_TOOLS_H__
#include <math.h>
#include "MY_MATH.h"


//**************************************************************************************
// Generic Functions.
//**************************************************************************************

template <class T>   T 
Norm(T *x, int number=3)	//infinite norm
{
	T ret=0;
	for(int i=0; i<number; i++)	
		if(ret<fabsf(x[i]))	ret=fabsf(x[i]);
	return ret;
}

template <class T>   T 
Dot(T *x, T *y, int number=3)
{
	T ret=0;
	for(int i=0; i<number; i++)	ret+=x[i]*y[i];
	return ret;
}

template <class T>  
void Matrix_Transpose(T *A, T *R, int nx, int ny)				//R=A'
{
	for(int i=0; i<nx; i++)
	for(int j=0; j<ny; j++)
		R[j*nx+i]=A[i*ny+j];
}

template <class T>  
void Matrix_Product(T *A, T *B, T *R, int nx, int ny, int nz)	//R=A*B
{
	memset(R, 0, sizeof(float)*nx*nz);
	for(int i=0; i<nx; i++)
	for(int j=0; j<nz; j++)
		for(int k=0; k<ny; k++)
			R[i*nz+j]+=A[i*ny+k]*B[k*nz+j];
}

template <class T>  
void Matrix_Self_Product(T *A, T *R, int nx, int ny)			//R=A'*A
{
	memset(R, 0, sizeof(T)*ny*ny);
	for(int i=0; i<ny; i++)
	for(int j=i; j<ny; j++)
	{
		for(int k=0; k<nx; k++)
			R[i*ny+j]+=A[k*ny+i]*A[k*ny+j];
		if(i!=j)	R[j*ny+i]=R[i*ny+j];		
	}
}

//**************************************************************************************
// 2D Functions.
//**************************************************************************************
template <class T>  
void Matrix_Inverse_2(T *A, T *R)
{
	T inv_det=1/(A[0]*A[3]-A[1]*A[2]);
	R[0]= A[3]*inv_det;
	R[1]=-A[1]*inv_det;
	R[2]=-A[2]*inv_det;
	R[3]= A[0]*inv_det;
}

template <class T>  
void Matrix_Transpose_2(T *A, T *R)
{
	memcpy(R, A, sizeof(T)*4);
	Swap(R[1], R[2]);
}

template <class T>  
void Matrix_Product_2(T *A, T *B, T *R)
{
	T temp_R[4];
	temp_R[0]=A[0]*B[0]+A[1]*B[2];
	temp_R[1]=A[0]*B[1]+A[1]*B[3];
	temp_R[2]=A[2]*B[0]+A[3]*B[2];
	temp_R[3]=A[2]*B[1]+A[3]*B[3];
	R[0]=temp_R[0];
	R[1]=temp_R[1];
	R[2]=temp_R[2];
	R[3]=temp_R[3];
}

void ED_2(float *G, float *w, float *V) // G-> V'w V
{
	float a=1;
	float b=-(G[0]+G[3]);
	float c=G[0]*G[3]-G[1]*G[2];
	float delta=(b*b-4*c);
	if(delta<0)	{delta=0;}
	else		delta=sqrtf(delta);

	w[0]=(-b+delta)*0.5f;
	w[1]=(-b-delta)*0.5f;
	float inv_length;
	a=G[0]-w[0];
	b=G[1];

	if(fabsf(a)<1e-6f && fabsf(b)<1e-6f)
	{
		V[0]=1;
		V[1]=0;
		V[2]=0;
		V[3]=1;
	}
	else
	{
		inv_length=1/sqrtf(a*a+b*b);
		V[0]= b*inv_length;
		V[1]=-a*inv_length;
		V[2]=-V[1];
		V[3]= V[0];
	}
	if(V[0]<0)
	{
		V[0]=-V[0];
		V[1]=-V[1];
		V[2]=-V[2];
		V[3]=-V[3];
	}	
	if(w[0]<0)	w[0]=0;
	if(w[1]<0)	w[1]=0;
}

void SVD_2(float *F, float *R, float *w, float *V) // F -> R w V
{
	//Remove the rotation component in F, and store it in R.
	float G[4];
	G[0]=F[0]*F[0]+F[2]*F[2];
	G[2]=F[0]*F[1]+F[2]*F[3];
	G[1]=G[2];
	G[3]=F[1]*F[1]+F[3]*F[3];
	ED_2(G, w, V);
	w[0]=sqrtf(w[0]);
	w[1]=sqrtf(w[1]);

	float inv_D[4];
	inv_D[1]=inv_D[2]=0;
	if(w[0]<1e-10f)	inv_D[0]=0;
	else			inv_D[0]=1/w[0];
	if(w[1]<1e-10f)	inv_D[3]=0;
	else			inv_D[3]=1/w[1];
	float Vt[4];
	Vt[0]=V[0];
	Vt[1]=V[2];
	Vt[2]=V[1];
	Vt[3]=V[3];
	Matrix_Product_2(F, Vt, R);
	Matrix_Product_2(R, inv_D, R);
}

//**************************************************************************************
// 3D Vector Functions.
//**************************************************************************************
template <class T>   
void Cross(T *x, T *y, T *r)
{
	r[0]= x[1]*y[2]-y[1]*x[2];
	r[1]=-x[0]*y[2]+y[0]*x[2];
	r[2]= x[0]*y[1]-y[0]*x[1];
}

template <class T>   
T Magnitude_Squared(T *x)
{
	return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

template <class T>   
T Magnitude(T *x)
{
	return sqrtf(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

template <class T>   
T Normalize(T *x, T *r=0)
{
	if(r==0)	r=x;
	T m=Magnitude(x);
	if(m<1e-14f)	{printf("ERROR: vector cannot be normalized.\n"); return m;}
	T inv_m=1/m;
	r[0]=x[0]*inv_m;
	r[1]=x[1]*inv_m;
	r[2]=x[2]*inv_m;
	return m;
}


template <class T>  
void Matrix_T_Product_2(T *A, T *B, T *R)
{
	for(int i=0; i<2; i++)
	for(int j=0; j<2; j++)
		R[i*2+j]=A[i]*B[j]+A[i+2]*B[j+2];	
}


//**************************************************************************************
// 3D Matrix Functions.
//**************************************************************************************
template <class T>  
T Determinant_2(T *A)
{
	return A[0]*(A[4]*A[8]-A[7]*A[5])+A[3]*(A[7]*A[2]-A[1]*A[8])+A[6]*(A[1]*A[5]-A[4]*A[2]);
}

template <class T>  
void Matrix_Add_3(T *A, T a, T *B, T b, T *R)	//R=aA+bB
{
	for(int i=0; i<9; i++)
		R[i]=A[i]*a+B[i]*b;
}

template <class T>  
void Matrix_Add_3(T *A, T *B, T *R)				//R=A+B
{
	for(int i=0; i<9; i++)
		R[i]=A[i]+B[i];
}

template <class T>  
void Matrix_Substract_3(T *A, T *B, T *R)		//R=A-B
{
	for(int i=0; i<9; i++)
		R[i]=A[i]-B[i];
}

template <class T>  
void Matrix_Inverse_3(T *A, T *R)				//R=inv(A)
{
	R[0]=A[4]*A[8]-A[7]*A[5];
	R[1]=A[7]*A[2]-A[1]*A[8];
	R[2]=A[1]*A[5]-A[4]*A[2];
	R[3]=A[5]*A[6]-A[3]*A[8];
	R[4]=A[0]*A[8]-A[2]*A[6];
	R[5]=A[2]*A[3]-A[0]*A[5];
	R[6]=A[3]*A[7]-A[4]*A[6];
	R[7]=A[1]*A[6]-A[0]*A[7];
	R[8]=A[0]*A[4]-A[1]*A[3];
	T inv_det=1/(A[0]*R[0]+A[3]*R[1]+A[6]*R[2]);
	for(int i=0; i<9; i++)
		R[i]*=inv_det;
}

template <class T>  					//R=A'
void Matrix_Transpose_3(T *A, T *R)
{
	memcpy(R, A, sizeof(T)*9);
	Swap(R[1], R[3]);
	Swap(R[2], R[6]);
	Swap(R[5], R[7]);
}

template <class T>  
void Matrix_Factorization_3(T *A, T *R)			//R=chol(A), Chelosky factorization
{
	R[0]=sqrtf(A[0]);
	R[1]=A[1]/R[0];
	R[2]=A[2]/R[0];
	R[3]=0;
	R[4]=sqrtf(A[4]-R[1]*R[1]);
	R[5]=(A[5]-R[1]*R[2])/R[4];
	R[6]=0;
	R[7]=0;
	R[8]=sqrtf(A[8]-R[2]*R[2]-R[5]*R[5]);
}

template <class T>  
void Matrix_Product_3(T *A, T *B, T *R)		//R=A*B
{
	R[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
	R[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
	R[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
	R[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
	R[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
	R[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
	R[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
	R[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
	R[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}

template <class T>  
void Matrix_T_Product_3(T *A, T *B, T *R)	//R=A'*B
{

	R[0]=A[0]*B[0]+A[3]*B[3]+A[6]*B[6];
	R[1]=A[0]*B[1]+A[3]*B[4]+A[6]*B[7];
	R[2]=A[0]*B[2]+A[3]*B[5]+A[6]*B[8];
	R[3]=A[1]*B[0]+A[4]*B[3]+A[7]*B[6];
	R[4]=A[1]*B[1]+A[4]*B[4]+A[7]*B[7];
	R[5]=A[1]*B[2]+A[4]*B[5]+A[7]*B[8];
	R[6]=A[2]*B[0]+A[5]*B[3]+A[8]*B[6];
	R[7]=A[2]*B[1]+A[5]*B[4]+A[8]*B[7];
	R[8]=A[2]*B[2]+A[5]*B[5]+A[8]*B[8];
}

template <class T>  
void Matrix_Vector_Product_3(T *A, T *x, T *r)	//r=A*x
{
	r[0]=A[0]*x[0]+A[1]*x[1]+A[2]*x[2];
	r[1]=A[3]*x[0]+A[4]*x[1]+A[5]*x[2];
	r[2]=A[6]*x[0]+A[7]*x[1]+A[8]*x[2];
}

template <class T>  
void Matrix_T_Vector_Product_3(T *A, T *x, T *r)//r=A'*x
{
	r[0]=A[0]*x[0]+A[3]*x[1]+A[6]*x[2];
	r[1]=A[1]*x[0]+A[4]*x[1]+A[7]*x[2];
	r[2]=A[2]*x[0]+A[5]*x[1]+A[8]*x[2];
}

//**************************************************************************************
// 3D Utility Functions.
//**************************************************************************************
template <class T>  
T Normal(T *p0, T *p1, T *p2, T *normal)
{
	T e0[3], e1[3];
	for(int i=0; i<3; i++)
	{
		e0[i]=p1[i]-p0[i];
		e1[i]=p2[i]-p0[i];
	}		
	Cross(e0, e1, normal);
	return Normalize(normal);
}

float Cot(float *xa, float *x, float *xb)
{
	// return the contangent value for the angle xa<-x->xb
	float edge_a[3], edge_b[3];
	edge_a[0]=xa[0]-x[0];
	edge_a[1]=xa[1]-x[1];
	edge_a[2]=xa[2]-x[2];
	edge_b[0]=xb[0]-x[0];
	edge_b[1]=xb[1]-x[1];
	edge_b[2]=xb[2]-x[2];
	Normalize(edge_a);
	Normalize(edge_b);
	float dot=Dot(edge_a, edge_b);
	return dot/sqrtf(1-dot*dot);
}

template <class T>  
T Distance(T *x, T *y)
{
	float r[3];
	r[0]=x[0]-y[0];
	r[1]=x[1]-y[1];
	r[2]=x[2]-y[2];
	return Magnitude(r);
}

template <class T>  
T Distance_Squared(T *x, T *y)
{
	float r[3];
	r[0]=x[0]-y[0];
	r[1]=x[1]-y[1];
	r[2]=x[2]-y[2];
	return Magnitude_Squared(r);}



//**************************************************************************************
// 4D Functions.
//**************************************************************************************
template <class T>   
void Matrix_Product_4(T *A, T *B, T *R) // r=a*b
{
	memset(R, 0, sizeof(T)*16);
	for(int i=0; i<4; i++)	
	for(int j=0; j<4; j++)	
	for(int k=0; k<4; k++)
		R[i*4+j]+=A[i*4+k]*B[k*4+j];
}

template <class T>   
void Matrix_Times_Vector_4(T *A, T *x, T *r)
{
	if(r==x)	printf("ERROR: r cannot be equal to x in Matrix_Times_Vector_4.\n");
	r[0]=A[ 0]*x[0]+A[ 1]*x[1]+A[ 2]*x[2]+A[ 3]*x[3];
	r[1]=A[ 4]*x[0]+A[ 5]*x[1]+A[ 6]*x[2]+A[ 7]*x[3];
	r[2]=A[ 8]*x[0]+A[ 9]*x[1]+A[10]*x[2]+A[11]*x[3];
	r[3]=A[12]*x[0]+A[13]*x[1]+A[14]*x[2]+A[15]*x[3];
}

template <class T>   
void Matrix_T_Times_Vector_4(T *A, T *x, T *r)
{
	if(r==x)	printf("ERROR: r cannot be equal to x in Matrix_Times_Vector_4.\n");
	r[0]=A[0]*x[0]+A[4]*x[1]+A[ 8]*x[2]+A[12]*x[3];
	r[1]=A[1]*x[0]+A[5]*x[1]+A[ 9]*x[2]+A[13]*x[3];
	r[2]=A[2]*x[0]+A[6]*x[1]+A[10]*x[2]+A[14]*x[3];
	r[3]=A[3]*x[0]+A[7]*x[1]+A[11]*x[2]+A[15]*x[3];
}

template <class T>   
void Homogeneous_Projection(T *a, T *r=0)
{
	if(r==0)	r=a;
	r[0]=a[0]/a[3];
	r[1]=a[1]/a[3];
	r[2]=a[2]/a[3];
}



#endif