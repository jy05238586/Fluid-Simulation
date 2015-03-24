#ifndef __FLuid_H__
#define __FLuid_H__
#include<stdlib.h>
#include<time.h>
#include<array>
#include<vector>
#include<gl/glut.h>
#define MY_PI		3.14159265358979323846f
using namespace std;

int Clamp(int val, int low, int high) 
{
	if(val<low) return low;
	if(val>high) return high;
	return val;
}

struct Vect {
public:
	// Vect Public Methods
	Vect() { x = y = z = 0.f; }
	Vect(float xx, float yy, float zz)
		: x(xx), y(yy), z(zz) {
	}

	Vect operator+(const Vect &v) const {
		return Vect(x + v.x, y + v.y, z + v.z);
	}

	Vect& operator+=(const Vect &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Vect operator-(const Vect &v) const {
		return Vect(x - v.x, y - v.y, z - v.z);
	}

	Vect& operator-=(const Vect &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	Vect operator*(float f) const { return Vect(f*x, f*y, f*z); }

	Vect &operator*=(float f) {
		x *= f; y *= f; z *= f;
		return *this;
	}
	Vect operator/(float f) const {
		float inv = 1.f / f;
		return Vect(x * inv, y * inv, z * inv);
	}

	Vect &operator/=(float f) {
		float inv = 1.f / f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	Vect operator-() const { return Vect(-x, -y, -z); }
	float operator[](int i) const {
		return (&x)[i];
	}

	float &operator[](int i) {
		return (&x)[i];
	}
	float LengthSquared() const {
		return x*x + y*y + z*z; 
	}
	float Length() const { return sqrtf(LengthSquared()); }

	bool operator==(const Vect &v) const {
		return x == v.x && y == v.y && z == v.z;
	}
	bool operator!=(const Vect &v) const {
		return x != v.x || y != v.y || z != v.z;
	}

	// Vect Public Data
	float x, y, z;
};

struct Fluid_particle
{
public:
	Vect x;			//	position
	Vect v;			//	velocity
	Vect f;
	float rho;		//	density
	float m;		//  mass
	float miu;
	float rho0;
	float ci;
	float cs;
	float nl;
	Fluid_particle(float a, float b, float c) {
		v.x = v.y = v.z = 0.f;
		x.x = a; x.y = b; x.z = c;
		rho = 0;
	}
	~Fluid_particle() {}

};

struct Voxel {
	// Voxel Public Methods
	int size() const { return particles.size(); }
	Voxel() { }
	Voxel(Fluid_particle op) {
		particles.push_back(op);
	}
	void AddParticles(Fluid_particle part) {
		particles.push_back(part);
	}

	vector<Fluid_particle> particles;
};

struct BBox {
	int MaximumExtent() const {
		Vect diag = pMax - pMin;
		if (diag.x > diag.y && diag.x > diag.z)
			return 0;
		else if (diag.y > diag.z)
			return 1;
		else
			return 2;
	}

	Vect pMin, pMax;
};

class GridAccel {
public:
	// GridAccel Private Data
	vector<Fluid_particle> particles;
	int nVoxels[3];
	Vect width, invWidth;
	Voxel *voxels[100000000];
	BBox bounds;

	//const values
	float h;
	float m0;
	float rho0;
	float rho_air;
	float k;
	float miu;
	Vect f_g;
	float b;
	float dt;
	float threshold_p;
	float threshold_n;
	float threshold_s;
	float threshold_rho;
	float air_count;
	float threshold_aircount;
	float d;
	float sigma;
	float box_size;
	float ground;
	float test;
	float k_p6;
	float g_k_p6;
	float l_k_p6;
	float l_k_v;
	float g_k_s;

	

	//for rho & color field
	float Kernel_poly6(Vect r, float h) {
		return (k_p6 * powf((h*h - r.LengthSquared() ), 3) );
	}

	Vect Grad_Kernel_poly6(Vect r, float h){
		Vect v;
		float lr = h*h-r.LengthSquared();
		lr = lr*lr;
		for (int axis = 0; axis < 3; axis++)
		{
			v[axis] =  g_k_p6 * lr * r[axis];
		}
		return v;
	}

	float Lapl_Kernel_poly6(Vect r, float h) {
		float rls = r.LengthSquared();
		float hs_rs = h*h - rls;
		return l_k_p6*(hs_rs * rls - hs_rs * hs_rs);
	}

	//for viscosity
	float Lapl_Kernel_vis(Vect r, float h){
		return (l_k_v * (h - r.Length()));
	}

	//for pressure
	Vect Grad_Kernel_spiky(Vect r, float h){
		Vect v;
		float lr = r.Length();
		if (lr != 0) {
			for (int axis = 0; axis < 3; axis++)
			{
				v[axis] =  g_k_s * powf(h-lr, 2) / lr * r[axis];
			}
		}
		return v;
	}


	int posToVoxel(const Vect P, int axis) const {
		int vv = int((P[axis] - bounds.pMin[axis]) *
			invWidth[axis]);
		return Clamp(vv, 0, nVoxels[axis]-1);
	}
	float voxelToPos(int p, int axis) const {
		return bounds.pMin[axis] + p * width[axis];
	}
	inline int offset(int x, int y, int z) const {
		return z*nVoxels[0]*nVoxels[1] + y*nVoxels[0] + x;
	}

	void render();

	void solve();

	Vect solve_fpressure(int, int, bool&);

	void solve_rho(int o);

	void moveparticle(int o, int i);

	void refine_grid();

	GridAccel();

	~GridAccel() {}
};

GridAccel:: GridAccel() {
	//init particles& bounds
	bounds.pMax[0] = 1;
	bounds.pMax[1] = 1;
	bounds.pMax[2] = 1;
	bounds.pMin[0] = -1;
	bounds.pMin[1] = -1;
	bounds.pMin[2] = -1;
	h = 0.045f;
	k = 10;
	miu = 0.5;
	b = 5;
	f_g = Vect(0, -9.8f, 0);
	dt = 0.004f;	
	threshold_n = 0.1f;
	sigma = 0.6f;
	threshold_p = 35;
	threshold_s = 1;
	threshold_rho = 50;
	rho0 = 1000;
	rho_air = 100;
	air_count = 0;
	threshold_aircount = 1000;
	d = 0.005;
	

	ground = -0.2f;
	box_size = 0.1;
	k_p6 = 315 / (64 * MY_PI * powf(h, 9));
	g_k_p6 = (0-6)*k_p6;
	l_k_p6 = 24 * k_p6;
	l_k_v = 45/(MY_PI * powf(h, 6));
	g_k_s = 0.f - 45 / (MY_PI * powf(h, 6));
	


	float range_max = 0.1;
	int range_min = -0.1;
	srand(time(0));
	for (int i = 0; i < 1500; i++) {
		float x = ((double)rand() / (RAND_MAX + 1) * (range_max - range_min) + range_min);
		float z = ((double)rand() / (RAND_MAX + 1) * (range_max - range_min) + range_min);
		float y = ((double)rand() / (RAND_MAX + 1) * (range_max*2 - range_min*2) + range_min*2);
		y = y +0.2;
		Fluid_particle p = Fluid_particle(x,y,z);
		p.nl = 0;
		p.m = 0.006f;
		p.cs = 1;
		p.ci = -0.5;
		p.rho0 = 500.f;
		p.miu = miu;
		particles.push_back(p);
	}
	for (int i = 0; i < 1500; i++) {
		float x = ((double)rand() / (RAND_MAX + 1) * (range_max - range_min) + range_min);
		float z = ((double)rand() / (RAND_MAX + 1) * (range_max - range_min) + range_min);
		float y = ((double)rand() / (RAND_MAX + 1) * (range_max*2 - range_min*2) + range_min*2);
		y = y +0.4;
		Fluid_particle p = Fluid_particle(x,y,z);
		p.nl = 0;
		p.m = 0.012f;
		p.cs = 1;
		p.ci = 0.5;
		p.miu = miu;
		p.rho0 = 1000.f;
		particles.push_back(p);
	}
	
	//---------------------------------------------------------------------

	Vect delta = bounds.pMax - bounds.pMin;

	// Find _voxelsPerUnitDist_ for grid
	int maxAxis = bounds.MaximumExtent();
	float invMaxWidth = 1.f / delta[maxAxis];
	//float cubeRoot = 3.f * powf(float(particles.size()), 1.f/3.f);		//
	float cubeRoot = 2.f / h;
	float voxelsPerUnitDist = cubeRoot * invMaxWidth;					//
	for (int axis = 0; axis < 3; ++axis) {
		nVoxels[axis] = (int)floorf(delta[axis] * voxelsPerUnitDist + 0.5f); //RoundtoInt
	}

	// Compute voxel widths and allocate voxels
	for (int axis = 0; axis < 3; ++axis) {
		width[axis] = delta[axis] / nVoxels[axis];
		invWidth[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
	}
	int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];

	// Add particles to grid voxels
	for (int i = 0; i < particles.size(); ++i) {
		int voxelpos[3];
		for (int axis = 0; axis < 3; ++axis) {
			voxelpos[axis] = posToVoxel(particles[i].x, axis);
		}
		// Add particles to voxels
		int o = offset(voxelpos[0],voxelpos[1],voxelpos[2]);
		if (!voxels[o]) {
			// Allocate new voxel and store particles in it
			voxels[o] = new Voxel(particles[i]);
		}
		else {
			// Add primitive to already-allocated voxel
			voxels[o]->AddParticles(particles[i]);
		}				
	}
}

void GridAccel::render() {
	int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];
	float diffuse_color[3]={0, 0, 0};
		

	float x, y ,z;
	float rend_color;
	float b_min = 90;
	float b_max = 2100;

	float gap = (b_max - b_min)/3;

	for (int i = 0; i < nv; i++) {
		if (voxels[i])
		{
			for (int j = 0 ; j < voxels[i]->particles.size(); j++) {
				x = voxels[i]->particles[j].x.x;
				y = voxels[i]->particles[j].x.y;
				z = voxels[i]->particles[j].x.z;
				diffuse_color[0] = 0;
				diffuse_color[1] = 0;
				diffuse_color[2] = 0;
				
				rend_color = voxels[i]->particles[j].rho0;
				if (rend_color < b_min+gap) {
					diffuse_color[2] = (rend_color - b_min)/gap;
					diffuse_color[2] = 1;
				}
				else if (rend_color < b_max-gap) {
					diffuse_color[1] = (rend_color - b_min-gap)/gap;
					diffuse_color[1] = 1;
				}
				else {
					diffuse_color[0] = (rend_color - b_max+gap)/gap;
					diffuse_color[0] = 1;
				}
				glPushMatrix();
				glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_color);
				glTranslated(x, y, z);
				glutSolidSphere(0.005, 5, 5);
				glPopMatrix();
			}
		}
	}

	float diffuse_color2[3]={1, 1, 1};
	glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_color2);

	glBegin(GL_POLYGON);
	glNormal3f(0,1,0);
	glVertex3f(-1,ground,-1);
	glVertex3f( 1,ground,-1);
	glVertex3f( 1,ground, 1);
	glVertex3f(-1,ground, 1);
	glEnd();	
}

void GridAccel::solve_rho(int o) {
	Vect dr;
	int o2 = o;
	int z = (int)(o2 / (nVoxels[0]*nVoxels[1]));
	o2 = o2 - z * nVoxels[0]*nVoxels[1];
	int y = (int)(o2/nVoxels[0]);
	int x = o2 - y * nVoxels[0];
	if (voxels[o])
		for (int i = 0; i < voxels[o]->size(); i++) {
			voxels[o]->particles[i].rho = 0;
			for (int ix = max(0, x - 1); ix < min (nVoxels[0], x+2); ix++)
				for (int iy = max(0, y - 1); iy < min (nVoxels[1], y+2); iy++)					
					for (int iz = max(0, z - 1); iz < min (nVoxels[2], z+2); iz++)
					{
						o2 = offset(ix, iy, iz);
						if (voxels[o2]) {
							for (int j = 0; j < voxels[o2]->size(); j++){
								dr = voxels[o]->particles[i].x - voxels[o2]->particles[j].x;
								if (h > dr.Length()){
									voxels[o]->particles[i].rho += voxels[o2]->particles[j].m * Kernel_poly6(dr, h);
								}
							}
						}
					}
			if (voxels[o]->particles[i].rho0 == rho_air && voxels[o]->particles[i].rho < threshold_rho) {
				vector<Fluid_particle>::iterator itr = voxels[o]->particles.begin(); 
				voxels[o]->particles.erase(itr+i);
				air_count--;
				i--;
			}
		}
}

Vect GridAccel::solve_fpressure(int o, int i, bool &flag) {
	Vect f;
	Vect f_v;
	Vect dr;
	Vect n;
	Vect ni;
	Vect ns;
	Vect np;
	float cs = 0;
	float lap_c = 0;
	float lap_ci = 0;
	float lap_cs = 0;
	int o2 = o;
	int z = (int)(o2 / (nVoxels[0]*nVoxels[1]));
	o2 = o2 - z * nVoxels[0]*nVoxels[1];
	int y = (int)(o2/nVoxels[0]);
	int x = o2 - y * nVoxels[0];

	if (voxels[o]) {
		for (int xx = max(0, x - 1); xx < min (nVoxels[0], x+2); xx++)
			for (int yy = max(0, y - 1); yy < min (nVoxels[1], y+2); yy++)					
				for (int zz = max(0, z - 1); zz < min (nVoxels[2], z+2); zz++)
				{
					o2 = offset(xx, yy, zz);
					if (voxels[o2])
						for (int j = 0; j < voxels[o2]->size(); j++) {
							dr = voxels[o]->particles[i].x - voxels[o2]->particles[j].x;
							if (h > dr.Length())
							{
								float m = voxels[o2]->particles[j].m;
								float pi = k * (voxels[o]->particles[i].rho - voxels[o]->particles[i].rho0);
								float pj = k * (voxels[o2]->particles[j].rho - voxels[o2]->particles[j].rho0);
								float rhoj_2 = 2 * voxels[o2]->particles[j].rho;
								f -= Grad_Kernel_spiky(dr, h) * (m * (pi+pj)/ rhoj_2);
								f_v += (voxels[o2]->particles[j].v - voxels[o]->particles[i].v) * m /voxels[o2]->particles[j].rho * Lapl_Kernel_vis(dr, h) * (voxels[o2]->particles[j].miu + voxels[o]->particles[i].miu)/2;
								n = Grad_Kernel_poly6(dr,h)*m/voxels[o2]->particles[j].rho;
								lap_c = Lapl_Kernel_poly6(dr, h)*m/voxels[o2]->particles[j].rho;
								cs += Kernel_poly6(dr, h)*m/voxels[o2]->particles[j].rho * voxels[o2]->particles[j].cs;
								ns += n * voxels[o2]->particles[j].cs;
								ni += n * voxels[o2]->particles[j].ci;
								np += n;
								lap_cs += lap_c * voxels[o2]->particles[j].cs;
								lap_ci += lap_c * voxels[o2]->particles[j].ci;								
							}							
						}
				}
				voxels[o]->particles[i].nl = cs;
				if (voxels[o]->particles[i].ci == 0)
				{
					if ((np.Length() > threshold_p && ns.Length() < threshold_s)) {
						vector<Fluid_particle>::iterator itr = voxels[o]->particles.begin(); 
						voxels[o]->particles.erase(itr+i);
						air_count--;
						flag = true;
					}
					else if (voxels[o]->particles[i].rho > rho0/3*2) {
						f += f_g * b * (voxels[o]->particles[i].rho0 - rho0);	
					}
				}
				if (flag == false) {			
					if (ns.Length() > threshold_n) {
						f -= ns * sigma * lap_cs / ns.Length();
					}
					if (ni.Length() > threshold_n) {
						f -= ni * sigma * lap_ci / ni.Length();
					}

					//for air particles add additional buoyancy force
						
					/*
					if (np.Length() > threshold_p && 
						np.y > 0 && 
						air_count < threshold_aircount &&
						voxels[o]->particles[i].rho < rho0 &&
						voxels[o]->particles[i].ci != 0) {
						//generate air particles
						Vect pos_x;
						pos_x = voxels[o]->particles[i].x - np * d / np.Length();
						Fluid_particle p = Fluid_particle(pos_x.x, pos_x.y, pos_x.z);
						p.v = voxels[o]->particles[i].v;
						p.nl = 0;
						p.m = 0.0012f;
						p.cs = 0;
						p.ci = 0;
						p.miu = 0;
						p.rho0 = 100.f;
						p.rho = 100.f;
						voxels[o]->particles.push_back(p);
						air_count++;
					}
					*/
				}
	}

	return  f_v + f;
}

void GridAccel::refine_grid() {
	int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];
	for (int i = 0 ; i < nv; i++) {
		if (voxels[i]) {
			for (int j = 0; j < voxels[i]->size(); j++) {
				int voxelpos[3];
				for (int axis = 0; axis < 3; ++axis) {
					voxelpos[axis] = posToVoxel(voxels[i]->particles[j].x, axis);
				}
				int o = offset(voxelpos[0],voxelpos[1],voxelpos[2]);
				if ( o != i) {
					if (!voxels[o]) {
						// Allocate new voxel and store particles in it
						voxels[o] = new Voxel(voxels[i]->particles[j]);
					}
					else {
						// Add primitive to already-allocated voxel
						voxels[o]->AddParticles(voxels[i]->particles[j]);
					}

					voxels[i]->particles.erase(voxels[i]->particles.begin()+j);
					j--;
				}
			}
		}
	}
}

void GridAccel::moveparticle(int o, int i) {
	//box bound
	Fluid_particle p = Fluid_particle(0,0,0);
	p.x = voxels[o]->particles[i].x;
	p.v = voxels[o]->particles[i].v;
	p.x = p.x + p.v * dt;

	if (p.x[0] > box_size) {
		p.x[0] = box_size * 2 - p.x[0];
		p.v[0] = 0.f - p.v[0];
		p.v = p.v * 0.5f;
	}

	if (p.x[0] < 0-box_size) {
		p.x[0] = 0.f - box_size * 2 - p.x[0];
		p.v[0] = 0.f - p.v[0];
		p.v = p.v * 0.5f;
	}

	if (p.x[2] > box_size) {
		p.x[2] = box_size * 2 - p.x[2];
		p.v[2] = 0.f - p.v[2];
		p.v = p.v * 0.5f;
	}

	if (p.x[2] < 0-box_size) {
		p.x[2] = 0- box_size*2 - p.x[2];
		p.v[2] = 0.f - p.v[2];
		p.v = p.v * 0.5f;
	}

	if (p.x[1] < ground) {
		p.x[1] = 2 * ground - p.x[1];
		p.v[1] = 0.f - p.v[1];
		p.v = p.v * 0.5f;
	}

	voxels[o]->particles[i].x = p.x;
	voxels[o]->particles[i].v = p.v;
 
}

void GridAccel::solve() {
	int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];
	Vect f;
	for (int i = 0; i < nv; i++) {		
		if (voxels[i])
		{
			solve_rho(i);
		}
	}


	//for each particle solve F
		for (int i = 0; i < nv; i++) {
		if (voxels[i]){
			for (int j = 0 ; j < voxels[i]->particles.size(); j++) {
				f.x = 0;
				f.y = 0;
				f.z = 0;
				bool flag = false;
				f = solve_fpressure(i, j, flag);
				if (flag) {
					j--;
				}
				else {
					f += f_g * voxels[i]->particles[j].rho;
					voxels[i]->particles[j].f = f;
				}
			}
		}
	}

	for (int i = 0; i < nv; i++) {
		if (voxels[i]){
			for (int j = 0; j < voxels[i]->particles.size(); j++) {
				voxels[i]->particles[j].v += voxels[i]->particles[j].f / voxels[i]->particles[j].rho * dt;
				moveparticle(i, j);
			}
		}
	}

	refine_grid();
}


#endif