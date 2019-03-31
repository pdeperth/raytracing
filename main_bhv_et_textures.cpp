#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <random>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <mutex>
#include <thread>

std::mutex mtex;

using namespace std;


thread_local default_random_engine seed;
thread_local uniform_real_distribution<double> uniforme(0,1);

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) :x(x), y(y), z(z) {};
	double norm2() { return x*x + y*y + z*z; }
	void normalize() {
		double n = sqrt(norm2());
		x /= n;
		y /= n;
		z /= n;
	}
	double x, y, z;
};

Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vector operator-(const Vector& a, const Vector &b) {
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector operator-(const Vector& a) {
  return Vector(-a.x, -a.y, -a.z);
}

Vector operator*(const double& a, const Vector &b) {
	return Vector(a * b.x, a * b.y, a * b.z);
}
Vector operator*(const Vector& a, const double &b) {
	return b * a;
}
Vector operator/(const Vector& a, const double &b) {
	return Vector(a.x / b, a.y / b, a.z / b);
}
double dot(const Vector& a, const Vector& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
Vector operator*(const Vector& a, const Vector& b){
	return Vector(a.x*b.x, a.y*b.y, a.z*b.z);
}
Vector pd_vect(const Vector& a, const Vector& b) {
	return Vector(a.y*b.z-b.y*a.z, -a.x*b.z+a.z*b.x, a.x*b.y-a.y*b.x);
}

class Ray {
public:
	Ray(const Vector& C, const Vector& u) : C(C), u(u) {};
	Vector u, C;
};

class Objet {
public:
	Objet(const Vector& albedo = Vector(0, 0, 0), const bool& miroir = false, const bool &transparent = false, const double& n = 1) : albedo(albedo), miroir(miroir), transparent(transparent), n(n){};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N) = 0;
	Vector albedo;
	bool miroir;
	bool transparent;
	double n;
};

class Triangle {
public:
	Triangle(const std::vector<Vector> &sommets, const Vector& albedo = Vector(0, 0, 0), const bool& miroir = false, const bool &transparent = false, const double& n = 1) : sommets(sommets) {};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, double &beta, double &gamma, double &t) {
		Vector A, B, C;
		A = sommets[0];
		B = sommets[1];
		C = sommets[2];
		// std::cout << A.x << " " << A.y << " " << A.z << std::endl;
		N = pd_vect(B-A, C-A);
		if (dot(r.u, N) != 0 && (dot(A-r.C, N)/dot(r.u, N) >= 0)) {
			t = dot(A-r.C, N)/dot(r.u, N);
			P = r.C + t*r.u;
			Vector AP, AB, AC;
			AP = P-A;
			AB = B-A;
			AC = C-A;
			// double beta, gamma;
			beta = (dot(AP, AB)*AC.norm2() - dot(AC, AB)*dot(AP, AC))/(AB.norm2()*AC.norm2() - pow(dot(AB, AC), 2));
			if (0 <= beta && beta <= 1) {
					gamma = (dot(AP, AC)*AB.norm2() - dot(AC, AB)*dot(AP, AB))/(AB.norm2()*AC.norm2() - pow(dot(AB, AC), 2));
					if (0 <= gamma && gamma <= 1 && 0 <= 1-beta-gamma && 1-beta-gamma <= 1){
						// std::cout << "triangle" << std::endl;
						return true;
					} else {
						return false;
					}
			} else {
				return false;
			}

		} else {
			return false;
		}
		return false;
	};

	std::vector<Vector> sommets;

};


class TriangleIndices {
public:
	 TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	 };
	 int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	 int uvi, uvj, uvk;  // indices within the uv coordinates array
	 int ni, nj, nk;  // indices within the normals array
	 int group;       // face group
};

class Bbox {
public:
	Bbox() {
		vmin = Vector(1E9, 1E9, 1E9);
		vmax = Vector(-1E9, -1E9, -1E9);
	}

	bool inter(const Ray& r) {
		double t1x = (vmin.x - r.C.x) / r.u.x;
		double t2x = (vmax.x - r.C.x) / r.u.x;
		double tminx = std::min(t1x, t2x);
		double tmaxx = std::max(t1x, t2x);
		double t1y = (vmin.y - r.C.y) / r.u.y;
		double t2y = (vmax.y - r.C.y) / r.u.y;
		double tminy = std::min(t1y, t2y);
		double tmaxy = std::max(t1y, t2y);
		double t1z = (vmin.z - r.C.z) / r.u.z;
		double t2z = (vmax.z - r.C.z) / r.u.z;
		double tminz = std::min(t1z, t2z);
		double tmaxz = std::max(t1z, t2z);
		return (std::min(std::min(tmaxx, tmaxy), tmaxz) >= std::max(std::max(tminx, tminy), tminz)) ;
	}
	Vector vmin, vmax;
};

class BHVNode {
public:
	BHVNode(Bbox b = Bbox(), int debut = 0, int fin = 0) {
		this->fd = NULL;
		this->fg = NULL;
		this->b = b;
		this->debut = debut;
		this->fin = fin;
	};
	BHVNode* fg;
	BHVNode* fd;
	Bbox b;
	int debut, fin;
	// void setlimites(const int &d, const int &f) {
	// 	debut = d;
	// 	fin = f;
	// }
};

class Geometry : public Objet {
public:
 ~Geometry() {};
 	// échanger les y et les z
	 Geometry(const char* filename = "BeautifulGirl_original.obj", const double &echelle = 10, const Vector & translation = Vector(0,0,0), const Vector& albedo = Vector(0, 0, 0), const bool& miroir = false, const bool &transparent = false, const double& n = 1) {
		this->albedo = albedo;
 		this-> miroir = miroir;
 		this->transparent = transparent;
 		this->n = n;
		readOBJ(filename);
		for (int i=0; i < indices.size(); i++) {
			vertices[i] = echelle*(Vector(vertices[i].x, vertices[i].z, vertices[i].y) + translation);
		}
		this->enveloppe_ext = build_box(0,indices.size());
		this->BHV = BHVNode(enveloppe_ext, 0, indices.size());
		// BHVNode* enveloppe;
		cout << "about to build BHV box" << endl;
		build_bhv(&BHV, 0, indices.size());
		cout << "BHV box built" << endl;
	};

	 void readOBJ(const char* obj) {

			 // char matfile[255];
			 char grp[255];

			 FILE* f;
			 f = fopen(obj, "r");
			 if (!f) {
				 std::cout << "error: can't read OBJ file." << std::endl;
			 }
			 int curGroup = -1;
			 while (!feof(f)) {
					 char line[255];
					 if (!fgets(line, 255, f)) break;

					 std::string linetrim(line);
					 linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
					 strcpy(line, linetrim.c_str());

					 if (line[0] == 'u' && line[1] == 's') {
							 sscanf(line, "usemtl %[^\n]\n", grp);
							 curGroup++;
					 }

					 if (line[0] == 'v' && line[1] == ' ') {
							 Vector vec;

							 Vector col;
							 if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.y, &vec.z, &col.x, &col.y, &col.z) == 6) {
									 col.x = std::min(1., std::max(0., col.x));
									 col.y = std::min(1., std::max(0., col.y));
									 col.z = std::min(1., std::max(0., col.z));

									 vertices.push_back(vec);
									 vertexcolors.push_back(col);

							 } else {
									 sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
									 vertices.push_back(vec);
							 }
					 }
					 if (line[0] == 'v' && line[1] == 'n') {
							 Vector vec;
							 sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
							 normals.push_back(vec);
					 }
					 if (line[0] == 'v' && line[1] == 't') {
							 Vector vec;
							 sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
							 uvs.push_back(vec);
					 }
					 if (line[0] == 'f') {
							 TriangleIndices t;
							 int i0, i1, i2, i3;
							 int j0, j1, j2, j3;
							 int k0, k1, k2, k3;
							 int nn;
							 t.group = curGroup;

							 char* consumedline = line + 1;
							 int offset;

							 nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
							 if (nn == 9) {
									 if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
									 if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
									 if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
									 if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
									 if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
									 if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
									 if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
									 if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
									 if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
									 indices.push_back(t);
							 } else {
									 nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
									 if (nn == 6) {
											 if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
											 if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
											 if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
											 if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
											 if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
											 if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
											 indices.push_back(t);
									 } else {
											 nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
											 if (nn == 3) {
													 if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
													 if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
													 if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
													 indices.push_back(t);
											 } else {
													 nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
													 if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
													 if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
													 if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
													 if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
													 if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
													 if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
													 indices.push_back(t);
											 }
									 }
							 }

							 consumedline = consumedline + offset;

							 while (true) {
									 if (consumedline[0] == '\n') break;
									 if (consumedline[0] == '\0') break;
									 nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
									 TriangleIndices t2;
									 t2.group = curGroup;
									 if (nn == 3) {
											 if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
											 if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
											 if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
											 if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
											 if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
											 if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
											 if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
											 if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
											 if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
											 indices.push_back(t2);
											 consumedline = consumedline + offset;
											 i2 = i3;
											 j2 = j3;
											 k2 = k3;
									 } else {
											 nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
											 if (nn == 2) {
													 if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
													 if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
													 if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
													 if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
													 if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
													 if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
													 consumedline = consumedline + offset;
													 i2 = i3;
													 j2 = j3;
													 indices.push_back(t2);
											 } else {
													 nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
													 if (nn == 2) {
															 if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
															 if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
															 if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
															 if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
															 if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
															 if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
															 consumedline = consumedline + offset;
															 i2 = i3;
															 k2 = k3;
															 indices.push_back(t2);
													 } else {
															 nn = sscanf(consumedline, "%u%n", &i3, &offset);
															 if (nn == 1) {
																	 if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
																	 if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
																	 if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
																	 consumedline = consumedline + offset;
																	 i2 = i3;
																	 indices.push_back(t2);
															 } else {
																	 consumedline = consumedline + 1;
															 }
													 }
											 }
									 }
							 }

					 }

			 }
			 fclose(f);
			 // std::cout << indices.size() << std::endl;
			 // std::cout << vertices[0].x << " " << vertices[0].y << " " << vertices[0].z << std::endl;

	 }

	 Bbox build_box(const int &i0, const int &i1) {
 		Bbox box;
 		for (int i = i0; i < i1; i++) {
			cout << "build_box, itération " << i <<endl;
 			box.vmin.x = std::min(box.vmin.x, vertices[indices[i].vtxi].x);
 			box.vmin.y = std::min(box.vmin.y, vertices[indices[i].vtxi].y);
 			box.vmin.z = std::min(box.vmin.z, vertices[indices[i].vtxi].z);
 			box.vmax.x = std::max(box.vmax.x, vertices[indices[i].vtxi].x);  //itérer pour i, j, k
 			box.vmax.y = std::max(box.vmax.y, vertices[indices[i].vtxi].y);
 			box.vmax.z = std::max(box.vmax.z, vertices[indices[i].vtxi].z);

 			box.vmin.x = std::min(box.vmin.x, vertices[indices[i].vtxj].x);
 			box.vmin.y = std::min(box.vmin.y, vertices[indices[i].vtxj].y);
 			box.vmin.z = std::min(box.vmin.z, vertices[indices[i].vtxj].z);
 			box.vmax.x = std::max(box.vmax.x, vertices[indices[i].vtxj].x); //itérer pour i, j, k
 			box.vmax.y = std::max(box.vmax.y, vertices[indices[i].vtxj].y);
 			box.vmax.z = std::max(box.vmax.z, vertices[indices[i].vtxj].z);

 			box.vmin.x = std::min(box.vmin.x, vertices[indices[i].vtxk].x);
 			box.vmin.y = std::min(box.vmin.y, vertices[indices[i].vtxk].y);
 			box.vmin.z = std::min(box.vmin.z, vertices[indices[i].vtxk].z);
 			box.vmax.x = std::max(box.vmax.x, vertices[indices[i].vtxk].x); //itérer pour i, j, k
 			box.vmax.y = std::max(box.vmax.y, vertices[indices[i].vtxk].y);
 			box.vmax.z = std::max(box.vmax.z, vertices[indices[i].vtxk].z);
 		}
		cout << "build_box a tourné. \n";
 		return box;
 	 }

	 double splitletter(Vector v, int splitdim) {
			if (splitdim == 0) {
				return v.x;
			} else if (splitdim == 1) {
				return v.y;
			} else { return v.z; }
		}

	 void build_bhv(BHVNode* node, const int &debut, const int &fin) {
		 if (fin-debut > 50) {
			 cout << "starting to build" << endl;
			 cout << "début : " << debut << " fin : " << fin << endl;
			 // mtex.lock();
			 // node->setlimites(debut, fin);
				node->debut = debut;
				node->fin = fin;
				Bbox b = build_box(debut, fin);
				cout << "first bbox built" << endl;
				node->b = b;
				Vector diagBbox = b.vmax - b.vmin;
				int splitdim;
				if (diagBbox.x > diagBbox.y && diagBbox.x > diagBbox.z ) {
					splitdim = 0;
				} else {
					if (diagBbox.y > diagBbox.x && diagBbox.y > diagBbox.z ) {
						splitdim = 1;
					} else {splitdim	= 2;}
				}
				Vector centerBbox = b.vmin + diagBbox*0.5;
				int p = debut;
				for (int i = debut; i < fin; i++) {
					double centerTri = (splitletter(vertices[indices[i].vtxi], splitdim) + splitletter(vertices[indices[i].vtxj], splitdim) + splitletter(vertices[indices[i].vtxk], splitdim))/3;
					if (centerTri < splitletter(centerBbox, splitdim)){
						std::swap(indices[i], indices[p]);
						p++;
					}
				}
				cout << "first loop ended" << endl;
				if (fin-debut>3 && debut<p && p<fin){
					Bbox b1 = build_box(debut, p);
					Bbox b2 = build_box(p+1, fin);
					node->fg = new BHVNode(b1, debut, p);
					node->fd = new BHVNode(b2, p+1, fin);
					cout << "entering second loop" << endl;
					// mtex.unlock();
					build_bhv(node->fg, debut, p);
					build_bhv(node->fd, p+1, fin);
				}
				cout << "built left and right boxes" << endl;
		 }

		}

	 void add_textures(const char* filename) {
	 	int w, h, c;
	 	textures.push_back(stbi_load(filename, &w, &h, &c, 3));
	 	textures_widths.push_back(w);
	 	textures_heights.push_back(h);
	 }

	 void build_matrix(const Vector &translation, const double &rz /*la rotation en radians */, const double &scaling) {
		 matrix.resize(3*4);
		 matrix[0] = sin(rz)*scaling;
		 matrix[1] = cos(rz)*scaling;
		 matrix[2] = 0;
		 matrix[3] = translation.x;

		 matrix[4] = -cos(rz)*scaling;
		 matrix[5] = sin(rz)*scaling;
		 matrix[6] = 0;
		 matrix[7] = translation.y;

		 matrix[8] = 0;
		 matrix[9] = 0;
		 matrix[10] = 0;
		 matrix[11] = translation.z;

		 inversematrix.resize(3*4);
		 inversematrix[0] = sin(rz)/scaling;
		 inversematrix[1] = -cos(rz)/scaling;
		 inversematrix[2] = 0;
		 inversematrix[3] = -(sin(rz)*translation.x - cos(rz)*translation.y) / scaling;

		 inversematrix[4] = cos(rz)/scaling;
		 inversematrix[5] = sin(rz)*scaling;
		 inversematrix[6] = 0;
		 inversematrix[7] = -(cos(rz)*translation.x + sin(rz)*translation.y) / scaling;

		 inversematrix[8] = 0;
		 inversematrix[9] = 0;
		 inversematrix[10] = 1/scaling;
		 inversematrix[11] = -1/scaling*translation.z;


		 inversetransposematrix.resize(4*3);
		 inversetransposematrix[0] = sin(rz)/scaling;
		 inversetransposematrix[1] = cos(rz)/scaling;
		 inversetransposematrix[2] = 0;

		 inversetransposematrix[3] = -cos(rz)/scaling;
		 inversetransposematrix[4] = sin(rz)*scaling;
		 inversetransposematrix[5] = 0;

		 inversetransposematrix[6] = 0;
		 inversetransposematrix[7] = 0;
		 inversetransposematrix[8] = 1/scaling;

		 inversetransposematrix[9] = -(sin(rz)*translation.x - cos(rz)*translation.y) / scaling;
		 inversetransposematrix[10] = -(cos(rz)*translation.x + sin(rz)*translation.y) / scaling;
		 inversetransposematrix[11] = -1/scaling*translation.z;


	 }

	 std::vector<unsigned char*> textures;
	 std::vector<int> textures_widths;
	 std::vector<int> textures_heights;
	 std::vector<double> matrix, inversematrix, inversetransposematrix;

	 void intersection_bhv(const Ray &r, BHVNode* &bhv, vector<int> &indicesTrianglesIntersectes) {
		 if (bhv->fg != NULL || bhv->fd != NULL) {
			 if(bhv->fg->b.inter(r)){
				 intersection_bhv(r, bhv->fg, indicesTrianglesIntersectes);
			 }
			 if (bhv->fd->b.inter(r)){
				 intersection_bhv(r, bhv->fd, indicesTrianglesIntersectes);
			 } else {
				 for (int i=bhv->debut; i<=bhv->fin; i++) {
					 indicesTrianglesIntersectes.push_back(i);
				 }
			 }
		 }
	 }

	 Ray ray_transform(vector<double> &matrix, const Ray r) {
 		// la matrice doit être de taille 3*4
 		Ray rt = r;
 		rt.u.x = matrix[0]*r.u.x + matrix[1]*r.u.y + matrix[2]*r.u.z + matrix[3];
		rt.u.y = matrix[4]*r.u.x + matrix[5]*r.u.y + matrix[6]*r.u.z + matrix[7];
		rt.u.z = matrix[8]*r.u.x + matrix[9]*r.u.y + matrix[10]*r.u.z + matrix[11];
		rt.C.x = matrix[0]*r.C.x + matrix[1]*r.C.y + matrix[2]*r.C.z + matrix[3];
		rt.C.y = matrix[4]*r.C.x + matrix[5]*r.C.y + matrix[6]*r.C.z + matrix[7];
		rt.C.z = matrix[8]*r.C.x + matrix[9]*r.C.y + matrix[10]*r.C.z + matrix[11];
		return rt;
 	}

	 virtual bool intersect(const Ray& r, Vector& P, Vector& N) {
		 	double smallestt = 1E15;
			bool has_inter = false;

			if (BHV.b.inter(r)) {
				// cout << "intersection avec l'enveloppe de la forme géométrique" << endl;
				BHVNode* boiteloc = new BHVNode; // new BHVNode éventuellement à supprimer.
				boiteloc = &BHV;
				// cout << " a attribué boiteloc " << endl;
				vector<int> indicesTrianglesIntersectes;
				Ray r_trans = ray_transform(inversematrix, r);
				intersection_bhv(r, boiteloc, indicesTrianglesIntersectes);

				// while (boiteloc->fg != NULL || boiteloc->fd != NULL) {
				// 	if (boiteloc->fg != NULL && boiteloc->fg->b.inter(r)) {
				// 		boiteloc = boiteloc->fg;
				// 	} else if (boiteloc->fd != NULL && boiteloc->fd->b.inter(r)) {
				// 		boiteloc = boiteloc->fd;
				// 	} else {
				// 		for (int i=boiteloc->debut; i<=boiteloc->fin; i++) {
				// 			indicesTrianglesIntersectes.push_back(i);
				// 		}
				// 		break;
				// 	}
				// }

				// cout << "nombre de triangles intersectés : " << indicesTrianglesIntersectes.size() << endl;
				// while (boiteloc->debut != boiteloc->fin) {
				// 	// cout << "entrée dans la première boucle " << endl;
				// 	if (boiteloc->fg != NULL) {
				// 		// cout << " fg n'est pas NULL" <<endl;
				// 		if (boiteloc->fg->b.inter(r)) {
				// 			// cout << "c'est dans la boite de gauche" <<endl;
				// 			boiteloc = boiteloc->fg;
				// 		} else {
				// 			// cout << "on part à droite" <<endl;
				// 			boiteloc = boiteloc->fd;
				// 		}
				// 	} else {
				// 		if (boiteloc->fd != NULL) {boiteloc = boiteloc->fd;}
				// 		else {
				// 			break;
				// 		}
				// 	}
				// }
				// cout << "début : " << boiteloc->debut << " fin : " << boiteloc->fin << endl;
				for (int k=0; k<indicesTrianglesIntersectes.size(); k++) {
					int i = indicesTrianglesIntersectes[k];
					// int i = boiteloc->debut;
					// cout << i <<endl;
					Triangle triloc = Triangle({vertices[indices[i].vtxi], vertices[indices[i].vtxk], vertices[indices[i].vtxj]});
					Vector Ploc, Nloc;
					double tloc, betaloc, gammaloc;
					bool inter = triloc.intersect(r, Ploc, Nloc, betaloc, gammaloc, tloc);
					// cout<< "test de triangle" <<endl;
					// if (inter) {cout << "intersection avec le triangle" <<endl;}
					// cout << "tloc : " << tloc <<endl;
					if (inter && tloc < smallestt) {
						smallestt = tloc;
						P = Ploc;
						N = (1-betaloc-gammaloc)*normals[indices[i].ni] + betaloc*normals[indices[i].nj] + gammaloc*normals[indices[i].nk];
						N.normalize();
						if (indices[i].group < 0) {
							this->albedo;
						} else {
							int tW = textures_widths[indices[i].group];
							int tH = textures_heights[indices[i].group];
							Vector UV = (1-betaloc-gammaloc)*uvs[indices[i].uvi] + betaloc*uvs[indices[i].uvj] + gammaloc*uvs[indices[i].uvi];
							int x = fabs(fmod(UV.x, 1.)*tW);
							int y = tH - fabs(fmod(UV.y, 1.)*tH);
							albedo.x = textures[indices[i].group][(y*tW+x)*3 + 0] / 255.;
							albedo.y = textures[indices[i].group][(y*tW+x)*3 + 1] / 255.;
							albedo.z = textures[indices[i].group][(y*tW+x)*3 + 2] / 255.;
						}
						has_inter = true;
						// cout << "le triangle est plus proche que les précédents" << endl;
					}
				}

			}

			return has_inter;
	 }

	 Bbox enveloppe_ext;
	 BHVNode BHV;
	 std::vector<TriangleIndices> indices;
	 std::vector<Vector> vertices;
	 std::vector<Vector> normals;
	 std::vector<Vector> uvs;
	 std::vector<Vector> vertexcolors;
};


class Sphere : public Objet {
public:
	Sphere(const Vector& O = Vector(0, 0, 0), double R = 0, const Vector& albedo = Vector(0, 0, 0), const bool& miroir = false, const bool &transparent = false, const double& n = 1) : O(O), R(R) {
		this->albedo = albedo;
		this-> miroir = miroir;
		this->transparent = transparent;
		this->n = n;
	};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N) {
		// r�solution de t�|u|� + 2t*(u scalaire C-O) + |C-O|� = R�
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - O).norm2() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta >= 0) {
			double t1 = (-b - sqrt(delta)) / (2 * a);
			double t2 = (-b + sqrt(delta)) / (2 * a);
			if (t1 > 0) {
				P = t1 * r.u + r.C;
				N = P - O;
				N.normalize();
			}
			else {
				if (t2 > 0) {
					P = t2 * r.u + r.C;
					N = P - O;
					N.normalize();
				}
				else {
					return false;
				}
			}
		}

		return delta >= 0;
	};
	Vector O;
	double R;
};

class Light {
public:
	Light(const Vector& pos, const double& I, const double & Rayon) : pos(pos), I(I), R(Rayon) {};
	Vector pos;
	double I, R;
};

class Scene {
public:
	Scene(const Light& L) : L(L) {};
	void addObject(Objet* o) {
		objets.push_back(o);
	}
	bool intersect(const Ray& r, Vector& P, Vector& N, Objet*& S) {
		double dist_min = 1E15;
		bool found = false;
		for (int i = 0; i < objets.size(); i++) {
			Objet* Sloc = objets[i];
			Vector Ploc, Nloc;
			if (Sloc->intersect(r, Ploc, Nloc)) {
				double dist_cam = (Ploc - r.C).norm2();
				if (dist_cam < dist_min) {
					dist_min = dist_cam;
					P = Ploc;
					N = Nloc;
					S = Sloc;
				}
				found = true;
			}
		}
		return found;
	}

	Vector cos_aleatoire(const Vector& N) {
		double r1 = uniforme(seed);
		double r2 = uniforme(seed);
		double s = sqrt(1-r2);
		Vector T1;
		if (abs(N.x)<=abs(N.y) && abs(N.x)<=abs(N.z)) {
			T1 = Vector(0, -N.z, N.y);
		} else {
			if (abs(N.y)<=abs(N.x) && abs(N.y)<=abs(N.z)) {
				T1 = Vector(-N.z, 0, N.x);
			} else {
				T1 = Vector(-N.y, N.x, 0);
			}
		}
		T1.normalize();
		Vector T2 = pd_vect(N, T1);
		Vector w1 = cos(2*M_PI*r1)*s*T1;
		Vector w2 = sin(2*M_PI*r1)*s*T2;
		Vector w3 = N*sqrt(r2);
		Vector w = w1+w2+w3;
		return w;
	}

	Vector getColor(Ray& r, Light& L, const int &nb_rebonds, const double &n_air) {
		Vector pixelColor(0, 0, 0);
		Vector P, N;
		Objet* S;
		bool epsilon = 0.05;
		if (intersect(r, P, N, S)) {
									if (S->miroir && nb_rebonds > 0) {
										r.u = r.u - 2*dot(r.u, N)*N;
										r.C = P + epsilon*N;
										return getColor(r, L, nb_rebonds - 1, n_air);
									}
									if (S->transparent  && nb_rebonds>0){
										if (dot(r.u, N)<0){
											r.u = n_air/S->n * (r.u - dot(r.u, N) * N) - sqrt(1 - pow(n_air/S->n,2)*(1 - pow(dot(r.u, N), 2))) * N;
											r.C = P - epsilon*N;
										} else {
											r.u = S->n/n_air * (r.u - dot(r.u, -N) * (-N)) + sqrt(1 - pow(S->n/n_air,2)*(1 - pow(dot(r.u, -N), 2))) * N;
											r.C = P + epsilon*N;
										}
										return getColor(r, L, nb_rebonds - 1, n_air);
									}

									else {
										// éclairage direct
										Vector PL = L.pos - P;
										double carre_dist_lum = PL.norm2();
										PL.normalize();

										Vector Nobjet = cos_aleatoire(-PL);
										Nobjet.normalize();
										Vector xprime = L.pos + L.R*Nobjet;
										Vector xprimex = P - xprime;
										Vector omegai = -xprimex;
										omegai.normalize();
										Ray rayon_lumiere(P + epsilon*omegai, omegai);

										Vector Pprime, Nprime;
										Objet* Sprime;
										if (intersect(rayon_lumiere, Pprime, Nprime, Sprime) && (P - Pprime).norm2() < xprimex.norm2()) {
										} else {
											pixelColor = L.I/(4*pow(M_PI,2)*pow(L.R,2)) * (S->albedo / M_PI) * dot(omegai, PL)*dot(-omegai,Nobjet) / (xprimex.norm2()*dot(-PL,Nobjet)/(M_PI*pow(L.R, 2)));
										}

										// Calcul de l'ombre
										// Pour eviter le bruit dans l'image on d�cale l�g�rement le point d'origine
										// Ray rayon_lumiere(P + epsilon * PL, PL);
										// Vector Pprime, Nprime;
										// Sphere Sprime;
										// // // // On a une ombre si le rayon qui part du point d'intersection intersecte un objet
										// // // // et si cet objet est plus proche que la lumi�re
										// bool ombre = intersect(rayon_lumiere, Pprime, Nprime, Sprime);
										// if (ombre && (P - Pprime).norm2() < carre_dist_lum && Sprime.transparent == false) {
										// 	// pixelColor = Vector(0, 0, 0);
										// 	}
										// else {
										// 		pixelColor = L.I * (S.albedo / M_PI) * std::max(0., dot(N, PL)) / (4 * M_PI*carre_dist_lum);
										// 	}


										// éclairage indirect
										if (nb_rebonds>0) {
											Vector reflechi = cos_aleatoire(N);
											Ray ray_reflechi(P + 0.01*N, reflechi);
											pixelColor = pixelColor + S->albedo*getColor(ray_reflechi, L, nb_rebonds - 1, n_air);
										}


									}
			}

		else {
			pixelColor = Vector(0, 0, 0);
		}
		return pixelColor;
	}
	std::vector<Objet*> objets;
	Light L;
};

int main() {

	// std::default_random_engine e;
	// uniform_real_distribution<double>  U(0,1);
	// double r1, r2, theta1, theta2, theta3, integrale, p, ecart_type, f;
	// ecart_type = 0.2;
	// integrale = 0;
	// int N,i;
	// N = 1000000;
	// for (i=0; i<N; i++) {
	// 	// r1 = U(e);
	// 	// r2 = U(e);
	// 	// theta1 = ecart_type*cos(2*M_PI*r1)*sqrt(-2*log(r2));
	// 	// theta2 = ecart_type*sin(2*M_PI*r1)*sqrt(-2*log(r2));
	// 	// r1 = U(e);
	// 	// r2 = U(e);
	// 	theta1 = M_PI*U(e)-M_PI/2;
	// 	theta2 = M_PI*U(e)-M_PI/2;
	// 	theta3 = M_PI*U(e)-M_PI/2;
	// 	// theta3 = ecart_type*cos(2*M_PI*r1)*sqrt(-2*log(r2));
	// 	// p = 1/(pow(ecart_type*sqrt(2*M_PI),3))*exp(-(pow(theta1,2)+pow(theta2,2)+pow(theta3,2))/(2*pow(ecart_type,2)));
	// 	p = 1/pow(M_PI,3);
	// 	f = pow(cos(theta1*theta2),30)*pow(cos(theta1*theta3),30);
	// 	integrale = integrale + f/(p*N);
	// }
	// cout << integrale << endl;



	int W = 100; // 512;
	int H = 100; // 512;
	int NRays = 2;
	// Camera
	double fov = 60 * M_PI / 180.;
	std::vector<unsigned char> image(W*H * 3, 0);

	// Light
	Light L(Vector(-10, 20, 40), 5000000000, 5);

	// Scene
	// Triangle t1({Vector(10,5,10), Vector(-10,2,-10), Vector(5,-10, 0)}, Vector(1., 1., 1.), false, false, 1);
	// Sphere s1(Vector(0, 0, 0), 5, Vector(1., 0.5, 0.), false, true, 1.6);
	// Sphere s2(Vector(15, 15, -15), 10, Vector(0.5, 1., 1.), true, false, 0.5);
	// Sphere s3(Vector(-40, 0, -50), 40, Vector(0., 0.7, 1.), false, false, 1);
	// Sphere s4(Vector(-10, 0, 0), 3, Vector(0, 0.1, 1), false, false, 1.6);
	// Sphere s5(Vector(0, -5, 0), 3, Vector(1, 0.3, 0), false, false, 1.6);
	// Sphere s6(Vector(10, -5, 3), 5, Vector(1, 0.5, 1), false, false, 1.6);
	Sphere ssol(Vector(0, -1000, 0), 1000-30, Vector(0.7, 0., 0.4), false, false, 1);
	Sphere sfond(Vector(0, 0, -1000), 1000-60, Vector(1, 1, 1), true, false, 1);
	Sphere sarriere(Vector(0, 0, 1000), 1000-60, Vector(0, 0.1, 1), false, false, 1);
	Sphere splafond(Vector(0, 1000, 0), 1000-60, Vector(1, 1, 1), false, false, 1);
	Sphere sdroite(Vector(1000, 0, 0), 1000-60, Vector(1, 0.4, 1), false, false, 1);
	Sphere sgauche(Vector(-1000, 0, 0), 1000-60, Vector(1, 0.5, 0.8), false, false, 1);
	const char* filename = "BeautifulGirl_original.obj";
	Geometry beautifulGirl(filename, 20., Vector(0,-1,0), Vector(1,1,1), false, false, 1);
	beautifulGirl.add_textures("Girltextures/visage.bmp");
	beautifulGirl.add_textures("Girltextures/cheveux.bmp");
	beautifulGirl.add_textures("Girltextures/corps.bmp");
	beautifulGirl.add_textures("Girltextures/pantalon.bmp");
	beautifulGirl.add_textures("Girltextures/accessoires.bmp");
	beautifulGirl.add_textures("Girltextures/mains.bmp");
	beautifulGirl.build_matrix(Vector(10.,0.,10.), 1.5, 1.2); // translation, rotation d'axe z en radians, mise à l'échelle.
	Scene scene(L);
	// scene.addObject(&t1);
	// scene.addObject(&s1);
	// scene.addObject(&s2);
	// scene.addObject(&s3);
	// scene.addObject(&s4);
	// scene.addObject(&s5);
	// scene.addObject(&s6);
	scene.addObject(&ssol);
	scene.addObject(&sfond);
	scene.addObject(&sarriere);
	scene.addObject(&splafond);
	scene.addObject(&sgauche);
	scene.addObject(&sdroite);
	scene.addObject(&beautifulGirl);


	double epsilon = 0.01;
	double n_air = 1;
	double nb_rebonds = 3;
	double pos_objet_focal = 0; // coordonnée z de l'objet focal net.

	cout << "Creating the pixels" << endl;
	// #pragma omp parallel for
	for (int i = 0; i < W; i++) {
		if (i % 10 == 0) {std::cout << "ligne : " << i << std::endl;}
		for (int j = 0; j < H; j++) {
			// std::cout << "colonne : " << j << std::endl;
			Vector pixelColor;
			// pixelColor = scene.getColor(r, L, nb_rebonds, n_air);
			for (int rayid=0; rayid<NRays; rayid++) {
				double r1 = uniforme(seed);
				double r2 = uniforme(seed);
				double offsetx = cos(2*M_PI*r1)*sqrt(-2*log(r2))*0.5;  // on multiplie par l'écart type 0.25 de pixel
				double offsety = sin(2*M_PI*r1)*sqrt(-2*log(r2))*0.5;
				Vector u = Vector(j-W/2 + offsetx, -i+H/2 + offsety, -W/(2*tan(fov/2)));
				u.normalize();
				r1 = uniforme(seed);
				r2 = uniforme(seed);
				offsetx = cos(2*M_PI*r1)*sqrt(-2*log(r2))*0.5;  // on multiplie par l'écart type
				offsety = sin(2*M_PI*r1)*sqrt(-2*log(r2))*0.5;
				Vector C(0, 0, 55);
				C.x = C.x + offsetx;
				C.y = C.y + offsety;
				u.x = u.x + atan(offsetx/(pos_objet_focal - C.z));
				u.y = u.y + atan(offsety/(pos_objet_focal - C.z));
				// u = u + Vector(atan(offsetx/(pos_objet_focal - C.z)), atan(offsety/(pos_objet_focal - C.z)), 0);
				u.normalize();
				Ray r = Ray(C, u);
				pixelColor = pixelColor + scene.getColor(r, L, nb_rebonds, n_air)/NRays;
			}
			image[(i*W + j) * 3 + 0] = std::min(255., pow(pixelColor.x, 0.45));
			image[(i*W + j) * 3 + 1] = std::min(255., pow(pixelColor.y, 0.45));
			image[(i*W + j) * 3 + 2] = std::min(255., pow(pixelColor.z, 0.45));
		}
	}
	cout << "Writing image" << endl;
	stbi_write_png("image2.png", W, H, 3, &image[0], 0);
	cout << "Finished writing image" << endl;
	// delete []; // ceci est un test
	return 0;
}
