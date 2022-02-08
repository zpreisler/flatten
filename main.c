#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define SQR(x) ((x)*(x))
#define NEDGES 8

typedef struct edge{
	double length;
	void **vertices;
}edge;

typedef struct vec{
	double x,y,z;
}vec;

typedef struct vertex{
	struct vertex **connected;
	vec *pos; //pointer to the value; position
	vec *force; //pointer to the value; force
	int nedges;
	edge **edges;
}vertex;

double dist(vec *v1,vec *v2){
	double d;
	d = sqrt(SQR(v1->x-v2->x) + SQR(v1->y-v2->y) + SQR(v1->z-v2->z));
	return d;
}

int mesh(vec *pos,int nx,int ny){

	int i,j;
	vec *p;

	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			 p = (pos+i+j*nx);
			 p->x = i;
			 p->y = j;
			 p->z = 0.5 * (SQR(i-(nx-1)*0.5) + SQR(j-(ny-1)*0.5));
		}
	}
	return 0;
}

int emesh(vertex *v,edge *edges,int nx,int ny){
	int i,j,k=0;
	vertex *v1,*v2;
	edge *e;

	for(i=0;i<nx-1;i++){
		for(j=0;j<ny;j++){
			v1 = v + i + j * nx;
			v2 = v + i + 1 + j * nx;
			//printf("%p %p %lf\n",v1,v2,(v2->pos->x) - (v1->pos->x));
			e = edges + k++;
			e->length = dist(v1->pos,v2->pos);

			*(e->vertices) = v1;
			*(e->vertices+1) = v2;

			*(v1->connected + v1->nedges) = v2;
			*(v2->connected + v2->nedges) = v1;

			*(v1->edges+v1->nedges++) = e;
			*(v2->edges+v2->nedges++) = e;
		}
	}

	for(i=0;i<nx;i++){
		for(j=0;j<ny-1;j++){
			v1 = v + i + j * nx;
			v2 = v + i + (j+1) * nx;
			//printf("%p %p %lf\n",v1,v2,(v2->pos->x) - (v1->pos->x));
			e = edges + k++;
			e->length = dist(v1->pos,v2->pos);

			*(e->vertices) = v1;
			*(e->vertices+1) = v2;

			*(v1->connected + v1->nedges) = v2;
			*(v2->connected + v2->nedges) = v1;

			*(v1->edges+v1->nedges++) = e;
			*(v2->edges+v2->nedges++) = e;
		}
	}

	for(i=0;i<nx-1;i++){
		for(j=0;j<ny-1;j++){
			v1 = v + i + 1 + j * nx;
			v2 = v + i + (j+1) * nx;
			//printf("%p %p %lf\n",v1,v2,(v2->pos->x) - (v1->pos->x));
			e = edges + k++;
			e->length = dist(v1->pos,v2->pos);

			*(e->vertices) = v1;
			*(e->vertices+1) = v2;

			*(v1->connected + v1->nedges) = v2;
			*(v2->connected + v2->nedges) = v1;

			*(v1->edges+v1->nedges++) = e;
			*(v2->edges+v2->nedges++) = e;
		}
	}

	return 0;
}

int dump_vec(vec *v){
	printf("%f %f %f\n",v->x,v->y,v->z);
	return 0;
}
int dump_edges(vertex *v){
	int i;
	edge *e;
	vertex *v1,*v2;
	for(i=0;i<v->nedges;i++){
		e = *(v->edges + i);
		v1 = *e->vertices;
		v2 = *(e->vertices+1);
		//printf("v:%p e:%p v1:%p v2:%p\n",v,e,v1,v2);
		printf("# %lf dist:%lf\n",e->length,dist(v1->pos,v2->pos));
		dump_vec(v1->pos);
		dump_vec(v2->pos);
		printf("\n\n");
	}
}

int interactions(vertex *vertices,int nvertices){
	int i,j;
	vertex *v,*c;
	vertex *v1,*v2;
	edge *e;
	double d,d2,dd,dh;
	double dx,dy,dz;
	double rep;

	for(i=0;i<nvertices;i++){
		v = vertices + i;

		//v->force->z += -0.5 * v->pos->z; //Gravity

		for(j=0;j<v->nedges;j++){
			e = *(v->edges + j);

			c = *(v->connected + j);
			//printf("%d %lf %lf\n",j,v->pos->x,c->pos->x);

			dx = v->pos->x - c->pos->x;
			dy = v->pos->y - c->pos->y;
			dz = v->pos->z - c->pos->z;

			d2 = SQR(dx) + SQR(dy) + SQR(dz);
			d = sqrt(d2);

			dd = (d - e->length); //Harmonic spring
			
			dh = -dd/(sqrt(dd*dd+1));
			//printf("rep: %lf\n",rep);
			
			//dx = dx/d * dh + dx/d * rep; 
			//dy = dy/d * dh + dy/d * rep;

			dx = dx/d * dh; 
			dy = dy/d * dh;
			dz = dz/d * dh;

			v->force->x += dx;
			v->force->y += dy;
			v->force->z += dz;
		}
	}

	for(i=0;i<nvertices-1;i++){
		v1 = vertices + i;
		for(j=i+1;j<nvertices;j++){
			v2 = vertices + j;

			dx = v1->pos->x - v2->pos->x;
			dy = v1->pos->y - v2->pos->y;
			dz = v1->pos->z - v2->pos->z;

			d2 = SQR(dx) + SQR(dy) + SQR(dz);
			d = sqrt(d2);

			rep = 1/(dx*dx + dy*dy + dz*dz) * 0.5 * 0.1;
			//printf("rep: %lf d:%lf\n",rep,d);

			dx = dx/d * rep; 
			dy = dy/d * rep;

			v1->force->x += dx;
			v1->force->y += dy;

			v2->force->x -= dx;
			v2->force->y -= dy;
		}
	}

	return 0;
}

void reset_forces(vec *forces,int nvertices){
	int i;
	for(i=0;i<nvertices;i++){
		(forces+i)->x = 0.0;
		(forces+i)->y = 0.0;
		(forces+i)->z = 0.0;
	}
}

void move(vertex *vertices,int nvertices,double dt){
	int i;
	vertex *v;

	for(i=0;i<nvertices;i++){
		v = vertices + i;
		v->pos->x += v->force->x * dt;
		v->pos->y += v->force->y * dt;
		v->pos->z += v->force->z * dt;
	}
		
}

void dump_forces(vertex *vertices,int nvertices){
	int i;
	vertex *v;
	for(i=0; i<nvertices; i++){
		v = (vertices + i);
		printf("%lf %lf %lf\n",v->force->x,v->force->y,v->force->z);
	}
}

int main(int argc,char *argv[]){

	int i;
	int nx = 6;
	int ny = 6;
	int nvertices = nx * ny;
	vertex *v;
	vertex *vertices = (vertex*)malloc(sizeof(vertex) * nvertices);
	vec *pos = (vec*)malloc(sizeof(vec) * nvertices);
	vec *forces = (vec*)malloc(sizeof(vec) * nvertices);
	edge *edges = (edge*)malloc(sizeof(edge) * nvertices * nvertices);
	for(i=0;i<nvertices*nvertices;i++){
		(edges + i)->vertices = (void**)malloc(sizeof(vertex*) * 2);
	}

	for(i=0; i<nvertices; i++){
		v = (vertices + i);
		v->edges = (edge**)malloc(sizeof(edge*) * NEDGES);
		v->connected = (vertex**)malloc(sizeof(vertex*) * NEDGES);
		v->nedges = 0;
		v->pos = pos + i;
		v->force = forces + i;
	}

	mesh(pos,nx,ny);
	emesh(vertices,edges,nx,ny);

	//for(i=0; i<nvertices; i++){
		//dump_vec(pos+i);
	//	dump_edges(vertices+i);
	//}
	//
	for(i=0; i<nvertices; i++){
		v = vertices + i;
		v->pos->z = 0;
	}

	
	for(i=0;i<10000;i++){
		interactions(vertices,nvertices);
		move(vertices,nvertices,0.001);
		reset_forces(forces,nvertices);
	}

	for(i=0; i<nvertices; i++){
		//dump_vec(pos+i);
		dump_edges(vertices+i);
	}

	//dump_forces(vertices,nvertices);

	return 0;
}

