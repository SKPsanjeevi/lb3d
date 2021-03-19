
/*
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <unistd.h>
#include <math.h>

#include <netinet/in.h>
#include "vol3d.h"

int load_voltex(char *fname, int nx, int ny, int nz) ;

int window; /* GLUT window ID */
int leftbutstate=0;
int middlebutstate=0;
int rightbutstate=0;


struct quat {
	GLfloat q0,q1,q2,q3; /* Quaternion structure to hold orientation */
};

typedef struct quat quat_t;



/*********** globals for main code ************/

const GLfloat nearclip=0.1; /* Distance to near clipping plane */
int screenx,screeny;	/* Screen dimensions */
int clickx=0,clicky=0;	/* Coords of last mouse click */

GLfloat sphererad=0.0; /* Radius of trackball sphere */

int glinit=0; /* Set to 1 once GL has been set up */
unsigned int frameno=0;

/*** Things which are global but shouldn't be ***/

/* Transfer function lookup tables */
uint8_t r[256],g[256],b[256],a[256];


/* Quaternion for global orientation */
quat_t orquat;

GLfloat geomx=0,geomy=0,geomz=-5;

/************** Globals for mouse state **************/

int mousex=0,mousey=0;		/* Mouse status */
int leftbutdown=0;
int rightbutdown=0;
int midbutdown=0;
int bx=0,by=0;


/*********** Quaternion stuff ****************/

/* Turn a quaternion into a matrix */

void quatmatrix(quat_t q,GLfloat *m )
{
	GLfloat w,x,y,z;

	w = q.q0; x = q.q1; y=q.q2; z=q.q3;

	m[ 0] = 1 - 2*y*y - 2*z*z;
	m[ 1] = 2*x*y + 2*w*z;
	m[ 2] = 2*x*z - 2*w*y;
	m[ 3] = 0;

	m[ 4] = 2*x*y - 2*w*z;
	m[ 5] = 1 - 2*x*x - 2*z*z;
	m[ 6] = 2*y*z + 2*w*x;
	m[ 7] = 0;

	m[ 8] = 2*x*z + 2*w*y;
	m[ 9] = 2*y*z - 2*w*x;
	m[10] = 1 - 2*x*x - 2*y*y;
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;

	return;

}

/* Turn a quaternion into a matrix and its inverse (ie transpose) */

void quatmatrixandinv(quat_t q,GLfloat *m,GLfloat *mi )
{
	GLfloat w,x,y,z;

	w = q.q0; x = q.q1; y=q.q2; z=q.q3;

	mi[ 0] = m[ 0] = 1 - 2*y*y - 2*z*z;
	mi[ 4] = m[ 1] = 2*x*y + 2*w*z;
	mi[ 8] = m[ 2] = 2*x*z - 2*w*y;
	mi[12] = m[ 3] = 0;

	mi[ 1] = m[ 4] = 2*x*y - 2*w*z;
	mi[ 5] = m[ 5] = 1 - 2*x*x - 2*z*z;
	mi[ 9] = m[ 6] = 2*y*z + 2*w*x;
	mi[13] = m[ 7] = 0;

	mi[ 2] = m[ 8] = 2*x*z + 2*w*y;
	mi[ 6] = m[ 9] = 2*y*z - 2*w*x;
	mi[10] = m[10] = 1 - 2*x*x - 2*y*y;
	mi[14] = m[11] = 0;

	mi[ 3] = m[12] = 0;
	mi[ 7] = m[13] = 0;
	mi[11] = m[14] = 0;
	mi[15] = m[15] = 1;

	return;

}

/* Construct the quaternion corresponding to a right-handed rotation of
 * theta radians about the given vector.
 */

quat_t quatrotation(GLfloat theta, GLfloat x, GLfloat y, GLfloat z)
{
	GLfloat st,ct;
	quat_t q;

	st = sin(0.5*theta);
	ct = cos(0.5*theta);

	q.q0 = ct;
	q.q1 = x*st;  
	q.q2 = y*st;  
	q.q3 = z*st;  
	return q;
}
/*
 * Return the quaternion product of two quaternions.
 *
 * ( a0, a ) ( b0, b) =  ( a0*b0 - a.b , a0*a + b0*b + a x b )
 *
 */

quat_t quatmultiply(quat_t a, quat_t b)
{
	quat_t q;

	q.q0 = a.q0 * b.q0 
		- a.q1*b.q1
		- a.q2*b.q2
		- a.q3*b.q3;
	
	q.q1 =    a.q0 * b.q1
		+ b.q0 * a.q1
		+ a.q2 * b.q3 - a.q3 * b.q2;

	q.q2 =    a.q0 * b.q2
		+ b.q0 * a.q2
		+ a.q3 * b.q1 - a.q1 * b.q3;

	q.q3 =    a.q0 * b.q3
		+ b.q0 * a.q3
		+ a.q1 * b.q2 - a.q2 * b.q1;
	return q;
}



/********************* Main code *********************/

void init_gl()
{

	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_EMISSION);
	glColor4f(1,0,0,1);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glColor4f(1,1,1,1);
	
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	
	load_voltex("32x32x32.dat",32,32,32);

	glEnable(GL_TEXTURE_3D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.2,0.2,0.2,0);


	glinit=1;
}







/* Process mouse motion from x0,y0 to x1,y1 */

void mouserotatemotion(int x0,int y0, int x1,int y1) {

	GLfloat mx,my,m,s;
	GLfloat theta;
	quat_t rotquat;

	s = sphererad;
	my = (GLfloat)(x1-x0);
	mx = (GLfloat)(y1-y0);
	m=sqrt(mx*mx+my*my);

	if ((m>0) && (m<s)) {
		theta = m/s;

		mx /= m;
		my /= m;

		rotquat = quatrotation(theta,mx,my,0.0);
		orquat = quatmultiply(rotquat,orquat);
	}

}

void mousezoommotion(int dz)
{
	const GLfloat zoomscale=0.1;
	geomz-=zoomscale*dz;
	return;
}

void mousetransmotion(int dx,int dy)
{
	const GLfloat transscale=0.06;
	geomx-=transscale*dx;
	geomy+=transscale*dy; /* GL uses left-handed coordinates */
	return;
}


/* Clear up and quit. */
void quit_all() {

	exit(0);
}

void reshapefunc(int width, int height)
{
	printf("Reshaped to %d x %d\n",width,height);
	screenx = width;
	screeny = height;
	if (screenx<screeny) {
		sphererad=0.5*screeny;
	} else {
		sphererad=0.5*screenx;
	}
	glViewport(0,0,width,height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, width/height, 1.0, 100.0 );
	glMatrixMode(GL_MODELVIEW);
}

void displayfunc(void)
{
	int nslices = 256;
	GLfloat m[16],minv[16];
	GLfloat r,dr,z,dz;
	int i;

	glColor4f(1, 1, 1, 1.4/nslices);
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(geomx,geomy,geomz);

	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glTranslatef(.5,.5,.5);

	quatmatrixandinv(orquat,m,minv);
	glMultMatrixf(minv);

	r = -0.5*sqrt(3); 
	dr = sqrt(3)/nslices;
	z = -1.0;
	dz = 2.0/nslices;

	for (i=0;i<nslices;i++) {
		glBegin(GL_QUADS);
			glTexCoord3f(-0.87,-0.87,r); glVertex3f(-1,-1,z); 
			glTexCoord3f(-0.87, 0.87,r); glVertex3f(-1, 1,z); 
			glTexCoord3f( 0.87, 0.87,r); glVertex3f( 1, 1,z); 
			glTexCoord3f( 0.87,-0.87,r); glVertex3f( 1,-1,z); 
		glEnd();
		r += dr;
		z += dz;
		
	}

	glFlush();
	glutSwapBuffers();
}

void mousefunc(int button, int state, int x, int y)
{
	clickx = x;
	clicky = y;

	switch (button) {
		case GLUT_LEFT_BUTTON:
			leftbutstate = (GLUT_DOWN == state) ? 1 : 0;
			break;
		case GLUT_MIDDLE_BUTTON:
			middlebutstate = (GLUT_DOWN == state) ? 1 : 0;
			break;
		case GLUT_RIGHT_BUTTON:
			rightbutstate = (GLUT_DOWN == state) ? 1 : 0;
			break;
	}
}

void motionfunc(int x, int y)
{
	if (leftbutstate) {
		mouserotatemotion(clickx,clicky,x,y);
	} else if (middlebutstate) {
		mousetransmotion(clickx-x,clicky-y);
	} else {
		mousezoommotion(y-clicky);
	}
	clickx=x;clicky=y;
	glutPostRedisplay();
}


/* Main volume rendering code */

/* Load an 8-bit dataset, convert to RGBA, return pointer to RGBA volume */

struct vol3d *load_volume(char *fname,
	unsigned int nx, unsigned int ny, unsigned int nz,
 	uint8_t *rt, uint8_t *gt, uint8_t *bt, uint8_t *at)
{
	struct vol3d *v8bit=NULL,*vrgba=NULL;

	/* Allocate new volume for 8-bit scalar data */
	if (NULL==(v8bit=vol3d_new(nx,ny,nz,VOL3D_UINT8))) {
		fprintf(stderr,"vol3d_new failed\n");
		return NULL;
	}

	/* Load it from file */

	if (0!=vol3d_read_raw(v8bit,fname)) {
		fprintf(stderr,"vol3d_read_raw() failed\n");
		vol3d_destroy(v8bit);
		return NULL;
	}

	if (NULL==(vrgba=vol3d_new_uc2rgba_lookup(v8bit,rt,gt,bt,at))) {
		fprintf(stderr,"vol3d_new_uc2rgba_lookup failed\n");
		vol3d_destroy(v8bit);
		return NULL;
	}

	vol3d_destroy(v8bit);

	return vrgba;
}

/* Initialise lookup tables for transfer function */

void init_tables(void)
{
	int i;
	for (i=0;i<256;i++) {
		if (i<128) {
		r[i] = g[i]=b[i]=a[i]=0;
		} else {
		r[i] = g[i]=b[i]=a[i]=255;
		}
	}
}


/* Load the given volume and upload to graphics card.
 * Return 0 on success.
 */
int load_voltex(char *fname, int nx, int ny, int nz) 
{
	struct vol3d *v=NULL;

	init_tables();

	if (NULL==(v=load_volume(fname,nx,ny,nz,r,g,b,a))) {
		fprintf(stderr,"load_volume() failed\n");
		return -1;
	}
	printf("Loaded %s\n",fname);

	glTexImage3D(
		GL_TEXTURE_3D,
		0, /* level */
		4, /* internalformat (no of channels) */
		nx,ny,nz,
		0, /* border */
		GL_RGBA, /* format */
		GL_UNSIGNED_BYTE, /* type */
		v->data
	);

	vol3d_destroy(v);
	printf("done\n");
	return 0;
}

int main(int argc, char *argv[])
{

	screenx=500;
	screeny=500;
	sphererad = screenx/2;


	/* Initialize GLUT */

	glutInitWindowSize(screenx,screeny);
	glutInitWindowPosition(screenx,0);
	glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);
	glutInit(&argc,argv);

	window = glutCreateWindow("OpenGL window");
	glutDisplayFunc(displayfunc);
	glutReshapeFunc(reshapefunc);
	glutMouseFunc(mousefunc);
	glutMotionFunc(motionfunc);


	/* Initialize GL */

	init_gl();
	orquat = quatrotation(0.0,0.0,0.0,1.0);

	/* Main loop */

	glutMainLoop();

	return 0;
}	




