//**************************************************************************************
// OPENGL_DRIVER
//**************************************************************************************
#ifndef __OPENGL_DRIVER_H__
#define __OPENGL_DRIVER_H__

#include "fluid.h"
#include <gl/glut.h>
#include "VECTOR_TOOLS.h"

#define NO_MOTION			0
#define ZOOM_MOTION			1
#define ROTATE_MOTION		2
#define TRANSLATE_MOTION	3
#define	SELECT_MOTION		4

using namespace std;

GridAccel body;
//CGLToMovie g_MovieRecorder("Output.Avi", VIEWPORTWIDTH, VIEWPORTHEIGHT);

class OPENGL_DRIVER
{
public:
	static int screen_width, screen_height;
	static float zoom, swing_angle, elevate_angle;
	static float target[3];
	static int motion_mode, mouse_x, mouse_y;

	static float fov, znear, zfar;

	OPENGL_DRIVER(int *argc,char **argv)
	{
		glutInit(argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
		glutInitWindowPosition(50,0);
		glutInitWindowSize(screen_width, screen_height);
		glutCreateWindow ("Rigid Body");
		glutDisplayFunc(Handle_Display);
		glutReshapeFunc(Handle_Reshape);
		glutKeyboardFunc(Handle_Keypress);
		glutMouseFunc(Handle_Mouse_Click);
		glutMotionFunc(Handle_Mouse_Move);
		glutSpecialFunc(Handle_SpecialKeypress);
		glutIdleFunc(Handle_Idle);
		Handle_Reshape(screen_width, screen_height); 
		glutMainLoop();
	}

	~OPENGL_DRIVER()
	{}

	static void Handle_Display()
	{	
		glLoadIdentity();
		glClearColor(0,0,0,0);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

		gluLookAt(0, 0, zoom, 0, 0, 0, 0, 1, 0);
		glRotatef(elevate_angle, 1, 0, 0);
		glRotatef(swing_angle, 0, 1, 0);
		glTranslatef(-target[0], -target[1], -target[2]);

		//Draw the axes
		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
		glColor3f(1,0,0);
		glVertex3f(-0.6, 0, 0);
		glVertex3f(0.6, 0, 0);
		glColor3f(0,1,0);
		glVertex3f(0, -0.6, 0);
		glVertex3f(0, 0.6, 0);
		glColor3f(0,0,1);
		glVertex3f(0, 0, -0.6);
		glVertex3f(0, 0, 0.6);
		glEnd();
		glEnable(GL_LIGHTING);

		body.render();

		//g_MovieRecorder.RecordFrame();	
		glutSwapBuffers();
	}

	static void Handle_Idle()
	{
		if(motion_mode!=SELECT_MOTION)
		{
			body.solve();
		}
		glutPostRedisplay();
	}
	
	static void Handle_Reshape(int w,int h)
	{
		screen_width=w,screen_height=h;
		glViewport(0,0,w, h);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(4, (float)w/(float)h, 0.2, 100);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glEnable(GL_DEPTH_TEST);
		{
			GLfloat LightDiffuse[] = { 1.0, 1.0, 1.0, 1};
			GLfloat LightPosition[] = { 0, 0, -100};
			glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
			glLightfv(GL_LIGHT0, GL_POSITION,LightPosition);
			glEnable(GL_LIGHT0);
		}
		{
			GLfloat LightDiffuse[] = { 1.0, 1.0, 1.0, 1};
			GLfloat LightPosition[] = { 0, 0, 100};
			glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
			glLightfv(GL_LIGHT1, GL_POSITION,LightPosition);
			glEnable(GL_LIGHT1);
		}
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		glutPostRedisplay();
	}
	

	static void Handle_Keypress(unsigned char key,int mousex,int mousey)
	{
		switch(key)
		{
			case 27: 
			break;

		}
		glutPostRedisplay();
	}

	static void Handle_SpecialKeypress(int key, int x, int y)
	{		
			 if(key==100)	swing_angle+=3;
		else if(key==102)	swing_angle-=3;
		else if(key==103)	elevate_angle-=3;
		else if(key==101)	elevate_angle+=3;		
		Handle_Reshape(screen_width, screen_height); 
		glutPostRedisplay();
	}

	static void Get_Selection_Ray(int x, int y, float *_p, float *_q)
	{
		float xx=x*2.0f/screen_width-1.0f;
		float yy=1.0f-y*2.0f/screen_height;	
		
		float m[16], temp;
		glGetFloatv(GL_PROJECTION_MATRIX, m);
		float inv_m[16];
		memset(inv_m, 0, sizeof(float)*16);
		inv_m[0]=1/m[0];
		inv_m[5]=1/m[5];
		inv_m[14]=1/m[14];
		inv_m[11]=-1;
		inv_m[15]=m[10]/m[14];

		float p[4]={xx, yy, 0, 1};
		float q[4]={xx, yy, 1, 1};
		float pp[4], qq[4];

		Matrix_Times_Vector_4(inv_m, p, pp);
		Matrix_Times_Vector_4(inv_m, q, qq);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(target[0], target[1], target[2]);
		glRotatef(-swing_angle, 0, 1, 0);
		glRotatef(-elevate_angle, 1, 0, 0);
		glTranslatef(0, 0, zoom);
		float inv_n[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, inv_n);

		Matrix_T_Times_Vector_4(inv_n, pp, p);
		Matrix_T_Times_Vector_4(inv_n, qq, q);

		Homogeneous_Projection(p, _p);
		Homogeneous_Projection(q, _q);
	}

	static void Handle_Mouse_Move(int x, int y)
	{
		if( motion_mode!=NO_MOTION)
		{
			if(motion_mode==ROTATE_MOTION) 
			{
				swing_angle   += (float)(x - mouse_x)*360/(float)screen_width;
				elevate_angle += (float)(y - mouse_y)*180/(float)screen_height;
		        if     (elevate_angle> 90)	elevate_angle =  90;
				else if(elevate_angle<-90)	elevate_angle = -90;
			}
			if(motion_mode==ZOOM_MOTION)	zoom+=0.05 * (y-mouse_y);
			if(motion_mode==TRANSLATE_MOTION)
			{
				target[0] -= 0.1*(mouse_x - x);
				target[2] += 0.1*(mouse_y - y);
			}
			mouse_x=x;
			mouse_y=y;
			glutPostRedisplay();
		}
	}

	static void Handle_Mouse_Click(int button, int state, int x, int y)
	{
		if(state==GLUT_UP) motion_mode=NO_MOTION;
		else
		{
			float p[3], q[3];
			Get_Selection_Ray(x, y, p, q);

				int modif = glutGetModifiers();
				if (modif & GLUT_ACTIVE_SHIFT)		motion_mode = ZOOM_MOTION;
				else if (modif & GLUT_ACTIVE_CTRL)	motion_mode = TRANSLATE_MOTION;
				else								motion_mode = ROTATE_MOTION;
				mouse_x=x;
				mouse_y=y;
		}
		glutPostRedisplay();
	}
};


int		OPENGL_DRIVER::screen_width=800;
int		OPENGL_DRIVER::screen_height=800;

// 3D Display Configuration
float	OPENGL_DRIVER::zoom=20;//10
float	OPENGL_DRIVER::swing_angle=-0; //-25;
float	OPENGL_DRIVER::elevate_angle=0; //0;
float	OPENGL_DRIVER::target[3]={0,0,0};//0,0,0
int		OPENGL_DRIVER::motion_mode=NO_MOTION;
int		OPENGL_DRIVER::mouse_x=0;
int		OPENGL_DRIVER::mouse_y=0;



#endif
