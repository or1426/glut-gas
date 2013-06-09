/*THE NUMBER OF CALCULATIONS SCALES WITH N^2 SO IT CAN SLOW DOWN FAST*/

#define _PI  3.1415926535

#define _FOREVER -1

#include <iostream>
#include <sstream>

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <chrono>

#include <GL/gl.h>
#include <GL/freeglut.h>

using namespace std;

class vect3d
{
private:
    double array[3];
public:
    vect3d()
    {
        array[0] = 0;
        array[1] = 0;
        array[2] = 0;
    }
    vect3d(double d1,double d2,double d3)
    {
        array[0] = d1;
        array[1] = d2;
        array[2] = d3;
    }

    double& operator[](int i)
    {
        if(i<3)
        {
            return array[i];
        }
        else
        {
            return array[0];
        }
    }

    vect3d operator +(vect3d param)
    {
        return vect3d(array[0] + param[0],array[1] + param[1],array[2] + param[2]);
    }

    void operator +=(vect3d param)
    {
        (*this) = vect3d(array[0] + param[0],array[1] + param[1],array[2] + param[2]);
    }

    vect3d operator - (vect3d param)
    {
        return vect3d(array[0] - param[0],array[1] - param[1],array[2] - param[2]);
    }

    vect3d operator *(double param)
    {
        return vect3d(array[0] * param,array[1] * param,array[2] * param);
    }
    vect3d operator /(double param)
    {
        return (*this)*(1/param);
    }

    double dot(vect3d param)
    {
        return (array[0]*param[0]) + (array[1]*param[1]) + (array[2]*param[2]);
    }
    double mod()
    {
        return sqrt((array[0]*array[0]) + (array[1]*array[1]) + (array[2]*array[2]));
    }
};

class particle
{
public:
    vect3d location;
    vect3d velocity;
    double radius;
    double mass;
    float R,G,B;
    float KE;
    GLfloat scaledX, scaledY, scaledZ,scaledRadius;
    int hasColidedWith;
    particle(double radius,vect3d location, vect3d velocity,double mass, float R, float G, float B)
    {
        (*this).radius = radius;
        (*this).R = R;
        (*this).G = G;
        (*this).B = B;
        (*this).location = location;
        (*this).velocity = velocity;
        (*this).mass = mass;
        hasColidedWith = -1;
    }

    particle(vect3d location)
    {
        *this = particle(1.0,location, vect3d(0.0, 0.0, 0.0), 1.0, 0.0, 0.1, 1.0);

    }

    particle()
    {
        *this = particle(1.0,vect3d(0.0, 0.0, 0.0), vect3d(0.0, 0.0, 0.0), 1.0, 0.0, 0.1, 1.0);

    }

    void draw(double scale,vect3d centre)
    {
        glColor3f(R,G,B);
        scaledX = ((location[0] - centre[0])/scale) ;
        scaledY = ((location[1] - centre[1])/scale) ;
        scaledZ = ((location[2] - centre[2])/scale) ;
        scaledRadius = radius/scale;
        int numPoints = 20;
        glBegin(GL_TRIANGLE_FAN);
        for(int i = 0; i<numPoints; ++i)
        {
            glVertex3f(scaledX + scaledRadius*cos(2.0*_PI*i/(float)numPoints),scaledY + scaledRadius*sin(2.0*_PI*i/(float)numPoints),scaledZ);

        }
        glEnd();
    }

    float updateKE()
    {
        KE = 0.5 * mass * velocity.dot(velocity);
        return KE;
    }

    void setColour(float red,float green,float blue)
    {
        R = red;
        G = green;
        B = blue;
    }


    particle operator *(double param)
    {
        return particle(radius*param,location*param,velocity*param,mass*param,R,G,B);
    }
    particle operator +(particle param)
    {
        return particle(radius+ param.radius,location +param.location,velocity+param.velocity,mass+param.mass,R,G,B);
    }

};

class world
{
private:
    particle * particleArray;
    double scale;
    vect3d center;
    int numParticles;
    double ** wallContribution;

    std::chrono::system_clock::time_point nowTime;
    std::chrono::system_clock::time_point prevTime;
    std::chrono::duration<int,std::milli> waitTime;

    void renderString(float x, float y, void *font, const char* string)
    {
        glColor3f(0.0, 0.0, 0.0);
        glRasterPos2f(x, y);
        glutBitmapString(font, (unsigned char *)string);
    }

    void drawAxis()
    {
        /*
        ostringstream lables[2][2];
        ostringstream origin;

        glColor3f(0.0,0.0,0.0);
        glLineWidth(3);

        glBegin(GL_LINES);
        glVertex2f(-1.00,0.00);
        glVertex2f( 1.00,0.00);

        glVertex2f( 0.00,-1.00);
        glVertex2f( 0.00, 1.00);

        glVertex2f(-1.00, 0.05);
        glVertex2f(-1.00,-0.05);

        glVertex2f( 1.00, 0.05);
        glVertex2f( 1.00,-0.05);

        glVertex2f(-0.05,-1.00);
        glVertex2f( 0.05,-1.00);

        glVertex2f(-0.05, 1.00);
        glVertex2f( 0.05, 1.00);
        glEnd();

        lables[0][0] << xMin;
        lables[0][1] << xMax;
        lables[1][0] << yMin;
        lables[1][1] << yMax;

        origin << '(' << 0.5*(xMax + xMin) << ',' << 0.5*(yMax + yMin) << ')';

        renderString(-1.0                , 0.05, GLUT_BITMAP_TIMES_ROMAN_24, lables[0][0].str().c_str());
        renderString(0.98-0.005*log(xMax), 0.05, GLUT_BITMAP_TIMES_ROMAN_24, lables[0][1].str().c_str());
        renderString(0.0                 , -1.0, GLUT_BITMAP_TIMES_ROMAN_24, lables[1][0].str().c_str());
        renderString(0.0                 , 0.95, GLUT_BITMAP_TIMES_ROMAN_24, lables[1][1].str().c_str());

        renderString(0.0                 , 0.05, GLUT_BITMAP_TIMES_ROMAN_24, origin.str().c_str());
        */
    }

    void updateScale()
    {
        center = vect3d(0,0,0);
        for(int i = 0; i<numParticles; ++i)
        {
            center += particleArray[i].location;

        }
        center = center/(double)numParticles;
        for(int i = 0; i<numParticles; ++i)
        {
            for(int j = 0; j<3; ++j)
            {
                if(fabs(particleArray[i].location[j] + particleArray[i].radius - center[j]) > scale)
                {
                    scale = 1.5*fabs(particleArray[i].location[j] + particleArray[i].radius- center[j]);
                }
            }
        }

    }

    int detectCollision(particle p)
    {
        for(int i = 0; i<numParticles; ++i)
        {
            if((p.location - particleArray[i].location).mod() < p.radius+particleArray[i].radius)
            {
                return -1;
            }
        }
        return 0;
    }

public:

    world()
    {
        *this = *new world(800,800);
    }

    world(int xSize,int ySize)
    {
        char *myargv [1];
        int myargc=1;
        myargv [0]=strdup ("Graph");
        glutInit(&myargc, myargv);

        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
        glutInitWindowSize(xSize,ySize);

        glutCreateWindow(" ");

        glMatrixMode ( GL_PROJECTION );
        glClearColor ( 1, 1, 1, 0.0 );
        glClear(GL_COLOR_BUFFER_BIT);
        glPointSize(6.0);

        (*this).scale = 0;
        numParticles = 0;

        particleArray = (particle *)malloc(0);
        waitTime = std::chrono::duration<int,std::milli>(20);
        wallContribution = (double **)malloc(2*sizeof(double));
        for(int i = 0; i<2; ++i)
        {
            wallContribution[i] = (double *)malloc(3*sizeof(double));
        }

        wallContribution[0][0] = 10;
        wallContribution[1][0] = 80;
        wallContribution[0][1] = 10;
        wallContribution[1][1] = 10;
        wallContribution[0][2] = 10;
        wallContribution[1][2] = 10;
    }

    int addParticle(particle point)
    {
        if(detectCollision(point) == 0)
        {
            if(scale < 5*point.radius)
            {
                scale = 5*point.radius;
            }
            ++numParticles;
            particleArray = (particle *)realloc(particleArray,numParticles * sizeof(particle));
            particleArray[numParticles - 1] = point;
            updateScale();
            return 0;
        }
        return -1;
    }

    void display()
    {
        float KE_max = 0;
        for(int i = 0; i<numParticles; ++i)
        {
            if(KE_max < particleArray[i].updateKE())
            {
                KE_max = particleArray[i].updateKE();
            }
        }

        for(int i = 0; i<numParticles; ++i)
        {
            particleArray[i].setColour(particleArray[i].KE/KE_max,0, 1 - particleArray[i].KE/KE_max);
        }
        glClear(GL_COLOR_BUFFER_BIT);
        //updateScale();
        drawAxis();
        for(int i = 0; i<numParticles; ++i)
        {
            particleArray[i].draw(scale,center);
        }
        glFlush();
        glutSwapBuffers();
    }

    void updateParticles(double h)
    {
        for(int i = 0; i<numParticles; ++i)
        {
            particleArray[i].location += particleArray[i].velocity*h;
            for(int j = 0; j<numParticles; ++j)
            {
                if(i != j && particleArray[i].hasColidedWith != j)
                {
                    if( (particleArray[i].location - particleArray[j].location).mod() < particleArray[i].radius+particleArray[j].radius )
                    {
                        particleArray[j].hasColidedWith = i;
                        vect3d relativeVelocity = particleArray[i].velocity - particleArray[j].velocity;
                        vect3d unitVect = (particleArray[i].location - particleArray[j].location)/(particleArray[i].location - particleArray[j].location).mod();
                        double parallelComponent = unitVect.dot(relativeVelocity);
                        particleArray[i].velocity = particleArray[i].velocity - unitVect*parallelComponent;
                        particleArray[j].velocity = particleArray[j].velocity + unitVect*parallelComponent;
                    }
                }
            }
            for(int j = 0; j<3; ++j)
            {
                if(particleArray[i].location[j] + particleArray[i].radius > center[j] + scale)
                {
                    particleArray[i].velocity[j] = - fabs((particleArray[i].velocity[j] + wallContribution[0][j])/(double)2);
                }
                else if(particleArray[i].location[j] - particleArray[i].radius < center[j] - scale)
                {
                    particleArray[i].velocity[j] =  fabs( (particleArray[i].velocity[j] - wallContribution[1][j])/(double)2);
                }
            }
            particleArray[i].hasColidedWith = -1;
        }
    }

    void animate(double h,double t_max)
    {
        /*The chrono shit makes this loads faster! It looks like drawing is slooooooooow*/
        nowTime = std::chrono::system_clock::now();
        prevTime = nowTime;
        if(t_max == _FOREVER)
        {
            while(1)
            {
                nowTime = std::chrono::system_clock::now();
                if(nowTime - prevTime > waitTime)
                {
                    prevTime = nowTime;
                    display();
                }
                updateParticles(h);
            }
        }
        else
        {
            for(double t = 0; t<t_max; t += h)
            {
                nowTime = std::chrono::system_clock::now();
                if(nowTime - prevTime > waitTime)
                {
                    prevTime = nowTime;
                    display();
                }
                updateParticles(h);
            }
        }
    }
};




