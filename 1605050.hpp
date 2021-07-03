#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<iostream>
#include <windows.h>
#include <GL/glut.h>
#include<vector>
#include <bits/stdc++.h>
#define pi (2*acos(0.0))
extern const int Near_value = 1;
extern const int Far_value = 1000;
extern int levelOfReflection = 0;
using namespace std;

void drawSquare(GLint x1, GLint y1, GLint x2, GLint y2, GLint x3, GLint y3, GLint x4, GLint y4,int c)
{

if (c == 0)
{
glColor3f(1, 1, 1);
}
else
{
glColor3f(0, 0, 0);

}

// Draw Square
glPushMatrix();
glTranslated(-200,-200,0);
glBegin(GL_QUADS);
glVertex3i(x1, y1,0);
glVertex3i(x2, y2,0);
glVertex3i(x3, y3,0);
glVertex3i(x4, y4,0);
glEnd();
glPopMatrix();
}
void chessboard(int length)
{

GLint x, y;
int dec= 0;
for (x = 0; x <= length; x += 10)
{
for (y = 0; y <= length; y += 10)
{
drawSquare(x, y + 10, x + 10, y + 10, x + 10, y, x, y,dec);
dec = 1 -dec;
}
}
// Process all OpenGL routine s as quickly as possible
glFlush();
}
struct point
{
	double arr[3];
		point(){
	arr[0]=0;
	arr[1]=0;
	arr[2]=0;
	}
	point(double a,double b,double c){
	    arr[0] = a;
	    arr[1] = b;
	    arr[2] = c;
	}
	void print(){
        cout<<"X: "<<arr[0]<<" Y : "<<arr[1]<<" Z: "<<arr[2]<<endl;
    }
	point operator+ (point p) {
	    point res;
	    res.arr[0] = arr[0] + p.arr[0];
	    res.arr[1] = arr[1] + p.arr[1];
	    res.arr[2] = arr[2] + p.arr[2];
	    return res; }
	    point operator*(double k) {
	    point res;
	    res.arr[0] = arr[0] * k;
	    res.arr[1] = arr[1] * k;
	    res.arr[2] = arr[2] * k;
	    return res; }
	    point operator-(point p) {
	    point res;
	    res.arr[0] = arr[0] - p.arr[0];
	    res.arr[1] = arr[1] - p.arr[1];
	    res.arr[2] = arr[2] - p.arr[2];
	    return res; }
	    bool operator==(point p) {
	   if((arr[0] == p.arr[0])&&(arr[1] == p.arr[1])&&(arr[2] == p.arr[2]))
        return true;
        else{
            return false;
        }}
    void normalize()
    {
        double value = sqrt(arr[0] * arr[0] + arr[1] * arr[1] + arr[2] * arr[2]);

        arr[0] /= value;
        arr[0] *= 1.0;
        arr[1] /= value;
        arr[1] *= 1.0;
        arr[2] /= value;
        arr[2] *= 1.0;
    }
};
double dotProduct(point vect_A, point vect_B)
{
    double product = 0;
    product += vect_A.arr[0] * vect_B.arr[0];
    product += vect_A.arr[1] * vect_B.arr[1];
    product += vect_A.arr[2] * vect_B.arr[2];
    product+=1;
    product -=1;
    return product;
}

point crossProduct(point vect_A, point vect_B)
  {
     point cross_P;
     cross_P.arr[0] = vect_A.arr[1] * vect_B.arr[2] - vect_A.arr[2] * vect_B.arr[1];
     cross_P.arr[1] = vect_A.arr[2] * vect_B.arr[0] - vect_A.arr[0] * vect_B.arr[2];
     cross_P.arr[2] = vect_A.arr[0] * vect_B.arr[1] - vect_A.arr[1] * vect_B.arr[0];
     return cross_P;
}
point getR(point L,point N){
        //R = 2 (L.N)N – L
        double t = 2 * dotProduct(L,N);
        point R = (N * t) - L;
        R.normalize();
        return R;
    }
point getReflectedRay(point L,point N){
        //R = 2 (L.N)N – L
        double t = 2 * dotProduct(L,N);
        point R = L - (N * t) ;
        //point R = (N * t) - L;
        R.normalize();
        return R;
    }
    point getReflection(point original_vec, point normal)
{
    double coeff = dotProduct(original_vec, normal) * 2;
    point reflected_vec = original_vec - (normal * coeff);

    reflected_vec.normalize();

    return reflected_vec;
}
class Ray{
    public:
        point start;
        point dir;
        Ray(point p1,point p2){
        start = p1;
        dir = p2;
        dir.normalize();}

};
class Light{
	public:
		point light_pos;
		double color[3];
		Light(point p){
            light_pos = p;
            setColor(255,255,255);
		}
        void setColor(int r,int g,int b){
            color[0] = r;
            color[1] = g;
            color[2] = b;
        }
        void print(){
            cout<<"Light Source"<<endl;
            light_pos.print();
            cout<<"r: "<<color[0]<<" g: "<<color[1]<<" b: "<<color[2]<<endl;
            cout<<endl;
        }
    void draw(){
        glTranslated(0,0,light_pos.arr[2]);
        int slices =25;
        int stacks=25;
		struct point points[100][100];
		int i,j;
		double h,r;
		int length = 1;
		glPushMatrix();
		glTranslated(light_pos.arr[0],light_pos.arr[1],light_pos.arr[2]);
		//generate points
		for(i=0;i<=stacks;i++)
		{
			h=length*sin(((double)i/(double)stacks)*(pi/2));
			r=length*cos(((double)i/(double)stacks)*(pi/2));
			for(j=0;j<=slices;j++)
			{
				points[i][j].arr[0]=r*cos(((double)j/(double)slices)*2*pi);
				points[i][j].arr[1]= r*sin(((double)j/(double)slices)*2*pi);
				points[i][j].arr[2]= h;

			}
		}
		//draw quads using generated points
		for(i=0;i<stacks;i++)
		{
			glColor3f(color[0],color[1],color[2]);
			for(j=0;j<slices;j++)
			{
				glBegin(GL_QUADS);{
					//upper hemisphere
					glVertex3f(points[i][j].arr[0],points[i][j].arr[1],points[i][j].arr[2]);
					glVertex3f(points[i][j+1].arr[0],points[i][j+1].arr[1],points[i][j+1].arr[2]);
					glVertex3f(points[i+1][j+1].arr[0],points[i+1][j+1].arr[1],points[i+1][j+1].arr[2]);
					glVertex3f(points[i+1][j].arr[0],points[i+1][j].arr[1],points[i+1][j].arr[2]);
					//lower hemisphere
					glVertex3f(points[i][j].arr[0],points[i][j].arr[1],-points[i][j].arr[2]);
					glVertex3f(points[i][j+1].arr[0],points[i][j+1].arr[1],-points[i][j+1].arr[2]);
					glVertex3f(points[i+1][j+1].arr[0],points[i+1][j+1].arr[1],-points[i+1][j+1].arr[2]);
					glVertex3f(points[i+1][j].arr[0],points[i+1][j].arr[1],-points[i+1][j].arr[2]);
				}glEnd();
			}
		}
		glPopMatrix();
    }
};
extern vector<Light>lightList;
class Object{
    public:
		point reference_point; // should have x, y, z
		double height, width, length;
		double color[3];
		double coEfficients[4]; // reflection coefficients
		int shine ;// exponent term of specular component
		object(){};
		 virtual void draw(){
            }
		virtual void setColor(){
		cout<<"Checking"<<endl;};
		virtual void setShine(){
		cout<<"Checking"<<endl;};
		virtual void setCoEfficients(){
		cout<<"Checking"<<endl;};
		virtual void print(){
		}
        virtual double intersect(Ray r, double *color, int level){
            return -1.0;
        }
};
extern vector<Object *> objectList;
class Sphere : public Object
{
    public:
    Sphere(point center,double radius){
        reference_point = center;
        length = radius;
        setColor(255,255,255);
    }

    void setColor(int r,int g,int b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void setCoEfficients(double a,double b,double c,double d){
        coEfficients[0] = a ;
        coEfficients[1] = b;
        coEfficients[2] = c;
        coEfficients[3] = d;
    }
    void setShine(double a){
        shine = a;
    }

    void draw()
    {
        glPushMatrix();
        glTranslated(reference_point.arr[0], reference_point.arr[1], reference_point.arr[2]);
        glColor3f(color[0], color[1], color[2]);
        glutSolidSphere(length, 100, 100);
        glPopMatrix();
    }
        void print(){
    cout<<"type : Sphere"<<endl;

        for(int j=0;j<3;j++){
            cout<<reference_point.arr[j]<<" ";
        }
    cout<<endl;
    cout<<"r: "<<color[0]<<" g: "<<color[1]<<" b: "<<color[2]<<endl;
    cout<<"Coefficient: "<<coEfficients[0]<<" "<<coEfficients[1]<<" "<<coEfficients[2]<<" "<<coEfficients[3]<<endl;
    cout<<"Shineness: "<<shine<<endl<<endl;
    }
    point getNormal(point x, point y)
    {
        point temp = x - y;
        temp.normalize();
        return temp;
    }

    double intersect(Ray r, double *final_color, int level){

        double a = 1.0;
        double b = 2 * dotProduct(r.dir,(r.start - reference_point));
        double c = dotProduct((r.start - reference_point),(r.start - reference_point)) - (length * length);

        double d = (b * b) - (4 * a * c);

        if(d < 0)
            return -1;

        d = sqrt(d);
        double t1 = (- b - d) / (2 * a);
        double t2 = (- b + d) / (2 * a);

        double intersect_value =  min(t1, t2);

        if(intersect_value <= 0)
            return -1;

        if((intersect_value < Near_value) || (intersect_value > Far_value) )
            return -1;

        if(level == 0)
            return intersect_value;

        point intersectionPoint = r.start + (r.dir * intersect_value);
        //calculate normal at intersectionPoint
        for(int c = 0 ; c < 3; c++)
            final_color[c] = (color[c] * coEfficients[0]);

        for(int i = 0; i < lightList.size(); i++)
        {
            Light lightSource = lightList.at(i);

            point L = lightSource.light_pos - intersectionPoint;
            L.normalize();
            point N = getNormal(intersectionPoint, reference_point);
            point R = getR(L, N);
            point V = r.start - intersectionPoint;
            V.normalize();

            bool is_obscured = false;

            point start = intersectionPoint + (L * 0.001);
            Ray r1(start,L);
            double * temp2;
            //check if obscured
            for(int j = 0; j < objectList.size(); j++)
            {
                Object * ob = objectList[j];

                double temp = (*ob).intersect(r1,temp2,0);

                if(temp > 0)
                {
                    is_obscured = true;
                    break;
                }
            }

            if(!is_obscured)
            {
                double lambartValue = coEfficients[1] * max(0.0, dotProduct(L, N));
                double phongValue = pow(max(0.0, dotProduct(R, V)), shine) * coEfficients[2];

                for(int c = 0; c < 3; c++){
                    final_color[c] += (lambartValue * color[c]) + (phongValue * 1.0);
                }

            }
        }
        point N = getNormal(intersectionPoint, reference_point);
        point reflect = getReflectedRay(r.dir, N);
        int nearest =-1;
        int tMin,t3;
        if(level < levelOfReflection){
            point p2 = intersectionPoint + reflect ;
            Ray r2(p2,reflect);
            nearest = -1;
            tMin = 1000000;
            double *new_color = new double[3];
            new_color[0] = 0;
            new_color[1] = 0;
            new_color[2] = 0;

            for(int k = 0; k < objectList.size(); k++)
            {
                Object * ob = objectList.at(k);
                t3= (*ob).intersect(r2, new_color, 0);

                if(t3> 0 && t3< tMin){
                    tMin = t3;
                    nearest = k;
                }
            }
            if(nearest != -1)
            {
                t3= objectList.at(nearest)->intersect(r2, new_color, level + 1);

                for(int c = 0; c < 3; c++)
                    final_color[c] += (new_color[c] * coEfficients[3]);
            }

            delete[] new_color;
        }
    }
};



class Triangle : public Object
{
    public:
    point points[3];
    Triangle(point a,point b,point c){
        points[0] = a;
        points[1] = b;
        points[2] = c;
        setColor(255,255,255);
    }
    void print(){
    cout<<"type:triangle"<<endl;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            cout<<points[i].arr[j]<<" ";
        }
        cout<<endl;
    }
    cout<<"r: "<<color[0]<<" g: "<<color[1]<<" b: "<<color[2]<<endl;
    cout<<"Coefficient: "<<coEfficients[0]<<" "<<coEfficients[1]<<" "<<coEfficients[2]<<" "<<coEfficients[3]<<endl;
    cout<<"Shineness: "<<shine<<endl<<endl;
    }
    void setColor(int r,int g,int b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void setCoEfficients(double a,double b,double c,double d){
        coEfficients[0] = a ;
        coEfficients[1] = b;
        coEfficients[2] = c;
        coEfficients[3] = d;
    }
    void setShine(double a){
        shine = a;
    }
    void draw(){
        glColor3f(color[0],color[1],color[2]);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(points[0].arr[0],points[0].arr[1],points[0].arr[2]);
			glVertex3f(points[1].arr[0],points[1].arr[1],points[1].arr[2]);
			glVertex3f(points[2].arr[0],points[2].arr[1],points[2].arr[1]);
        }
        glEnd();
    }
     point getNormal()
    {
        point temp = crossProduct((points[1] - points[0]),(points[2] - points[0]));
        temp.normalize();

        return temp;
    }
    double intersecting_point(Ray ray)
    {
        point v0 = points[0];
        point v1 = points[1];
        point v2 = points[2];

        point edge1 = v1 - v0;
        point edge2 = v2 - v0;

        point h = crossProduct(ray.dir, edge2);
        double a = dotProduct(h, edge1);

        //this ray is parallel to the triangle
        if(a > -0.00001 && a < 0.00001)
            return -1;

        double f = 1.0 / a;
        point s = ray.start - v0;
        double u = f * dotProduct(s, h);

        if(u < 0.0 || u > 1.0)
            return -1;

        point q = crossProduct(s, edge1);
        double v = f * dotProduct(ray.dir, q);

        if(v < 0.0 || u + v > 1.0)
            return -1;

        //now we compute t
        float t = f * dotProduct(edge2, q);
        if(t > 0.00001)
            return t;
        else
            return -1;
    }

    double intersect(Ray ray, double *current_color, int level)
    {
        point v0 = points[0];
        point v1 = points[1];
        point v2 = points[2];

        point edge1 = v1 - v0;
        point edge2 = v2 - v0;

        point h = crossProduct(ray.dir, edge2);
        double a = dotProduct(h, edge1);

        if(a > -0.00001 && a < 0.00001)
            return -1;

        double f = 1.0 / a;
        point s = ray.start - v0;
        double u = f * dotProduct(s, h);

        if(u < 0.0 || u > 1.0)
            return -1;

        point q = crossProduct(s, edge1);
        double v = f * dotProduct(ray.dir, q);

        if(v < 0.0 || u + v > 1.0)
            return -1;

        //now we compute t
        float t = f * dotProduct(edge2, q);

        if(t <= 0)
            return -1;

        if(t < Near_value || t > Far_value)
            return -1;

        if(level == 0)
            return t;

        for(int c = 0; c < 3; c++)
            current_color[c] = color[c] * coEfficients[0];

        //intersection point is => (r0 + t * rd)
        point intersectionPoint = ray.start + (ray.dir * t);
        point normal = getNormal();
        point reflection = getReflection(ray.dir, normal);

        //Illumination
        for(int i = 0; i < lightList.size(); i++)
        {
            point L = lightList.at(i).light_pos - intersectionPoint;
            L.normalize();

            //point start = add(intersectionPoint, multiplyWithScaler(L, EPSILON * 100));
            point start = intersectionPoint + L;
            Ray sunLight(start, L);

            point N = getNormal();
            point R = getR(L, N);

            point V = ray.start - intersectionPoint;
            V.normalize();

            //check if obscured
            bool obscured = false;
            double * hudai;
            for(int j = 0; j < objectList.size(); j++)
            {
                Object * ob = objectList.at(j);
                double temp = (*ob).intersect(sunLight,hudai,0);

                if(temp > 0)
                {
                    obscured = true;
                    break;
                }
            }

            if(!obscured)
            {
                double cosTheta = max(0.0, dotProduct(L, N));
                double cosPhi = max(0.0, dotProduct(R, V));

                double lambart = coEfficients[1] * cosTheta;
                double phong = pow(cosPhi, this->shine) * coEfficients[2];

                for(int c = 0; c < 3; c++)
                    current_color[c] += (lambart * color[c]) + (phong * 1.0);
            }
        }

        int nearest, t_min, t2;

        //Reflection
        if(level < levelOfReflection)
        {
            point start = intersectionPoint + reflection;
            Ray reflectionRay(start, reflection);

            nearest = -1;
            t_min = 1e4;

            double *reflected_color = new double[3];
            reflected_color[0] = reflected_color[1] = reflected_color[2] = 0.0;

            for(int k = 0; k < objectList.size(); k++)
            {

                Object * ob = objectList.at(k);
                double t2 = (*ob).intersect(reflectionRay, reflected_color, 0);

                if(t2 > 0 && t2 < t_min)
                    t_min = t2, nearest = k;
            }

            if(nearest != -1)
            {
                Object * ob = objectList.at(nearest);
                double t2 = (*ob).intersect(reflectionRay,reflected_color,level + 1);


                for(int c = 0; c < 3; c++)
                    current_color[c] += (reflected_color[c] * coEfficients[3]);
            }

            delete[] reflected_color;
        }

        return t;
    }
};

class General : public Object
{
    public:
        double * arr;
    General(double * value,point p,double l,double w,double h){
        arr=value;
        reference_point = p;
        length = l;
        width =w;
        height = h;
        setColor(255,255,255);
    }

    void setColor(int r,int g,int b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void setCoEfficients(double a,double b,double c,double d){
        coEfficients[0] = a ;
        coEfficients[1] = b;
        coEfficients[2] = c;
        coEfficients[3] = d;
    }
    void setShine(double a){
        shine = a;
    }
        void print(){
    cout<<"type : General"<<endl;
        for(int i=0 ; i<10 ; i++){
            cout<<arr[i]<<" ";
        }
        cout<<endl;
        for(int j=0;j<3;j++){
            cout<<reference_point.arr[j]<<" ";
        }
    cout<<length<<" "<<width<<" "<<height<<endl;
    cout<<"r: "<<color[0]<<" g: "<<color[1]<<" b: "<<color[2]<<endl;
    cout<<"Coefficient: "<<coEfficients[0]<<" "<<coEfficients[1]<<" "<<coEfficients[2]<<" "<<coEfficients[3]<<endl;
    cout<<"Shineness: "<<shine<<endl<<endl;
    }
    void draw(){
        cout<<"General "<<endl;
        }
    double intersect(Ray r, double *color, int level){
        return -1;
    }
};


class Floor:public Object{
public:
    double colors[6];
    Floor(int floorWidth,int tileWidth){
        reference_point = point(-floorWidth/2,-floorWidth/2,0);
        length=tileWidth;
        }
    void setColor(int r1,int g1,int b1,int r2,int g2,int b2){
        colors[0] = r1;
        colors[1] = g1;
        colors[2] = b1;
        colors[3] = r2;
        colors[4] = g2;
        colors[5] = b2;
    }

    void setCoEfficients(double a,double b,double c,double d){
        coEfficients[0] = a ;
        coEfficients[1] = b;
        coEfficients[2] = c;
        coEfficients[3] = d;
    }
    void setShine(double a){
        shine = a;
    }

void draw(){

    chessboard(1000);
}
void print(){
    cout<<"Floor: "<<colors[0]<<" "<<colors[2] <<" "<<colors[2]<<" "<<colors[3]<<" "<<colors[4] <<" "<<colors[5]<<endl;
    }

double intersect(Ray r, double *final_color, int level){

       //this->print();
        if(r.dir.arr[2] == 0)
            return -1;

        double t = (-r.start.arr[2] / r.dir.arr[2]);

        if(t <= 0)
            return -1;

        if((t < Near_value) || (t > Far_value) )
            return -1;

        if(level == 0)
            return t;

        point intersectionPoint = r.start + (r.dir * t);

        for(int c = 0 ; c < 3; c++)
            final_color[c] = (color[c] * coEfficients[0]);

        for(int i = 0; i < lightList.size(); i++)
        {
            Light lightSource = lightList.at(i);

            point L = lightSource.light_pos - intersectionPoint;
            L.normalize();
            point N = point(0,0,1);
            point R = getR(L, N);
            point V = r.start - intersectionPoint;
            V.normalize();

            bool is_obscured = false;

            point start = intersectionPoint + (L * 0.001);
            Ray r1(start,L);
            double * temp2;
            //check if obscured
            for(int j = 0; j < objectList.size(); j++)
            {
                Object * ob = objectList[j];

                double temp = (*ob).intersect(r1,temp2,0);

                if(temp > 0)
                {
                    is_obscured = true;
                    break;
                }
            }

            if(!is_obscured)
            {
                double lambartValue = coEfficients[1] * max(0.0, dotProduct(L, N));
                double phongValue = pow(max(0.0, dotProduct(R, V)), shine) * coEfficients[2];
                int xValue = floor(abs(intersectionPoint.arr[0]) / 10.0) ;
		        int yValue = floor(abs(intersectionPoint.arr[1]) / 10.0) ;
		        bool choice = ((xValue + yValue) % 2);
               // cout<<lambartValue<<" "<<phongValue<<endl;
                for(int c = 0; c < 3; c++){
                if(choice){
                    final_color[c] += (lambartValue * colors[c]) + (phongValue * 1.0);
                    }
                else{
                     final_color[c] += (lambartValue * colors[c+3]) + (phongValue * 1.0);
                    }

                }
            }
        }
        point N = point(0,0,1);
        point reflect = getReflectedRay(r.dir, N);
        int nearest =-1;
        int tMin,t3;
        if(level < levelOfReflection){
            point p2 = intersectionPoint + reflect ;
            Ray r2(p2,reflect);
            nearest = -1;
            tMin = 1000000;
            double *new_color = new double[3];
            new_color[0] = 0;
            new_color[1] = 0;
            new_color[2] = 0;

            for(int k = 0; k < objectList.size(); k++)
            {
                Object * ob = objectList.at(k);
                t3= (*ob).intersect(r2, new_color, 0);

                if(t3> 0 && t3< tMin){
                    tMin = t3;
                    nearest = k;
                }
            }
            if(nearest != -1)
            {
                t3= objectList.at(nearest)->intersect(r2, new_color, level + 1);

                for(int c = 0; c < 3; c++)
                    final_color[c] += (new_color[c] * coEfficients[3]);
            }

            delete[] new_color;
        }
    }
};


