#include "bitmap_image.hpp"
#include <iostream>
#include <string>
#include <fstream>

#include "1605050.hpp"
#include "draw.hpp"

double cameraHeight;
double cameraAngle;
int totalObject=-1;
double angle;
bool flag = false;
 vector<Light> lightList;
 vector<Object *> objectList;

point pos,u,r,l,center;

double k =2.0 ;
double a = 4;
double b,d,e;
extern int levelOfReflection;
int pixelCount;




point point_init(string s)
{
    stringstream ss(s);
    double word;
    vector<double> values;
    while (ss >> word) {
        values.push_back(word);
    }
   point temp(values[0],values[1],values[2]);
   return temp;
}
double *getColor(string str){
    double * arr =new double[3];
    stringstream ss(str);
    ss >> arr[0];
    ss >> arr[1];
    ss >> arr[2];
    //cout<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<endl;
    return arr;
}

double * getCoeff(string str){
    double * arr =new double[4];
    stringstream ss(str);
    ss >> arr[0];
    ss >> arr[1];
    ss >> arr[2];
    ss>> arr[3];
    //cout<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<" "<<arr[3]<<endl;
    return arr;
}

double * getGen2nd(string str){
    double * arr =new double[6];
    stringstream ss(str);
    for(int i=0;i<6;i++){
        ss >> arr[i];
    }
    //cout<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<" "<<arr[3]<<endl;
    return arr;
}
double * getGenCoeff(string str){
    double * arr =new double[10];
    stringstream ss(str);
    for(int i=0;i<10;i++){
        ss>>arr[i];
    }
    //cout<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<" "<<arr[3]<<endl;
    return arr;
}
void loadData(){

    cout<<"LoadData called"<<endl;
	ifstream MyReadFile("E:\offline/CSE 410/offline 3/Offline/ray Tracing/scene.txt",ios::in);
	string myText;
	stringstream ss(myText);
	int lineCount = 1;
	if(MyReadFile.is_open()){
    while (getline (MyReadFile, myText) && !flag) {

		if(lineCount == 1){
			ss<<myText;
			ss >> levelOfReflection;
			cout<<levelOfReflection<<endl;

		}
		else if(lineCount == 2){
            stringstream(myText) >> pixelCount;
			cout<<pixelCount<<endl;

		}
		else if(lineCount == 4){
            stringstream(myText) >> totalObject;
			//cout<<totalObject<<endl;
		}
       myText = std::regex_replace(myText, std::regex("^ +| +$|( ) +"), "$1");

        if(myText== "sphere"){
         getline (MyReadFile, myText);
         point center = point_init(myText);
          getline (MyReadFile, myText);
         double radius;
         stringstream(myText) >> radius;
         Sphere * s = new Sphere(center,radius);
         getline (MyReadFile, myText);
          double * arr;
          arr = getColor(myText);
          (*s).setColor(arr[0],arr[1],arr[2]);
          getline(MyReadFile,myText);
          arr = getCoeff(myText);
          (*s).setCoEfficients(arr[0],arr[1],arr[2],arr[3]);
          getline(MyReadFile,myText);
          int shine;
          stringstream(myText)>>shine;
          (*s).setShine(shine);
          (*s).print();
          objectList.push_back(s);
          totalObject--;
        }
        if(myText == "triangle"){
            vector<point> temp;
            double * arr;
            for(int i=0;i<3;i++){
            getline (MyReadFile, myText);
                temp.push_back(point_init(myText));
            }
            Triangle * tri = new Triangle(temp[0],temp[1],temp[2]);
            getline (MyReadFile, myText);
            arr= getColor(myText);
            //cout<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<endl;
            (*tri).setColor(arr[0],arr[1],arr[2]);
            getline(MyReadFile,myText);
            arr = getCoeff(myText);
            (*tri).setCoEfficients(arr[0],arr[1],arr[2],arr[3]);
            getline(MyReadFile,myText);
            int shine;
            stringstream(myText)>>shine;
            (*tri).setShine(shine);
            (*tri).print();
            objectList.push_back(tri);
            totalObject--;
            }
        if(myText == "general")
        {
            double * values;
            double * a;
            getline(MyReadFile,myText);
            values = getGenCoeff(myText);
            getline(MyReadFile,myText);
            a = getGen2nd(myText);
            point p(a[0],a[1],a[2]);
            General * g  = new General(values,p,a[3],a[4],a[5]);
             getline (MyReadFile, myText);
            values= getColor(myText);
            (*g).setColor(values[0],values[1],values[2]);
            getline(MyReadFile,myText);
            values = getCoeff(myText);
            (*g).setCoEfficients(values[0],values[1],values[2],values[3]);
            getline(MyReadFile,myText);
            int shine;
            stringstream(myText)>>shine;
            (*g).setShine(shine);
            objectList.push_back(g);
            (*g).print();
            totalObject--;
        }
        cout<<"TotalObject: "<<totalObject<<endl;
        if(totalObject == 0){

            getline(MyReadFile,myText);
            getline(MyReadFile,myText);
            cout<<"Text: "<<myText<<endl;
            int lightSource;
            stringstream(myText) >> lightSource;
            cout<<"LightSource: "<<lightSource<<endl;
            for(int i=0;i<lightSource;i++){
                getline(MyReadFile,myText);
                Light l(point_init(myText));
                getline(MyReadFile,myText);
                double * arr = getColor(myText);
                l.setColor(arr[0],arr[1],arr[2]);
                l.print();
                lightList.push_back(l);
            }
            flag = true;
        }
        lineCount++;
	}
	}else{
	cout<<"Error in opening file"<<endl;}
    Floor * temp = new Floor(1000,20);
    (*temp).setColor(0,0,0,1,1,1);
    (*temp).setCoEfficients(0.4,0.2,0.1,0.3);
    (*temp).setShine(5);
    objectList.push_back(temp);
}



point rotateAxis(point k,point v,double a){
    point res;
    res = crossProduct(k,v) * sin(pi/180*a);
    res = res + v * cos(pi/180*a);
    return res;
}

void capture(){
    bitmap_image image(pixelCount,pixelCount);

    for(int i=0;i<pixelCount;i++){
        for(int j=0;j<pixelCount;j++){
            image.set_pixel(i,j,0,0,0);
        }
    }
    double planeDistance = (500/2.0) / tan(3.1416*80/360);
    l.print(); r.print(); u.print(); pos.print();

    point topleft = (((pos + (l * planeDistance)) + (u * 250)) - (r*250));
    topleft.print();
    double du = (500.0/pixelCount)*1.0;
    double dv = (500.0/pixelCount)*1.0;
    topleft = topleft + r*(0.5*du) - u*(0.5*dv);
    //topleft.print();
    cout<< planeDistance<<" | "<<du<<"  | "<<dv<<" | "<<pixelCount<< endl;
    int nearest;
    double t, tMin;
    point corner , a;
    double *dummy_color = new double[3];
    for(int i=0;i<pixelCount;i++){
        for(int j=0;j<pixelCount;j++){
            a = (r * (i * du)) - (u * (j * dv));
            corner = topleft + a;
            a = corner - pos;
            Ray ray(pos, a);
            tMin =10000;
            nearest = -1;
            for(int k=0;k<objectList.size();k++){
                Object * ob = objectList.at(k);
                t = (*ob).intersect(ray,dummy_color,0);
                //cout<<t<<endl;
                if(t > 0 && t < tMin){
                    tMin = t;
                    nearest = k;
                }
            }
            if(nearest != -1)
            {
                //cout<<nearest<<endl;
                Object * ob = objectList.at(nearest);
                t  = (*ob).intersect(ray,dummy_color,1);
                //t = (*ob).intersect(ray,dummy_color,1);
                if(nearest == 2){
                    //cout<<dummy_color[0]<<" "<<dummy_color[1]<<endl;
                }
                for(int c = 0; c < 3; c++)
                {
                    if(dummy_color[c] < 0.0)
                        dummy_color[c] = 0.0;

                    else if(dummy_color[c] > 1.0)
                        dummy_color[c] = 1.0;
                }
            }

            else
            {
                dummy_color[0] = 0;
                dummy_color[1] = 0;
                dummy_color[2] = 0;
            }
            image.set_pixel(i, j, 255 * dummy_color[0], 255 * dummy_color[1], 255 * dummy_color[2]);
        }
    }


    image.save_image("E:\offline/CSE 410/offline 3/Offline/ray Tracing/test.bmp");
    cout<<"done"<<endl;
}
void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
		    r = rotateAxis(u,r,a);
		    l = rotateAxis(u,l,a);
		    //printf("%f %f %f",crossProduct(r,l).x,crossProduct(r,l).y,crossProduct(r,l).z);
			break;
        case '2':
            r = rotateAxis(u,r,-a);
		    l = rotateAxis(u,l,-a);
			break;
        case '3':
            l = rotateAxis(r,l,a);
		    u = rotateAxis(r,u,a);
			break;
        case '4':
             l = rotateAxis(r,l,-a);
		    u = rotateAxis(r,u,-a);
			break;
        case '5':
             r = rotateAxis(l,r,-a);
		    u = rotateAxis(l,u,-a);
			break;
        case '6':
		     r = rotateAxis(l,r,a);
		    u = rotateAxis(l,u,a);
			break;
        case '0':
            capture();
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;

                pos = pos - l * k;
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;

                pos = pos + l * k;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;

			pos = pos + r * k;
			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;

                pos = pos - r * k;
			break;

		case GLUT_KEY_PAGE_UP:

                pos = pos + u * k;
			break;
		case GLUT_KEY_PAGE_DOWN:

			pos = pos - u*k;
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP



			}
			break;

		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP

			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}




void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
    gluLookAt(pos.arr[0], pos.arr[1], pos.arr[2], pos.arr[0] + l.arr[0], pos.arr[1] + l.arr[1], pos.arr[2] + l.arr[2], u.arr[0], u.arr[1], u.arr[2]);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects


    //glTranslated(20,0,0);
    point p (0,0,0);
    for(int i=0;i<objectList.size();i++){
     Object * ptr = objectList.at(i);
    (*ptr).draw();
    }
    //cout<<lightList.size()<<endl;
    Light lig = lightList.at(0);
    lig.draw();

    //sp1.draw(10,10);
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}



void init(){
	//codes for initialization

    loadData();
   // capture();
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	Near_value,	Far_value);

	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
	pos.arr[0] = 100;
	pos.arr[1] = 100;
	pos.arr[2] = 0;
    u.arr[0] = 0;
    u.arr[1] = 0;
    u.arr[2] = 1;
    r.arr[0] = -(1.0/sqrt(2));
    r.arr[1] = (1.0 / sqrt(2));
    r.arr[2] = 0;
    l.arr[0] = -(1.0/sqrt(2));
    l.arr[1] = -(1.0/sqrt(2));
    l.arr[2] = 0;

}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}

