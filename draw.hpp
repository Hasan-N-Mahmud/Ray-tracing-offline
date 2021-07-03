void drawFunel(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=2*radius-radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].arr[0]=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].arr[1]=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].arr[2]=h;
		}
	}
	//draw quads using generated points
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
                glColor3f(0,0,0);
            stripe = true;
        }
			     glVertex3f(points[i][j].arr[0],points[i][j].arr[1],-points[i][j].arr[2]);
				glVertex3f(points[i][j+1].arr[0],points[i][j+1].arr[1],-points[i][j+1].arr[2]);
				glVertex3f(points[i+1][j+1].arr[0],points[i+1][j+1].arr[1],-points[i+1][j+1].arr[2]);
				glVertex3f(points[i+1][j].arr[0],points[i+1][j].arr[1],-points[i+1][j].arr[2]);
			}glEnd();
		}
	}
}

void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].arr[0]=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].arr[1]=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].arr[0],points[i].arr[1],0);
			glVertex3f(points[i+1].arr[0],points[i+1].arr[1],0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].arr[0]=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].arr[1]=radius*sin(((double)i/(double)segments)*2*pi);

    }
    //draw triangles using generated points
    for(i=0;i<1;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(0.5,1,0);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].arr[0],points[i].arr[1],0);
			glVertex3f(points[i+1].arr[0],points[i+1].arr[1],1);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].arr[0]=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].arr[1]=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].arr[2]=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
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
}

void drawCylinder(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius;
		for(j=0;j<=slices;j++)
		{
			points[i][j].arr[0]=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].arr[1]=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].arr[2]=h;
		}
	}
	//draw quads using generated point
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);
			{
                //lower hemisphere
                if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
            glColor3f(0,0,0);
            stripe = true;
        }
                glVertex3f(points[i][j].arr[0],points[i][j].arr[1],points[i][j].arr[2]);
				glVertex3f(points[i][j+1].arr[0],points[i][j+1].arr[1],points[i][j+1].arr[2]);
				glVertex3f(points[i+1][j+1].arr[0],points[i+1][j+1].arr[1],points[i+1][j+1].arr[2]);
				glVertex3f(points[i+1][j].arr[0],points[i+1][j].arr[1],points[i+1][j].arr[2]);

				glVertex3f(points[i][j].arr[0],points[i][j].arr[1],-points[i][j].arr[2]);
				glVertex3f(points[i][j+1].arr[0],points[i][j+1].arr[1],-points[i][j+1].arr[2]);
				glVertex3f(points[i+1][j+1].arr[0],points[i+1][j+1].arr[1],-points[i+1][j+1].arr[2]);
				glVertex3f(points[i+1][j].arr[0],points[i+1][j].arr[1],-points[i+1][j].arr[2]);
			}glEnd();
		}
	}
}

void drawUpperSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].arr[0]=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].arr[1]=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].arr[2]=h;
		}
	}
	//draw quads using generated points
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
            glColor3f(0,0,0);
            stripe = true;
        }
			    //upper hemisphere
				glVertex3f(points[i][j].arr[0],points[i][j].arr[1],points[i][j].arr[2]);
				glVertex3f(points[i][j+1].arr[0],points[i][j+1].arr[1],points[i][j+1].arr[2]);
				glVertex3f(points[i+1][j+1].arr[0],points[i+1][j+1].arr[1],points[i+1][j+1].arr[2]);
				glVertex3f(points[i+1][j].arr[0],points[i+1][j].arr[1],points[i+1][j].arr[2]);

			}glEnd();
		}
	}
}

void drawLowerSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].arr[0]=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].arr[1]=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].arr[2]=h;
		}
	}
	//draw quads using generated points
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);
			{
                //lower hemisphere
                if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
            glColor3f(0,0,0);
            stripe = true;
        }
                glVertex3f(points[i][j].arr[0],points[i][j].arr[1],-points[i][j].arr[2]);
				glVertex3f(points[i][j+1].arr[0],points[i][j+1].arr[1],-points[i][j+1].arr[2]);
				glVertex3f(points[i+1][j+1].arr[0],points[i+1][j+1].arr[1],-points[i+1][j+1].arr[2]);
				glVertex3f(points[i+1][j].arr[0],points[i+1][j].arr[1],-points[i+1][j].arr[2]);
			}glEnd();
		}
	}
}
