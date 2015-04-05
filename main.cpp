/*
 * main.cpp
 *
 * 2D
 *  Created on: Jul 20, 2010
 *      Author: jean-philippePeraud
 
 
 * applies the linearized method to the periodic problem
 * terminates the particles after a fixed number of collision events (chosen)
 * 
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <math.h>
#include "randomClass.h"


using namespace std;





string date()
{
	time_t     now;
		    struct tm  *ts;
		    char       buf[80];

		    /* Get the current time */
		    now = time(NULL);

		    /* Format and print the time, "ddd yyyy-mm-dd hh:mm:ss zzz" */
		    ts = localtime(&now);
		    strftime(buf, sizeof(buf), "%m%d%H%M", ts);
		    return buf;
		   // printf("%s\n", buf);
}

int choice_omega_norm(double norm_cumul[],int N_interv,RandomClass r)
{
	//assign a frequency
	
	
	double R=(double)rand()/(double)RAND_MAX;//r.randu();
	int i1 = 0;
	int i2 = (N_interv+1)/2;
	int i3 = N_interv+1;
	int count = 0;
	while (i2-i1>0){
	    count=count+1;
	    if (R<norm_cumul[i2-1]){
	        i3 = i2;
	        i2 = (i2+i1)/2;}
	    else{
	        i1 = i2;
	        i2 = (i2+i3)/2;
	    }
	//    cout << norm_cumul[i2-1] << ": " << i2 << endl;
	    
	}
//	cout << R << " " << norm_cumul[i2-1] << ": " << i2 << endl;
	return i2;
}

void calcul_contribution_MC_elem(double* temp_grid, double Xo1, double Yo1, double X11, double Y11, double Vx, double Vy, int sign1, double elapsed_time, double Peff, grid my_grid, int ncellX, int ncellY, double LX, double LY, int jX, Parameters Par)
{
	
	
	
	double DeltaX=LX/ncellX;
	double DeltaY=LY/ncellY;
	int indeXo=(int)floor(Yo1/DeltaY);
	int indeX1=(int)floor(Y11/DeltaY);
	double time;
	//cout << "Yo " << indeXo << endl;
//	cout << "Y1 " << Y11 << DeltaY << indeX1 << endl;
	if (indeX1==ncellY){
		indeX1=ncellY-1;
	}
	if (indeXo==ncellY){
		indeXo=ncellY-1;
	}
	
	if (indeXo==indeX1){
//		cout << temp_grid[jX+indeXo*ncellX] << endl;
		if (Yo1==Y11 && Xo1==X11){
			temp_grid[jX+indeXo*ncellX]=temp_grid[jX+indeXo*ncellX]+Peff*sign1*elapsed_time/my_grid.get_area(jX,indeXo,Par);
	    }else{
		    temp_grid[jX+indeXo*ncellX]=temp_grid[jX+indeXo*ncellX]+Peff*sign1*abs(Y11-Yo1)/abs(Vy)/my_grid.get_area(jX,indeXo,Par);
		}
//		cout << Peff*sign1*abs(Y11-Yo1)/abs(Vy)/my_grid.get_area(jX,indeXo,Par) << " " << temp_grid[jX+indeXo*ncellX] << " " << Peff << " " << Y11-LY/ncellY*indeX1 << endl;
		
		
	}else{
		if (indeX1>indeXo){
            for (int j=indeXo+1; j<indeX1; j++){
                temp_grid[jX+j*ncellX]=temp_grid[jX+j*ncellX]+Peff*(LY/ncellY)*sign1/abs(Vy)/my_grid.get_area(jX,j,Par);
		    }
			
			temp_grid[jX+indeXo*ncellX]=temp_grid[jX+indeXo*ncellX]+Peff*sign1*(LY/ncellY*(indeXo+1)-Yo1)/abs(Vy)/my_grid.get_area(jX,indeXo,Par);
			temp_grid[jX+indeX1*ncellX]=temp_grid[jX+indeX1*ncellX]+Peff*sign1*(Y11-LY/ncellY*indeX1)/abs(Vy)/my_grid.get_area(jX,indeX1,Par);
//			cout << Peff*Vx*sign1*(LY/ncellY*(indeXo+1)-Yo1)/abs(Vy)/DeltaY/my_grid.get_area(jX,indeXo,Par) << " " << temp_grid[jX+indeXo*ncellX] << " " << Vy << " " << Y11-LY/ncellY*indeX1 << endl;
//			cout << jX << " " << indeXo <<" " << my_grid.get_area(jX,indeXo,Par) << endl;
	    }
	    if(indeX1<indeXo){
            for (int j=indeX1+1; j<indeXo; j++){
                temp_grid[jX+j*ncellX]=temp_grid[jX+j*ncellX]+Peff*LY/ncellY*sign1/abs(Vy)/my_grid.get_area(jX,j,Par);
		    }
		    temp_grid[jX+indeX1*ncellX]=temp_grid[jX+indeX1*ncellX]+Peff*sign1*(LY/ncellY*(indeX1+1)-Yo1)/abs(Vy)/my_grid.get_area(jX,indeX1,Par);
		    temp_grid[jX+indeXo*ncellX]=temp_grid[jX+indeXo*ncellX]+Peff*sign1*(Y11-LY/ncellY*indeXo)/abs(Vy)/my_grid.get_area(jX,indeXo,Par);
//			cout << Peff*sign1*abs(Y11-Yo1)/abs(Vy)/my_grid.get_area(jX,indeXo,Par) << " " << temp_grid[jX+indeXo*ncellX] << " " << Vy << "sign " << sign1 << endl;
			if (isnan(temp_grid[jX+indeXo*ncellX])){cout << "NaN!" << endl; exit(1);}
//	for my memory, previous version	    temp_grid[jX+indeX1*ncellX]=temp_grid[jX+indeX1*ncellX]+Peff*Vx*sign1*(LY/ncellY*(indeX1+1)-Yo1)/abs(Vx)/DeltaY/my_grid.get_area(jX,indeX1,Par);
	    }
    }
}

void calcul_contribution_MC(double* temp_grid, double Xo, double Yo, double X1, double Y1, double Vx, double Vy, int sign1, double elapsed_time, double Peff, grid my_grid, int ncellX, int ncellY, double LX, double LY, Parameters Par)
{
double DeltaX=LX/ncellX;
double DeltaY=LY/ncellY;
int indeXo=(int)floor(Xo/DeltaX);
int indeX1=(int)floor(X1/DeltaX);
//cout << "Xo " << indeXo << endl;
//cout << "X1 " << indeX1 << endl;
double Xo1;
double Yo1;
double X11;
double Y11;
if (indeX1==ncellX){
    indeX1=ncellX-1;
}
if (indeXo==ncellX){
    indeXo=ncellX-1;
}

if (indeXo==indeX1){

    calcul_contribution_MC_elem(temp_grid, Xo, Yo, X1, Y1, Vx, Vy, sign1, elapsed_time, Peff, my_grid, ncellX, ncellY, LX, LY, indeXo, Par);

}else{
	if(indeX1>indeXo){
        for (int j=indeXo+1; j<indeX1; j++){
			Xo1=j*LX/ncellX;
			X11=(j+1)*LX/ncellX;
			Yo1=Yo+(Y1-Yo)/(X1-Xo)*(Xo1-Xo);
			Y11=Yo+(Y1-Yo)/(X1-Xo)*(X11-Xo);
            calcul_contribution_MC_elem(temp_grid, Xo1, Yo1, X11, Y11, Vx, Vy, sign1, elapsed_time, Peff, my_grid, ncellX, ncellY, LX, LY, j, Par);
        }
		Xo1=(indeXo+1)*LX/ncellX;
		X11=(indeX1)*LX/ncellX;
		Yo1=Yo+(Y1-Yo)/(X1-Xo)*(Xo1-Xo);
		Y11=Yo+(Y1-Yo)/(X1-Xo)*(X11-Xo);
//		cout << "Yo1 " << Yo1 << endl;
//					cout << "Y11 " << Y11 << endl;
		calcul_contribution_MC_elem(temp_grid, Xo, Yo, Xo1, Yo1, Vx, Vy, sign1, elapsed_time, Peff, my_grid, ncellX, ncellY, LX, LY, indeXo, Par);
		calcul_contribution_MC_elem(temp_grid, X11, Y11, X1, Y1,  Vx, Vy, sign1, elapsed_time, Peff, my_grid, ncellX, ncellY, LX, LY, indeX1, Par);

	}
	if(indeX1<indeXo){
        for (int j=indeX1+1; j<indeXo; j++){
			Xo1=j*LX/ncellX;
			X11=(j+1)*LX/ncellX;
			Yo1=Yo+(Y1-Yo)/(X1-Xo)*(Xo1-Xo);
			Y11=Yo+(Y1-Yo)/(X1-Xo)*(X11-Xo);
//			cout << "Yo1 " << Yo1 << " " << X1 << " " << Xo << " " << Xo1 << " " << Y1 << " " << Yo << endl;
//			cout << "Y11 " << Y11 << endl;
            calcul_contribution_MC_elem(temp_grid, Xo1, Yo1, X11, Y11, Vx, Vy, sign1, elapsed_time, Peff, my_grid, ncellX, ncellY, LX, LY, j, Par);
        }
		Xo1=(indeX1+1)*LX/ncellX;
		X11=(indeXo)*LX/ncellX;
		Yo1=Yo+(Y1-Yo)/(X1-Xo)*(Xo1-Xo);
		Y11=Yo+(Y1-Yo)/(X1-Xo)*(X11-Xo);
		calcul_contribution_MC_elem(temp_grid, Xo, Yo, Xo1, Yo1, Vx, Vy, sign1, elapsed_time, Peff, my_grid, ncellX, ncellY, LX, LY, indeX1, Par);
		calcul_contribution_MC_elem(temp_grid, X11, Y11, X1, Y1,  Vx, Vy, sign1, elapsed_time, Peff, my_grid, ncellX, ncellY, LX, LY, indeXo, Par);
		
	}
}
}

void calcul_contribution_MCT(double* temp_grid, double Xo, double Yo, double X1, double Y1, double Vx, double Vy, int sign1, double elapsed_time, double Peff, grid my_grid, double tau_0, Parameters Par, int p)
{
	int ncellX=Par.ncellX;
	int ncellY=Par.ncellY;
	double LX=Par.LX;
	double LY=Par.LY;
	int index0=(int)floor(Xo/LX*ncellX);
	int index1=(int)floor(X1/LX*ncellX);
	int indey0=(int)floor(Yo/LY*ncellY);
	int indey1=(int)floor(Y1/LY*ncellY);
	int incrementx;
	int incrementy;
	int startx;
	int starty;


	double pointsx_vertical[ncellX+2];
	double pointsy_vertical[ncellY+2];
	double pointsx_horiz[ncellX+2];
	double pointsy_horiz[ncellY+2];
	
	if (Xo!=X1 || Yo!=Y1){
	
	if (Xo<X1){
		incrementx=1;
		startx=(int)(floor(Xo/LX*ncellX))+1;
	}
	else
	{
		incrementx=-1;
		if (floor(Xo/LX*ncellX<Xo/LX*ncellX)){
    		startx=(int)(floor(Xo/LX*ncellX));
		}
		else{
			startx=(int)(floor(Xo/LX*ncellX))-1;
		}
	}
	if (Yo<Y1){
		incrementy=1;
		starty=(int)(floor(Yo/LY*ncellY))+1;
	}
	else
	{
		incrementy=-1;
		if (floor(Yo/LY*ncellY<Yo/LY*ncellY)){
		    starty=(int)(floor(Yo/LY*ncellY));
		}
		else
		{
			starty=(int)(floor(Yo/LY*ncellY))-1;
		}
	}
	bool overshoot=0;
	int nx=0;
	int size_v=0;
	while (overshoot==0)
	{
		if (((double)(startx+incrementx*nx)*LX/ncellX-X1)*incrementx<-1e-20)
		{
			 pointsx_vertical[size_v]=((double)(startx+incrementx*nx)*LX/ncellX);
			 pointsy_vertical[size_v]=(Yo+(Y1-Yo)/(X1-Xo)*((startx+incrementx*nx)*LX/ncellX-Xo));
			 if (((double)(startx+incrementx*nx)*LX/ncellX-X1)*incrementx>LX){
			 cout << "check01 " << ((double)(startx+incrementx*nx)*LX/ncellX-X1)*incrementx << endl;
			 }
			 size_v=size_v+1;
		}
		else
		{
			overshoot=1;
		}
		nx=nx+1;
	}
//	cout << "check1"<< endl;
	overshoot=0;
	int ny=0;
	int size_h=0;
	while (overshoot==0)
		{
			if (((double)(starty+incrementy*ny)*LY/ncellY-Y1)*incrementy<-1e-20)
			{
				 pointsx_horiz[size_h]=Xo+(X1-Xo)/(Y1-Yo)*((starty+incrementy*ny)*LY/ncellY-Yo);
				 if ((Xo+(X1-Xo)/(Y1-Yo)*((starty+incrementy*ny)*LY/ncellY-Yo))>LX){
					 cout << "horiz " << endl;
				 }
				 pointsy_horiz[size_h]=((double)(starty+incrementy*ny)*LY/ncellY);
				 size_h=size_h+1;
			}
			else
			{
				overshoot=1;
			}
			ny=ny+1;
		}
	
	double sortedx[size_h+size_v+2];
	double sortedy[size_h+size_v+2];
	sortedx[0]=Xo;
	sortedy[0]=Yo;
//	cout << "check15 " << endl;
	int itv=0;
	int ith=0;
	
	int itv2=size_v;
	
	int ith2=size_h;
	
	int it2=0;
	double d1;
	double d2;
//	cout << "check17" << Xo << " " << Yo << " " << X1 << " " << Y1 << endl;
	if (size_v==0 && size_h>0){
		for (int it=0; it<ith2; it++){
			sortedx[it+1]=pointsx_horiz[it];
			sortedy[it+1]=pointsy_horiz[it];
		}
		sortedx[ith2+1]=X1;
		sortedy[ith2+1]=Y1;
		it2=ith2;
	}
	
	if(size_h==0 && size_v>0){
		for (int it=0; it<itv2; it++){
			sortedx[it+1]=pointsx_vertical[it];
			sortedy[it+1]=pointsy_vertical[it];
		}
		sortedx[itv2+1]=X1;
		sortedy[itv2+1]=Y1;
		it2=itv2;
	}
	if (size_h==0 && size_v==0){
		sortedx[1]=X1;
		sortedy[1]=Y1;
		it2=0;
	}
	if (size_h>0 && size_v>0){
	while (itv<itv2 || ith<ith2)
	{
		if (itv<itv2){
		    d1=(pointsx_vertical[itv]-sortedx[it2])*(pointsx_vertical[itv]-sortedx[it2])+(pointsy_vertical[itv]-sortedy[it2])*(pointsy_vertical[itv]-sortedy[it2]);
		}
		else
		{
			d1=(X1-sortedx[it2])*(X1-sortedx[it2])+(Y1-sortedy[it2])*(Y1-sortedy[it2]);
		}
		
		if (ith<ith2)
		{
			d2=(pointsx_horiz[ith]-sortedx[it2])*(pointsx_horiz[ith]-sortedx[it2])+(pointsy_horiz[ith]-sortedy[it2])*(pointsy_horiz[ith]-sortedy[it2]);
		}
		else
		{
			d2=(X1-sortedx[it2])*(X1-sortedx[it2])+(Y1-sortedy[it2])*(Y1-sortedy[it2]);
		}
		if (d1<d2)
		{
	//		cout << pointsx_vertical[itv] << " " << pointsy_vertical[itv] << endl;
			sortedx[it2+1]=pointsx_vertical[itv];
			sortedy[it2+1]=pointsy_vertical[itv];
			itv=itv+1;
			
		}
		else
		{
		//	cout << *ithx << endl;
	//		cout << pointsx_horiz[itv] << " 2 " << pointsy_horiz[itv] << endl;
			sortedx[it2+1]=pointsx_horiz[ith];
			sortedy[it2+1]=pointsy_horiz[ith];
			ith=ith+1;
		}
		it2=it2+1;
	}
	
	sortedx[it2+1]=X1; //feb 7 2012 just changed from "it2" to "it2+1" in brackets. have to check
	sortedy[it2+1]=Y1; //feb 7 2012 just changed from "it2" to "it2+1" in brackets. have to check
	}
	
	double Xm;
	double Ym;
//	cout << "check2 "<<  endl;

    	for (int i=0; i<it2+1; i++){ //feb 7 2012 just changed from "it2" to "it2+1" in brackets. have to check
	    	Xm=sortedx[i+1]/2+sortedx[i]/2;
		    Ym=sortedy[i+1]/2+sortedy[i]/2;
		    if (Xm>LX || Ym>LX || Ym<0 || Xm <0)
		    {
			    for (int c1=0; c1< itv2; c1++){
				    cout << c1 << " X " << pointsx_vertical[c1] << " Y " << pointsy_vertical[c1] << endl;
			    }
			    for (int c2=0; c2< ith2; c2++){
				    			cout << c2 << " X " << pointsx_horiz[c2] << " Y " << pointsy_horiz[c2] << endl;
					    		cout << " " << (Xo+(X1-Xo)/(Y1-Yo)*((starty+incrementy*c2)*LY/ncellY-Yo)) << " " << (double)(starty+incrementy*c2)*LY/(double)ncellY-Y1 << " " << starty+incrementy*c2 << " " << LY/ncellY << " "<< Y1 << " " <<endl;
						    }
			    for (int j=0; j<it2+1; j++){
				    cout << j << " " << sortedx[j] << " " << sortedy[j] << endl;
				
			    }
			    cout << i << " " << it2+1 << " " << Xo << " " << Yo << " " << X1 << " " << Y1 << " " << Xm << " " << Ym << " " << Vx << " " << Vy << " " << p << endl;
		    }   
//		cout << sortedx[i] << " " << sortedy[i] << " " << sortedx[i+1] << " " << sortedy[i+1] << " " << Xm << " " << Ym << endl;
//		cout << sortedx[i+1] << " " << sortedy[i+1] << " " << Xm << " " << Ym << " " <<(int)floor(Xm/LX*ncellX)+(int)floor(Ym/LY*ncellY)*ncellX << " " << Xo << " " << Vx << endl;
		    temp_grid[(int)floor(Xm/LX*ncellX)+(int)floor(Ym/LY*ncellY)*ncellX]=temp_grid[(int)floor(Xm/LX*ncellX)+(int)floor(Ym/LY*ncellY)*ncellX]+Peff*sign1*sqrt(((sortedx[i+1]-sortedx[i])*(sortedx[i+1]-sortedx[i])+(sortedy[i+1]-sortedy[i])*(sortedy[i+1]-sortedy[i]))/(Vx*Vx+Vy*Vy))*ncellX*ncellY;///my_grid.get_area(jX,indeX1,Par);
    	}
	
//	cout << "check 3"<< endl;
	}
	else{
//		if (Xo==X1 && Yo==Y1){
			temp_grid[(int)floor(Xo/LX*ncellX)+(int)floor(Yo/LY*ncellY)*ncellX]=temp_grid[(int)floor(Xo/LX*ncellX)+(int)floor(Yo/LY*ncellY)*ncellX]+Peff*sign1*elapsed_time*ncellX*ncellY;///my_grid.get_area(jX,indeX1,Par);
//		}
//		else{
	}
	
}

void calcul_contribution_MCH(double* temp_grid, double Xo, double Yo, double X1, double Y1, double Vx, double Vy, int sign1, double elapsed_time, double Peff, grid my_grid, double tau_0, Parameters Par)
{
	int ncellX=Par.ncellX;
	int ncellY=Par.ncellY;
	double LX=Par.LX;
	double LY=Par.LY;
	int index0=(int)floor(Xo/LX*ncellX);
	int index1=(int)floor(X1/LX*ncellX);
	int indey0=(int)floor(Yo/LY*ncellY);
	int indey1=(int)floor(Y1/LY*ncellY);
	int incrementx;
	int incrementy;
	int startx;
	int starty;


	double pointsx_vertical[ncellX+2];
	double pointsy_vertical[ncellY+2];
	double pointsx_horiz[ncellX+2];
	double pointsy_horiz[ncellY+2];
	
	if (Xo!=X1 || Yo!=Y1){
	
	if (Xo<X1){
		incrementx=1;
		startx=(int)(floor(Xo/LX*ncellX))+1;
	}
	else
	{
		incrementx=-1;
		if (floor(Xo/LX*ncellX)<Xo/LX*ncellX){
    		startx=(int)(floor(Xo/LX*ncellX));
		}
		else{
			startx=(int)(floor(Xo/LX*ncellX))-1;
		}
	}
	if (Yo<Y1){
		incrementy=1;
		starty=(int)(floor(Yo/LY*ncellY))+1;
	}
	else
	{
		incrementy=-1;
		if (floor(Yo/LY*ncellY)<Yo/LY*ncellY){
		    starty=(int)(floor(Yo/LY*ncellY));
		}
		else
		{
			starty=(int)(floor(Yo/LY*ncellY))-1;
		}
	}
	bool overshoot=0;
	
	int nx=0;
	int size_v=0;
	while (overshoot==0)
	{
		if (((double)(startx+incrementx*nx)*LX/ncellX-X1)*incrementx<-1e-20)
		{
			 pointsx_vertical[size_v]=((double)(startx+incrementx*nx)*LX/ncellX);
			 pointsy_vertical[size_v]=(Yo+(Y1-Yo)/(X1-Xo)*((startx+incrementx*nx)*LX/ncellX-Xo));
			 if (((double)(startx+incrementx*nx)*LX/ncellX-X1)*incrementx>LX){
			 cout << "check01 " << ((double)(startx+incrementx*nx)*LX/ncellX-X1)*incrementx << endl;
			 }
			 size_v=size_v+1;
		}
		else
		{
			overshoot=1;
		}
		nx=nx+1;
	}
//	cout << "check1"<< endl;
	overshoot=0;
	int ny=0;
	int size_h=0;
	while (overshoot==0)
		{
			if (((double)(starty+incrementy*ny)*LY/ncellY-Y1)*incrementy<-1e-20)
			{
				 pointsx_horiz[size_h]=Xo+(X1-Xo)/(Y1-Yo)*((starty+incrementy*ny)*LY/ncellY-Yo);
				 if ((Xo+(X1-Xo)/(Y1-Yo)*((starty+incrementy*ny)*LY/ncellY-Yo))>LX){
					 cout << "horiz " << endl;
				 }
				 pointsy_horiz[size_h]=((double)(starty+incrementy*ny)*LY/ncellY);
				 size_h=size_h+1;
			}
			else
			{
				overshoot=1;
			}
			ny=ny+1;
		}
	
	double sortedx[size_h+size_v+2];
	double sortedy[size_h+size_v+2];
	sortedx[0]=Xo;
	sortedy[0]=Yo;
//	cout << "check15 " << endl;
	int itv=0;
	int ith=0;
	
	int itv2=size_v;
	
	int ith2=size_h;
	
	int it2=0;
	double d1;
	double d2;
//	cout << "check17" << Xo << " " << Yo << " " << X1 << " " << Y1 << endl;
	if (size_v==0 && size_h>0){
		for (int it=0; it<ith2; it++){
			sortedx[it+1]=pointsx_horiz[it];
			sortedy[it+1]=pointsy_horiz[it];
		}
		sortedx[ith2+1]=X1;
		sortedy[ith2+1]=Y1;
		it2=ith2;
	}
	
	if(size_h==0 && size_v>0){
		for (int it=0; it<itv2; it++){
			sortedx[it+1]=pointsx_vertical[it];
			sortedy[it+1]=pointsy_vertical[it];
		}
		sortedx[itv2+1]=X1;
		sortedy[itv2+1]=Y1;
		it2=itv2;
	}
	if(size_h==0 && size_v==0){
		sortedx[1]=X1;
		sortedy[1]=Y1;
		it2=0;
	}
	if (size_h>0 && size_v>0){
	while (itv<itv2 || ith<ith2)
	{
		if (itv<itv2){
		    d1=(pointsx_vertical[itv]-sortedx[it2])*(pointsx_vertical[itv]-sortedx[it2])+(pointsy_vertical[itv]-sortedy[it2])*(pointsy_vertical[itv]-sortedy[it2]);
		}
		else
		{
			d1=(X1-sortedx[it2])*(X1-sortedx[it2])+(Y1-sortedy[it2])*(Y1-sortedy[it2]);
		}
		
		if (ith<ith2)
		{
			d2=(pointsx_horiz[ith]-sortedx[it2])*(pointsx_horiz[ith]-sortedx[it2])+(pointsy_horiz[ith]-sortedy[it2])*(pointsy_horiz[ith]-sortedy[it2]);
		}
		else
		{
			d2=(X1-sortedx[it2])*(X1-sortedx[it2])+(Y1-sortedy[it2])*(Y1-sortedy[it2]);
		}
		if (d1<d2)
		{
	//		cout << pointsx_vertical[itv] << " " << pointsy_vertical[itv] << endl;
			sortedx[it2+1]=pointsx_vertical[itv];
			sortedy[it2+1]=pointsy_vertical[itv];
			itv=itv+1;
			
		}
		else
		{
		//	cout << *ithx << endl;
	//		cout << pointsx_horiz[itv] << " 2 " << pointsy_horiz[itv] << endl;
			sortedx[it2+1]=pointsx_horiz[ith];
			sortedy[it2+1]=pointsy_horiz[ith];
			ith=ith+1;
		}
		it2=it2+1;
	}
	
	sortedx[it2+1]=X1; //feb 7 2012 just changed from "it2" to "it2+1" in brackets. have to check
	sortedy[it2+1]=Y1; //feb 7 2012 just changed from "it2" to "it2+1" in brackets. have to check
	}
	double Xm;
	double Ym;
//	cout << "check2 "<<  endl;
	for (int i=0; i<it2+1; i++){ //feb 7 2012 just changed from "it2" to "it2+1" in brackets. have to check
		Xm=sortedx[i+1]/2+sortedx[i]/2;
		Ym=sortedy[i+1]/2+sortedy[i]/2;
		if (Xm>LX || Ym>LX || Ym<0 || Xm <0)
		{
			for (int c1=0; c1< itv2; c1++){
				cout << c1 << " X " << pointsx_vertical[c1] << " Y " << pointsy_vertical[c1] << endl;
			}
			for (int c2=0; c2< ith2; c2++){
							cout << c2 << " X " << pointsx_horiz[c2] << " Y " << pointsy_horiz[c2] << endl;
							cout << " " << (Xo+(X1-Xo)/(Y1-Yo)*((starty+incrementy*c2)*LY/ncellY-Yo)) << " " << (double)(starty+incrementy*c2)*LY/(double)ncellY-Y1 << " " << starty+incrementy*c2 << " " << LY/ncellY << " "<< Y1 << " " <<endl;
						}
			for (int j=0; j<it2+1; j++){
				cout << j << " " << sortedx[j] << " " << sortedy[j] << endl;
				
			}
			cout << i << " " << it2+1 << " " << Xo << " " << Yo << " " << X1 << " " << Y1 << " " << Xm << " " << Ym << " " << Vx << " " << Vy <<  endl;
		}
//		cout << sortedx[i] << " " << sortedy[i] << " " << sortedx[i+1] << " " << sortedy[i+1] << " " << Xm << " " << Ym << endl;
//		cout << sortedx[i+1] << " " << sortedy[i+1] << " " << Xm << " " << Ym << " " <<(int)floor(Xm/LX*ncellX)+(int)floor(Ym/LY*ncellY)*ncellX << " " << Xo << " " << Vx << endl;
		temp_grid[(int)floor(Xm/LX*ncellX)+(int)floor(Ym/LY*ncellY)*ncellX]=temp_grid[(int)floor(Xm/LX*ncellX)+(int)floor(Ym/LY*ncellY)*ncellX]+Peff*sign1*(sortedx[i+1]-sortedx[i])*ncellX*ncellY;///my_grid.get_area(jX,indeX1,Par);
		
	}
//	cout << "check 3"<< endl;
	}
	
}


string convertInt(int number)
{
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

int main () {
	string s1 = date(); //for naming the output
	RandomClass r;

	// initialize the random number generator
	time_t seconds;
    seconds = time(NULL);
	long seed=(long)seconds;
	r.initialize(seed);

    // simulation inputs
	
	
	
	// load material data
	
	
	
	// load file of prescribed temperature walls
	
	// load file of reflective walls
	
	// load file of periodic walls
	
	// load file of pores
	
	// load file of detectors
	
	
	
	cout << double(clock()-clo)/double(CLOCKS_PER_SEC) << endl;

	return 0;
	}
	
