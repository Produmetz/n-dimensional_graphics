#include "MyGraph.h"
#include <cstdio>

POINT center;
DOUBLE COS_ = cos(M_PI - M_PI / 12), SIN_ = sin(M_PI - M_PI / 12);
DOUBLE COS = cos(M_PI + M_PI / 12), SIN = sin(M_PI + M_PI / 12);
COLORREF color;
HPEN hPen;

void SetColor(COLORREF cl){
    color=cl;
}
void SetCenter(int w,int h){
    center.x=w/2;
    center.y=h/2;
}
void DrawPoint(HDC hdc,int x,int y){
    SetPixel(hdc,center.x+x,center.y-y,color);	
    
}
void DrawLine(HDC hdc,int x1,int y1,int x2,int y2){
    HPEN hOldPen;
    hPen = CreatePen(PS_SOLID, 1, color);
    hOldPen = (HPEN)SelectObject(hdc, hPen);
    
    MoveToEx(hdc, center.x + x1, center.y - y1, NULL);
    LineTo(hdc, center.x + x2, center.y - y2);
    
    SelectObject(hdc, hOldPen);
    DeleteObject(hPen);
}

void DrawSystemCoordinat(HDC hdc){
    int h=center.y*2/15,w=center.x*2/20;
    for(int i=1;i<20;i++){
        DrawLine( hdc,-w*10,i*h,w*10,i*h);
        DrawLine( hdc,-w*10,-i*h,w*10,-i*h);
        DrawLine(hdc,i*w,-h*15/2,i*w,h*15/2);
        DrawLine(hdc,-i*w,-h*15/2,-i*w,h*15/2);
    }
    SetColor(RGB(0,0,255));
    DrawLine( hdc,-w*10,0,w*10,0);
    DrawLine(hdc,0,-h*15/2,0,h*15/2);
}
void DrawCircle(HDC hdc,double Radius,double Center_x,double Center_y){
    SelectObject(hdc,GetStockObject(DC_PEN));
    SetDCPenColor(hdc,color);
    Ellipse(hdc,center.x+Center_x-Radius,center.y-Center_x-Radius,center.x+Center_x+Radius,center.y-Center_x+Radius);
}
void DrawParallelScheme(HDC hdc,double R,double h,int N,int K){
    double R_cos,R_sin,hcos,hsin,angle=M_PI/N,fi;
    for(int n=0;n<N;n++){
        fi=n*angle;
        R_cos=R*cos(fi),R_sin=R*sin(fi),hcos=h*cos(fi),hsin=h*sin(fi);
        for(int k=0;k<K;k++){
            DrawLine(hdc,R_cos-k*hsin,R_sin+k*hcos,-R_cos-k*hsin,-R_sin+k*hcos);
            DrawLine(hdc,R_cos+k*hsin,R_sin-k*hcos,-R_cos+k*hsin,-R_sin-k*hcos);
        }
    }
}
void DrawFanScheme(HDC hdc,double R,double alpha,int N,int K){
    double R_cos,R_sin,angle=2*M_PI/N,fi;
    for(int n=0;n<N;n++){
        fi=n*angle-M_PI;
        R_cos=R*cos(fi+M_PI),R_sin=R*sin(fi+M_PI);
        for(int k=0;k<K;k++){
            DrawLine(hdc,R_cos,R_sin,R*cos(fi+k*alpha),R*sin(fi+k*alpha));
            DrawLine(hdc,R_cos,R_sin,R*cos(fi-k*alpha),R*sin(fi-k*alpha));
        }
    }
}
void Draw3D(HDC hdc,double r1,double r2,double r3){
    double c1=-700,c2=0,c3=0,b1x=0,b1y=10,b1z=0,b2x=0,b2y=0,b2z=10,bx=0,by=0,bz=0,Del,Del1,Del2,Del3,alpha,bet1,bet2;
    Del =(r3-c3)*(b1x*b2y-b1y*b2x)-(r2-c2)*(b1x*b2z-b1z*b2x)+-(r1-c1)*(b1y*b2z-b1z*b2y);
    Del1=(bz-c3)*(b1x*b2y-b1y*b2x)-(by-c2)*(b1x*b2z-b1z*b2x)+-(bx-c1)*(b1y*b2z-b1z*b2y);
    Del2=(r3-c3)*((bx-c1)*b2y-(by-c2)*b2x)-(r2-c2)*((bx-c1)*b2z-(bz-c3)*b2x)+-(r1-c1)*((by-c2)*b2z-(bz-c3)*b2y);
    Del3=(r3-c3)*(b1x*(by-c2)-b1y*(bx-c1))-(r2-c2)*(b1x*(bz-c3)-b1z*(bx-c1))+-(r1-c1)*(b1y*(bz-c3)-b1z*(by-c2));
    
    alpha=Del1/Del;
    bet1=Del2/Del;
    bet2=Del3/Del;
    
    
    SetColor(RGB(255,0,0));
    DrawPoint(hdc,bet1,bet2);
}
void DrawSphere(HDC hdc, double radius){
    for(double fi1 = 0; fi1 <= 2*M_PI; fi1 += 0.25){
        for(double fi2 = 0; fi2 <= M_PI; fi2 += 0.25){
            Draw3D(hdc,radius*sin(fi2)*cos(fi1),radius*sin(fi2)*sin(fi1),radius*cos(fi2));
        }
    }
}
void DrawLineWithArrow(HDC hdc,int x1,int y1,int x2,int y2){
    double radius=sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    double sin_fi=(y2-y1)/radius;
    double cos_fi=(x2-x1)/radius;
    double scale=radius/3;
    DrawLine(hdc,x1,y1,x2,y2);
    DrawLine(hdc,x2,y2,x2+scale*(cos_fi*COS_-sin_fi*SIN_),y2+scale*(sin_fi*COS_+cos_fi*SIN_));
    DrawLine(hdc,x2,y2,x2+scale*(cos_fi*COS-sin_fi*SIN),y2+scale*(sin_fi*COS+cos_fi*SIN));
}
void DrawGrad(HDC hdc){
    for(double a=-400;a<=400;a+=5){
        for(double b=-400;b<=400;b+=5){
            DrawLineWithArrow(hdc,a,b,a+(a*cos(a*b)-a*sin(a*b))/25,b+(b*cos(a*b)-b*sin(a*b))/25);
        }
    }
}
void DrawPoints(HDC hdc,double *array_x,double *array_y,INT CountP){
    for(int i=0;i<CountP;i++){
        DrawPoint(hdc,array_x[i],array_y[i]);
    }
}

