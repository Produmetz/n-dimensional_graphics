#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <windows.h>
#include <cmath>

extern POINT center;
extern DOUBLE COS_;
extern DOUBLE SIN_;
extern DOUBLE COS;
extern DOUBLE SIN;
extern COLORREF color;
extern HPEN hPen;

void SetColor(COLORREF cl);
void SetCenter(int w, int h);
void DrawPoint(HDC hdc, int x, int y);
void DrawLine(HDC hdc, int x1, int y1, int x2, int y2);
void DrawSystemCoordinat(HDC hdc);
void DrawCircle(HDC hdc, double Radius, double Center_x, double Center_y);
void DrawParallelScheme(HDC hdc, double R, double h, int N, int K);
void DrawFanScheme(HDC hdc, double R, double alpha, int N, int K);
void Draw3D(HDC hdc, double r1, double r2, double r3);
void DrawSphere(HDC hdc, double radius);
void DrawLineWithArrow(HDC hdc, int x1, int y1, int x2, int y2);
void DrawGrad(HDC hdc);
void DrawAproxLine(HDC hdc, double k, double b, int w);
void DrawAproxLine2(HDC hdc, double A, double B, double C, int w);
void DrawPoints(HDC hdc, double* array_x, double* array_y, INT CountP);