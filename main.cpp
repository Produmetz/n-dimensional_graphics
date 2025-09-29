#define _USE_MATH_DEFINES

#include <iostream>
#include <random>
#include <windows.h>
#include <vector>
#include <cmath>
#include "MyGraph.h"
#include <cstdio>
#include <algorithm>

using namespace std;

HWND hwnd;
RECT rect = {0};
bool flagControl = false;
bool flagChangeAngles = false;
bool isNormalProjection = false;

#define for_md(type, name_index, initial_index, final_index, increment) \
    for (type name_index = initial_index; name_index < final_index; name_index += increment)

vector<int> slice(vector<int> v, int a, int len)
{
    vector<int> result;
    for (int i = 0; i < len; i++)
    {
        result.push_back(v[a + i]);
    }
    return result;
}

vector<vector<double>> matricesMultiplication(vector<vector<double>> m1, vector<vector<double>> m2)
{
    int h1 = m1.size(), w1 = m1[0].size(), h2 = m2.size(), w2 = m2[0].size();
    if (w1 == h2)
    {
        vector<vector<double>> result(h1, vector<double>(w2, 0));
        for (int i = 0; i < h1; i++)
        {
            for (int j = 0; j < w2; j++)
            {
                for (int k = 0; k < h2; k++)
                {
                    result[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return result;
    }
    return vector<vector<double>>();
}

double ScalarMultiplication(vector<double> v1, vector<double> v2)
{
    double result = 0;
    int len = max(v1.size(), v2.size());
    v1.resize(len, 0);
    v2.resize(len, 0);
    for (int i = 0; i < len; i++)
    {
        result += v1[i] * v2[i];
    }
    return result;
}

vector<double> operator-(vector<double> v1, vector<double> v2)
{
    int len = max(v1.size(), v2.size());
    v1.resize(len, 0);
    v2.resize(len, 0);
    vector<double> result(len);
    for (int i = 0; i < len; i++)
    {
        result[i] = v1[i] - v2[i];
    }
    return result;
}
vector<double> operator+(vector<double> v1, vector<double> v2)
{
    int len = max(v1.size(), v2.size());
    v1.resize(len, 0);
    v2.resize(len, 0);
    vector<double> result(len);
    for (int i = 0; i < len; i++)
    {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

vector<double> operator*(double a, vector<double> v2)
{
    for (int i = 0; i < v2.size(); i++)
    {
        v2[i] = a * v2[i];
    }
    return v2;
}

vector<double> operator*(double a, vector<int> v2)
{
    vector<double> result(v2.size(), 0);
    for (int i = 0; i < v2.size(); i++)
    {
        result[i] = a * v2[i];
    }
    return result;
}

vector<vector<double>> matrixRotation(int dimension, vector<double> angles)
{
    vector<vector<double>> result(dimension, vector<double>(dimension, 0));
    vector<vector<double>> data(dimension, vector<double>(dimension, 0));
    int n = (dimension * (dimension - 1)) / 2, p = 0;
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            result[i][j] = (i == j) ? 1 : 0;
        }
    }
    for (int i = 0; i < dimension; i++)
    {
        for (int l = i + 1; l < dimension; l++)
        {
            for (int j = 0; j < dimension; j++)
            {
                for (int k = 0; k < dimension; k++)
                {
                    if (j == k)
                    {
                        if ((j == i) || (k == l))
                        {
                            data[j][k] = cos(angles[p]);
                        }
                        else
                        {
                            data[j][k] = 1;
                        }
                    }
                    else if ((j == i) && (k == l))
                    {
                        data[i][l] = -sin(angles[p]);
                    }
                    else if ((j == l) && (k == i))
                    {
                        data[l][i] = sin(angles[p]);
                    }
                    else
                    {
                        data[j][k] = 0;
                    }
                }
            }
            p++;
            result = matricesMultiplication(result, data);
            data = vector<vector<double>>(dimension, vector<double>(dimension, 0));
        }
    }
    return result;
}

const double coef_scale_normal_camera = 0.5;
vector<vector<double>> Camera = {
    {100, 0, 0}, {0, 100, 0}, {0, 0, 100}, {0, 0, 0}, {0, 0, 0}};
vector<double> CameraShift = {0, 0, -5, -5, -5};
vector<double> angles = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double normal_by_normal = ScalarMultiplication({Camera[0][2], Camera[1][2], Camera[2][2], Camera[3][2], Camera[4][2]}, {Camera[0][2], Camera[1][2], Camera[2][2], Camera[3][2], Camera[4][2]});
vector<double> bounds{-1000, 1000, -600, 600};
double sqrt_cam1_by_cam1 = sqrt(ScalarMultiplication({Camera[0][0], Camera[1][0], Camera[2][0], Camera[3][0], Camera[4][0]},{Camera[0][0], Camera[1][0], Camera[2][0], Camera[3][0], Camera[4][0]}));
double sqrt_cam2_by_cam2 =  sqrt(ScalarMultiplication({Camera[0][1], Camera[1][1], Camera[2][1], Camera[3][1], Camera[4][1]}, {Camera[0][1], Camera[1][1], Camera[2][1], Camera[3][1], Camera[4][1]}));

struct Point{
    vector<double> point;
    bool isVisible = true;
};

Point Point2DByVector(Point n_dimen)
{

    Point result;
    result.point.resize(2, 0.0);
    vector<double> cam1 = {Camera[0][0], Camera[1][0], Camera[2][0], Camera[3][0], Camera[4][0]};
    vector<double> cam2 = {Camera[0][1], Camera[1][1], Camera[2][1], Camera[3][1], Camera[4][1]};
    if (isNormalProjection)
    {
        vector<double> normal = {Camera[0][2], Camera[1][2], Camera[2][2], Camera[3][2], Camera[4][2]};
        vector<double> point_ = n_dimen.point - CameraShift;
        double scalar_point_normal = ScalarMultiplication(point_, normal);
        double coef_dist =  normal_by_normal / scalar_point_normal;
        
        if(scalar_point_normal < 0){
            
           /*if(!(bounds[0] <= result.point[0] && result.point[0] <= bounds[1] && bounds[2] <= result.point[1] && result.point[1] <= bounds[3])){
                result.isVisible = false;
           }*/
            coef_dist = -coef_dist;
            result.isVisible = false;
        }
        point_ = coef_dist * point_ - normal;
        result.point[0] = ScalarMultiplication(cam1, point_) / sqrt_cam1_by_cam1;
        result.point[1] = ScalarMultiplication(cam2, point_) / sqrt_cam2_by_cam2;
    }
    else
    {
        result.point[0] = ScalarMultiplication(cam1, n_dimen.point - CameraShift) / sqrt_cam1_by_cam1;
        result.point[1] = ScalarMultiplication(cam2, n_dimen.point - CameraShift) / sqrt_cam2_by_cam2;
    }

    return result;
}

struct Line {
    Point start;
    Point end;
    bool isVisible = true;
};


Line Line2DByVector(Line line)
{
    Line result;
    result.start = Point2DByVector(line.start);
    result.end = Point2DByVector(line.end);
    
    if (!result.start.isVisible && !result.end.isVisible) {
        result.isVisible = false;
        return result;
    }

    return result;
}

vector<vector<double>> RotationCamera(int dimension, vector<double> angles)
{
    return matricesMultiplication(matrixRotation(dimension, angles), Camera);
}

double Module(vector<int> v1, int n)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        result += (v1[i] - v1[i + n]) * (v1[i] - v1[i + n]);
    }
    return sqrt(result);
}

void DrawHypercube(HDC hdc, double H)
{
    int n = 5;
    vector<Point> vertices(1 << n);
    for (int i = 0; i < (1 << n); i++)
    {
        vector<double> vertex(n);
        for (int j = 0; j < n; j++)
        {
            vertex[j] = ((i >> j) & 1) ? H/2 : -H/2;
        }
        vertices[i].point = vertex;
    }

    for (int i = 0; i < (1 << n); i++)
    {
        for (int j = 0; j < n; j++)
        {
            int neighbor = i ^ (1 << j);
            if (i < neighbor)
            {
                Line edge;
                edge.start = vertices[i];
                edge.end = vertices[neighbor];
                Line projectedEdge = Line2DByVector(edge);
                if (projectedEdge.isVisible)
                {
                    DrawLine(hdc, projectedEdge.start.point[0], projectedEdge.start.point[1], projectedEdge.end.point[0], projectedEdge.end.point[1]);
                }
            }
        }
    }
}

vector<double> retVector(vector<vector<double>> point)
{
    vector<double> result(point.size());
    for (int i = 0; i < result.size(); i++)
    {
        result[i] = point[i][0];
    }
    return result;
}

void DrawSphereN(HDC hdc, double R, int numPoints = 1000)
{
    int n = 5;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(0, 1);

    for (int i = 0; i < numPoints; ++i)
    {
        vector<double> point(n);
        double length_sq = 0.0;

        for (int j = 0; j < n; ++j)
        {
            point[j] = dis(gen);
            length_sq += point[j] * point[j];
        }

        double length = sqrt(length_sq);
        for (int j = 0; j < n; ++j)
        {
            point[j] = point[j] * R / length;
        }

        Point p;
        p.point = point;
        Point projected = Point2DByVector(p);
        if (projected.isVisible) {
            DrawPoint(hdc, projected.point[0], projected.point[1]);
        }
    }
}

void DrawArrayVertexes(HDC hdc, vector<vector<double>> Vertexes)
{
    for (size_t i = 0; i < Vertexes.size(); i++)
    {
        Point p;
        p.point = Vertexes[i];
        Point projected = Point2DByVector(p);
        if (projected.isVisible) {
            DrawPoint(hdc, projected.point[0], projected.point[1]);
        }
    }
}

void DrawGrid(HDC hdc, vector<double> startPoint,
              vector<double> direction1, vector<double> direction2,
              int steps1, int steps2)
{
    int dimensions = startPoint.size();

    vector<double> step1 = (1.0 / steps1) * direction1;
    vector<double> step2 = (1.0 / steps2) * direction2;

    for (int i = 0; i <= steps1; ++i)
    {
        vector<double> lineStart = startPoint + (double)i * step1;
        vector<double> lineEnd = lineStart + direction2;

        Line line;
        line.start.point = lineStart;
        line.end.point = lineEnd;
        Line projectedLine = Line2DByVector(line);
        if (projectedLine.isVisible) {
            DrawLine(hdc, projectedLine.start.point[0], projectedLine.start.point[1], projectedLine.end.point[0], projectedLine.end.point[1]);
        }
    }

    for (int j = 0; j <= steps2; ++j)
    {
        vector<double> lineStart = startPoint + (double)j * step2;
        vector<double> lineEnd = lineStart + direction1;

        Line line;
        line.start.point = lineStart;
        line.end.point = lineEnd;
        Line projectedLine = Line2DByVector(line);
        if (projectedLine.isVisible) {
            DrawLine(hdc, projectedLine.start.point[0], projectedLine.start.point[1], projectedLine.end.point[0], projectedLine.end.point[1]);
        }
    }
}

void WinShow(HDC dc)
{
    if (rect.right == 0 && rect.bottom == 0)
        return;

    HDC memDc = CreateCompatibleDC(dc);
    HBITMAP memBM = CreateCompatibleBitmap(dc, rect.right - rect.left, rect.bottom - rect.top);
    SelectObject(memDc, memBM);

    HBRUSH hBrush = CreateSolidBrush(RGB(0, 0, 0));
    FillRect(memDc, &rect, hBrush);
    DeleteObject(hBrush);

    int w = rect.right - rect.left, h = rect.bottom - rect.top;
    SetCenter(w, h);

    if (flagControl)
    {
        if (flagChangeAngles)
        {
            Camera = RotationCamera(5, angles);
            flagChangeAngles = false;
        }
        SetColor(RGB(255, 0, 0));
        Line xAxis;
        xAxis.start.point = {0,0,0,0,0};
        xAxis.end.point = {500,0,0,0,0};
        Line projectedXAxis = Line2DByVector(xAxis);
        if (projectedXAxis.isVisible) {
            DrawLine(memDc, projectedXAxis.start.point[0], projectedXAxis.start.point[1], projectedXAxis.end.point[0], projectedXAxis.end.point[1]);
        }

        Line yAxis;
        yAxis.start.point = {0,0,0,0,0};
        yAxis.end.point = {0,500,0,0,0};
        Line projectedYAxis = Line2DByVector(yAxis);
        if (projectedYAxis.isVisible) {
            DrawLine(memDc, projectedYAxis.start.point[0], projectedYAxis.start.point[1], projectedYAxis.end.point[0], projectedYAxis.end.point[1]);
        }

        Line zAxis;
        zAxis.start.point = {0,0,0,0,0};
        zAxis.end.point = {0,0,500,0,0};
        Line projectedZAxis = Line2DByVector(zAxis);
        if (projectedZAxis.isVisible) {
            DrawLine(memDc, projectedZAxis.start.point[0], projectedZAxis.start.point[1], projectedZAxis.end.point[0], projectedZAxis.end.point[1]);
        }

        Line wAxis;
        wAxis.start.point = {0,0,0,0,0};
        wAxis.end.point = {0,0,0,500,0};
        Line projectedWAxis = Line2DByVector(wAxis);
        if (projectedWAxis.isVisible) {
            DrawLine(memDc, projectedWAxis.start.point[0], projectedWAxis.start.point[1], projectedWAxis.end.point[0], projectedWAxis.end.point[1]);
        }

        Line vAxis;
        vAxis.start.point = {0,0,0,0,0};
        vAxis.end.point = {0,0,0,0,500};
        Line projectedVAxis = Line2DByVector(vAxis);
        if (projectedVAxis.isVisible) {
            DrawLine(memDc, projectedVAxis.start.point[0], projectedVAxis.start.point[1], projectedVAxis.end.point[0], projectedVAxis.end.point[1]);
        }

        SetColor(RGB(0, 250, 40));
        DrawHypercube(memDc, 250);
         //DrawSphereN(memDc, 250);
        //DrawSphereN(memDc, 25, 3000); // Рисуем сферу с 5000 точек
        //  Рисуем сетку между двумя точками
        //vector<double> start = {-100, -100, -100, -100, -100};
        //vector<double> dir1 = {200, 0, 0, 0, 0}; // Направление по первой оси
        //vector<double> dir2 = {0, 200, 0, 0, 0}; // Направление по второй оси

        // SetColor(RGB(255, 255, 0));
        //DrawGrid(memDc, start, dir1, dir2, 5, 5); // 10 частей в каждом направлении
    }
    else
    {
        SetColor(RGB(255, 255, 255));
        const wchar_t *msg = L"Click to enable control";
        TextOutW(memDc, 10, 10, msg, wcslen(msg));
    }

    BitBlt(dc, 0, 0, w, h, memDc, 0, 0, SRCCOPY);
    DeleteDC(memDc);
    DeleteObject(memBM);
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam)
{
    const int h_ = 50;
    switch (message)
    {
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    case WM_KEYDOWN:
        if (flagControl)
        {
            switch (wparam)
            {
            case '0':
            case VK_NUMPAD0:
                angles[0] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '1':
            case VK_NUMPAD1:
                angles[1] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '2':
            case VK_NUMPAD2:
                angles[2] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '3':
            case VK_NUMPAD3:
                angles[3] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '4':
            case VK_NUMPAD4:
                angles[4] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '5':
            case VK_NUMPAD5:
                angles[5] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '6':
            case VK_NUMPAD6:
                angles[6] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '7':
            case VK_NUMPAD7:
                angles[7] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '8':
            case VK_NUMPAD8:
                angles[8] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case '9':
            case VK_NUMPAD9:
                angles[9] += M_PI / h_;
                flagChangeAngles = true;
                break;
            case 'A':
                CameraShift[0] -= 10;
                break;
            case 'S':
                CameraShift[1] -= 10;
                break;
            case 'D':
                CameraShift[2] -= 10;
                break;
            case 'F':
                CameraShift[3] -= 10;
                break;
            case 'G':
                CameraShift[4] -= 10;
                break;
            case 'Q':
                CameraShift[0] += 10;
                break;
            case 'W':
                CameraShift[1] += 10;
                break;
            case 'E':
                CameraShift[2] += 10;
                break;
            case 'R':
                CameraShift[3] += 10;
                break;
            case 'T':
                CameraShift[4] += 10;
                break;

            case 'N': // Добавить этот case
                isNormalProjection = !isNormalProjection;
                InvalidateRect(hwnd, NULL, TRUE);
                break;
            }
            // Принудительная перерисовка после изменений
            InvalidateRect(hwnd, NULL, TRUE);
        }
        break;
    case WM_SIZE:
        GetClientRect(hwnd, &rect);
        break;
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);
        WinShow(hdc);
        EndPaint(hwnd, &ps);
        break;
    }
    case WM_LBUTTONDOWN:
        flagControl = !flagControl;
        InvalidateRect(hwnd, NULL, TRUE);
        break;
    default:
        return DefWindowProc(hwnd, message, wparam, lparam);
    }
    return 0;
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    WNDCLASSW wcl = {0};
    wcl.lpszClassName = L"myWind";
    wcl.lpfnWndProc = WndProc;
    wcl.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wcl.hCursor = LoadCursor(NULL, IDC_ARROW);
    RegisterClassW(&wcl);

    hwnd = CreateWindowW(L"myWind", L"My Program", WS_OVERLAPPEDWINDOW,
                         CW_USEDEFAULT, CW_USEDEFAULT, 800, 600,
                         NULL, NULL, hInstance, NULL);

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return msg.wParam;
}