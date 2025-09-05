#define _USE_MATH_DEFINES

#include <iostream>
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

vector<vector<double>> Camera = {
    {1, 0}, {0, 1}, {0, 0}, {0, 0}, {0, 0}};
vector<double> CameraShift = {0, 0, 0, 0, 0};
vector<double> angles = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

vector<double> Point2DByVector(vector<double> Point)
{
    vector<double> result(2);
    vector<double> cam1 = {Camera[0][0], Camera[1][0], Camera[2][0], Camera[3][0], Camera[4][0]};
    vector<double> cam2 = {Camera[0][1], Camera[1][1], Camera[2][1], Camera[3][1], Camera[4][1]};
    result[0] = ScalarMultiplication(cam1, Point - CameraShift);
    result[1] = ScalarMultiplication(cam2, Point - CameraShift);
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
    int n = 4;
    vector<double> vertex2d1(2, 0), vertex2d2(2, 0);
    
    // Проходим по всем вершинам гиперкуба
    for (int i = 0; i < (1 << n); i++)
    {
        vector<int> vertex1(n, 0);
        for (int j = 0; j < n; j++)
        {
            vertex1[j] = (i >> j) & 1;
        }
        
        // Рисуем связи между вершинами
        for (int j = 0; j < n; j++)
        {
            vector<int> vertex2 = vertex1;
            vertex2[j] = 1 - vertex2[j]; // Инвертируем j-й бит
            
            vertex2d1 = Point2DByVector(H * vector<double>(vertex1.begin(), vertex1.end()));
            vertex2d2 = Point2DByVector(H * vector<double>(vertex2.begin(), vertex2.end()));
            
            DrawLine(hdc, vertex2d1[0], vertex2d1[1], vertex2d2[0], vertex2d2[1]);
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

vector<vector<double>> DrawSphereN(double H)
{
    int n = 5, k = (n * n - n) / 2;
    vector<vector<double>> result;
    vector<vector<double>> point1 = {{H}, {0}, {0}, {0}, {0}};
    
    // Упрощенная версия без использования проблемного макроса
    for (double angle1 = 0; angle1 < 2 * M_PI; angle1 += 1.25)
    {
        for (double angle2 = 0; angle2 < 2 * M_PI; angle2 += 1.25)
        {
            vector<double> index = {angle1, angle2};
            vector<vector<double>> point2 = matricesMultiplication(matrixRotation(n, index), point1);
            result.push_back(retVector(point2));
        }
    }
    return result;
}

void DrawArrayVertexes(HDC hdc, vector<vector<double>> Vertexes)
{
    vector<double> vertex2d(2, 0);
    for (size_t i = 0; i < Vertexes.size(); i++)
    {
        vertex2d = Point2DByVector(Vertexes[i]);
        DrawPoint(hdc, vertex2d[0], vertex2d[1]);
    }
}

void WinShow(HDC dc)
{
    if (rect.right == 0 && rect.bottom == 0)
        return;

    HDC memDc = CreateCompatibleDC(dc);
    HBITMAP memBM = CreateCompatibleBitmap(dc, rect.right - rect.left, rect.bottom - rect.top);
    SelectObject(memDc, memBM);

    // Заливаем фон черным цветом
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

        vector<double> P2DX1 = Point2DByVector({-700, 0, 0, 0, 0});
        vector<double> P2DX2 = Point2DByVector({700, 0, 0, 0, 0});
        vector<double> P2DY1 = Point2DByVector({0, -700, 0, 0, 0});
        vector<double> P2DY2 = Point2DByVector({0, 700, 0, 0, 0});
        vector<double> P2DZ1 = Point2DByVector({0, 0, -700, 0, 0});
        vector<double> P2DZ2 = Point2DByVector({0, 0, 700, 0, 0});
        vector<double> P2DW1 = Point2DByVector({0, 0, 0, -700, 0});
        vector<double> P2DW2 = Point2DByVector({0, 0, 0, 700, 0});
        vector<double> P2DV1 = Point2DByVector({0, 0, 0, 0, -700});
        vector<double> P2DV2 = Point2DByVector({0, 0, 0, 0, 700});

        SetColor(RGB(255, 0, 0));
        DrawLine(memDc, P2DX1[0], P2DX1[1], P2DX2[0], P2DX2[1]);
        DrawLine(memDc, P2DY1[0], P2DY1[1], P2DY2[0], P2DY2[1]);
        DrawLine(memDc, P2DZ1[0], P2DZ1[1], P2DZ2[0], P2DZ2[1]);
        DrawLine(memDc, P2DW1[0], P2DW1[1], P2DW2[0], P2DW2[1]);
        DrawLine(memDc, P2DV1[0], P2DV1[1], P2DV2[0], P2DV2[1]);
        
        SetColor(RGB(0, 250, 40));
        DrawHypercube(memDc, 250);
    }
    else
    {
        // Рисуем приветственное сообщение
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
                angles[0] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '1':
            case VK_NUMPAD1:
                angles[1] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '2':
            case VK_NUMPAD2:
                angles[2] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '3':
            case VK_NUMPAD3:
                angles[3] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '4':
            case VK_NUMPAD4:
                angles[4] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '5':
            case VK_NUMPAD5:
                angles[5] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '6':
            case VK_NUMPAD6:
                angles[6] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '7':
            case VK_NUMPAD7:
                angles[7] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '8':
            case VK_NUMPAD8:
                angles[8] += M_PI / 12;
                flagChangeAngles = true;
                break;
            case '9':
            case VK_NUMPAD9:
                angles[9] += M_PI / 12;
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