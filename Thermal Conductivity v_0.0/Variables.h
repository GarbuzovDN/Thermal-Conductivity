#pragma once

/* ���������� ��������� */
int max_str;

/* ��������� ����� � ��������� ���������*/
struct Point
{

    /* �������������� ����� ���� */
    int Num_node;

    /*�������������� �������*/
    bool Boundary;

    /* ���������� ���� */
    double x, y;

};
struct Element
{
    /*�������������� ����� �������� */
    int Num_el;

    /* ��� �������� */
    int Geom_el;

    /* ��� ������� */
    int Num_bound;

    /* ������ ������ �������� */
    int Num_vert_1, Num_vert_2, Num_vert_3;

    /* ���������� ������ �������� */
    Point Coord_vert_1;
    Point Coord_vert_2;
    Point Coord_vert_3;

    /* ������ ��������� ������ �������� */
    int Num_bound_vert_1;
    int Num_bound_vert_2;

    /* ����� ����� �������� */
    double Length_face_el_1;
    double Length_face_el_2;
    double Length_face_el_3;

    /* ������� �������� */
    double Area_el;

    /* ���������� ������ �������� */
    Point Coord_center_el;

    /* ���������� �� ������ �� ����� �������� */
    double h_1, h_2, h_3;

    /* �������� ������� */
    int N1_e, N2_e, N3_e;

    /* �������� ����������� � �� �� ������ ���� �� ������� */
    double t;

    /* �������� ����������� � �� �� ����� ���� �� ������� */
    double T;

    /* ����������� �� ������ */
    double l;

};

/* ������ ����� � ������ ��������� */
vector<Point> vectorPoint;
vector<Element> vectorElement;

/* ������� �������� */
int Iter_Glob;

/* ������ ����������� */
double* T;

/* ��������� ������� */
double T0 = 1.0;
double T1 = 0.0;

/* ����� �������� ��� ������������ */
int num_el_1, num_el_2, num_el_3;

/* ���������� ������������ �������� */
double xx_1 = 0.5, yy_1 = 0.1;

/* �������� ������������ */
double E_T;

/* ���������� ����� � ������ */
string File_Mesh_Name =
"Documents/Mesh/Mesh_Coaxial_Cylinders_(El = 178).msh";
ifstream File_Mesh(File_Mesh_Name);

/* ��� � ������� ������� */
double dt = 0.001;
double _time = 0.0;
double final_time = 1.0;