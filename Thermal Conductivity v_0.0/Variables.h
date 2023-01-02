#pragma once

double Pi = 3.14159;

/* ���������� ��������� */
int max_str, max_node, max_el;

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

    /* �������� ��. ����������� � �� �� ������ ���� �� ������� */
    double alfa;

    /* �������� ��. ����������� � �� �� ����� ���� �� ������� */
    double Alfa;

    /* �������� ������. �� Alfa */
    double dalfa;

};

/* ������ ����� � ������ ��������� */
vector<Point> vectorPoint;
vector<Element> vectorElement;

/* ������� �������� */
int Iter_Glob;

/* ��������� ������� */
double Tw_0 = 100 + 273.15;
double Tw_1 = 230 + 273.15;
double T0 = (Tw_0 - Tw_0) / (Tw_1 - Tw_0);
double T1 = (Tw_1 - Tw_0) / (Tw_1 - Tw_0);

/* ����� �������� ��� ������������ */
int num_el_1, num_el_2, num_el_3;

/* ���������� ������������ �������� */
double xx_1 = 0.5, yy_1 = 0.1;

/* �������� ������������ */
double E_T;

/* ���������� ����� � ������ � Save*/
string File_Mesh_Name =
"Documents/Mesh/Mesh_Coaxial_Cylinders_(El=2870).msh";
ifstream File_Mesh(File_Mesh_Name);

bool Read_From_Save = false;
string File_Save_Name =
"Documents/Save/Save_El = 2870/Save_(El = 726).DAT";

/* ��������� */
double rho = 2080;
double L = 0.01;
double lambda = 0.432;
double Da;

/* ��� � ������� ������� */
double dt = 0.00001;
double _time = 0.0;
double final_time = 1.0;
double dem_time = 0.0;

/* ����������/���������� ������������ */
bool Constant_Heat_Capacity = false;

/* ������������� */
bool Thermal_Insulation = true;
