#pragma once

double Pi = 3.14159;

/* Количество элементов */
int max_str, max_node, max_el;

/* Структуры точек и структура элементов*/
struct Point
{

    /* Индивидуальный номер узла */
    int Num_node;

    /*Принадлежность границе*/
    bool Boundary;

    /* Координаты узла */
    double x, y;

};
struct Element
{
    /*Индивидуальный номер элемента */
    int Num_el;

    /* Тип элемента */
    int Geom_el;

    /* Тег границы */
    int Num_bound;

    /* Номера вершин элемента */
    int Num_vert_1, Num_vert_2, Num_vert_3;

    /* Координаты вершин элемента */
    Point Coord_vert_1;
    Point Coord_vert_2;
    Point Coord_vert_3;

    /* Номера граничных вершин элемента */
    int Num_bound_vert_1;
    int Num_bound_vert_2;

    /* Длина грани элемента */
    double Length_face_el_1;
    double Length_face_el_2;
    double Length_face_el_3;

    /* Площадь элемента */
    double Area_el;

    /* Координаты центра элемента */
    Point Coord_center_el;

    /* Расстояние от центра до грани элемента */
    double h_1, h_2, h_3;

    /* Соседний элемент */
    int N1_e, N2_e, N3_e;

    /* Значение температуры в КО на старом слое по времени */
    double t;

    /* Значение температуры в КО на новом слое по времени */
    double T;

    /* Значение ст. отверждения в КО на старом слое по времени */
    double alfa;

    /* Значение ст. отверждения в КО на новом слое по времени */
    double Alfa;

    /* Значение произв. от Alfa */
    double dalfa;

};

/* Вектор точек и вектор элементов */
vector<Point> vectorPoint;
vector<Element> vectorElement;

/* Счетчик итераций */
int Iter_Glob;

/* Граничные условия */
double Tw_0 = 100 + 273.15;
double Tw_1 = 230 + 273.15;
double T0 = (Tw_0 - Tw_0) / (Tw_1 - Tw_0);
double T1 = (Tw_1 - Tw_0) / (Tw_1 - Tw_0);

/* Номер элемента для отслеживания */
int num_el_1, num_el_2, num_el_3;

/* Координаты контрольного элемента */
double xx_1 = 0.5, yy_1 = 0.1;

/* Параметр установления */
double E_T;

/* Директория файла с сеткой и Save*/
string File_Mesh_Name =
"Documents/Mesh/Mesh_Coaxial_Cylinders_(El=2870).msh";
ifstream File_Mesh(File_Mesh_Name);

bool Read_From_Save = false;
string File_Save_Name =
"Documents/Save/Save_El = 2870/Save_(El = 726).DAT";

/* Параметры */
double rho = 2080;
double L = 0.01;
double lambda = 0.432;
double Da;

/* Шаг и счетчик времени */
double dt = 0.00001;
double _time = 0.0;
double final_time = 1.0;
double dem_time = 0.0;

/* Постоянная/переменная теплоемкость */
bool Constant_Heat_Capacity = false;

/* Теплоизоляция */
bool Thermal_Insulation = true;
