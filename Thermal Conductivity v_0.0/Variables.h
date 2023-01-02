#pragma once

/* Количество элементов */
int max_str;

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

    /* Удаленность от центра */
    double l;

};

/* Вектор точек и вектор элементов */
vector<Point> vectorPoint;
vector<Element> vectorElement;

/* Счетчик итераций */
int Iter_Glob;

/* Массив температуры */
double* T;

/* Граничные условия */
double T0 = 1.0;
double T1 = 0.0;

/* Номер элемента для отслеживания */
int num_el_1, num_el_2, num_el_3;

/* Координаты контрольного элемента */
double xx_1 = 0.5, yy_1 = 0.1;

/* Параметр установления */
double E_T;

/* Директория файла с сеткой */
string File_Mesh_Name =
"Documents/Mesh/Mesh_Coaxial_Cylinders_(El = 178).msh";
ifstream File_Mesh(File_Mesh_Name);

/* Шаг и счетчик времени */
double dt = 0.001;
double _time = 0.0;
double final_time = 1.0;