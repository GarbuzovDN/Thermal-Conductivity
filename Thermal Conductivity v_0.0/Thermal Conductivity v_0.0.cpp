#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <windows.h>

using namespace std;

#include "Variables.h"

void Arrays_Creation()
{

}

void Arrays_Remove()
{

}

void Find_String(string Str)
{

    if (!File_Mesh)
    {
        cout << "ERROR opening the file: " << File_Mesh_Name << endl << endl << endl;

        exit(1);
    }
    else
    {
        if (max_str != 0) cout << "File opened SUCCESSFULLY: " << File_Mesh_Name << endl;
        if (max_str != 0)
        {

            cout << fixed << setprecision(4) << "Time: " << _time << endl << "Mesh (Number of nodes): " << max_el << endl;

        }

    }

    int num_str;

    string line;

    num_str = 0;

    while (line != Str)
    {

        getline(File_Mesh, line);

        num_str++;

        int pp = 0;

    }

    File_Mesh >> max_str;

}

void Mesh_Init()
{

    Find_String("$Nodes");

    max_node = max_str;

    Point node;

    /* Заполнение структуры точек */
    for (int i = 1; i <= max_node; i++)
    {

        double temp;
        int temp_int;

        File_Mesh >> temp_int >> node.x >> node.y >> temp;

        node.Num_node = temp_int - 1;

        node.Boundary = false;

        vectorPoint.push_back(node);

        i = i;

    }

    Find_String("$Elements");

    max_el = max_str;

    Element el;

    /* Заполнение структуры элементов */
    for (int i = 1; i <= max_el + 1; i++)
    {

        double temp;
        int temp_int;

        File_Mesh >> temp_int >> el.Geom_el;

        el.Num_el = temp_int - 1;

        if (el.Geom_el == 15)
        {
            int temp_int_1;

            File_Mesh >> temp >> temp >> temp >> temp_int_1;
            el.Num_vert_1 = temp_int_1 - 1;
            el.Coord_vert_1 = vectorPoint[el.Num_vert_1];
            vectorPoint[el.Num_vert_1].Boundary = true;

        }
        if (el.Geom_el == 1)
        {

            int temp_int_1, temp_int_2;

            File_Mesh >> temp >> temp >> el.Num_bound >> temp_int_1 >> temp_int_2;
            el.Num_vert_1 = temp_int_1 - 1;
            el.Num_vert_2 = temp_int_2 - 1;

            el.Coord_vert_1 = vectorPoint[el.Num_vert_1];
            el.Coord_vert_2 = vectorPoint[el.Num_vert_2];
            vectorPoint[el.Num_vert_1].Boundary = true;
            vectorPoint[el.Num_vert_2].Boundary = true;

            el.Num_bound_vert_1 = el.Num_vert_1;
            el.Num_bound_vert_2 = el.Num_vert_2;

            i = i;

        }
        if (el.Geom_el == 2)
        {

            int temp_int_1, temp_int_2, temp_int_3;

            File_Mesh >> temp >> temp >> temp >> temp_int_1 >> temp_int_2 >> temp_int_3;
            el.Num_vert_1 = temp_int_1 - 1;
            el.Num_vert_2 = temp_int_2 - 1;
            el.Num_vert_3 = temp_int_3 - 1;

            el.Coord_vert_1 = vectorPoint[el.Num_vert_1];
            el.Coord_vert_2 = vectorPoint[el.Num_vert_2];
            el.Coord_vert_3 = vectorPoint[el.Num_vert_3];

            el.Num_bound = 0;

            double pp = 0.0;

        }

        el.Length_face_el_1 = sqrt(pow((el.Coord_vert_2.x - el.Coord_vert_1.x), 2) + pow((el.Coord_vert_2.y - el.Coord_vert_1.y), 2));
        el.Length_face_el_2 = sqrt(pow((el.Coord_vert_3.x - el.Coord_vert_2.x), 2) + pow((el.Coord_vert_3.y - el.Coord_vert_2.y), 2));
        el.Length_face_el_3 = sqrt(pow((el.Coord_vert_3.x - el.Coord_vert_1.x), 2) + pow((el.Coord_vert_3.y - el.Coord_vert_1.y), 2));

        double a = el.Length_face_el_1;
        double b = el.Length_face_el_2;
        double c = el.Length_face_el_3;
        double p = 0.5 * (a + b + c);
        el.Area_el = sqrt(p * (p - a) * (p - b) * (p - c));

        /* Блок нахождения центра элемента */
        {

            /* Уравнение серединного перпендикуляра 1 */
            double AB = el.Length_face_el_1, BC = el.Length_face_el_2, AC = el.Length_face_el_3;
            double x_a = el.Coord_vert_1.x, x_b = el.Coord_vert_2.x, x_c = el.Coord_vert_3.x;
            double y_a = el.Coord_vert_1.y, y_b = el.Coord_vert_2.y, y_c = el.Coord_vert_3.y;

            el.Coord_center_el.x = (BC * x_a + AC * x_b + AB * x_c) / (AB + BC + AC);
            el.Coord_center_el.y = (BC * y_a + AC * y_b + AB * y_c) / (AB + BC + AC);

        }

        el.h_1 = el.Area_el / p;
        el.h_2 = el.Area_el / p;
        el.h_3 = el.Area_el / p;

        vectorElement.push_back(el);

        i = i;

    }

    int p1, p2;

    /* Блок нахождения соседних и граничных элементов */
    {

        for (int i = 1; i < vectorElement.size() - 1; i++)
        {

            if (vectorElement[i].Geom_el == 2)
            {

                /* Cосед между точками 1 и 2 */
                p1 = vectorElement[i].Num_vert_1;
                p2 = vectorElement[i].Num_vert_2;
                if (vectorPoint[p1].Boundary && vectorPoint[p2].Boundary)
                {

                    vectorElement[i].N1_e = -1;

                    for (int j = 0; j < vectorElement.size(); j++)
                    {

                        if ((p1 == vectorElement[j].Num_bound_vert_1 && p2 == vectorElement[j].Num_bound_vert_2) || (p1 == vectorElement[j].Num_bound_vert_2 && p2 == vectorElement[j].Num_bound_vert_1))
                        {

                            vectorElement[i].Num_bound = vectorElement[j].Num_bound;

                            break;

                        }

                    }

                }
                else
                {

                    for (int j = 0; j < vectorElement.size(); j++)
                    {
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_2 || p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_1) && j != vectorElement[i].Num_el) { vectorElement[i].N1_e = j; break; }
                        if ((p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_2) && j != vectorElement[i].Num_el) { vectorElement[i].N1_e = j; break; }
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_1) && j != vectorElement[i].Num_el) { vectorElement[i].N1_e = j; break; }
                    }

                }

                /* Cосед между точками 2 и 3 */
                p1 = vectorElement[i].Num_vert_2;
                p2 = vectorElement[i].Num_vert_3;
                if (vectorPoint[p1].Boundary && vectorPoint[p2].Boundary)
                {

                    vectorElement[i].N2_e = -1;

                    for (int j = 0; j < vectorElement.size(); j++)
                    {

                        if ((p1 == vectorElement[j].Num_bound_vert_1 && p2 == vectorElement[j].Num_bound_vert_2) || (p1 == vectorElement[j].Num_bound_vert_2 && p2 == vectorElement[j].Num_bound_vert_1))
                        {

                            vectorElement[i].Num_bound = vectorElement[j].Num_bound;

                            break;

                        }

                    }


                }
                else
                {

                    for (int j = 0; j < vectorElement.size(); j++)
                    {
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_2 || p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_1) && j != vectorElement[i].Num_el) { vectorElement[i].N2_e = j; break; }
                        if ((p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_2) && j != vectorElement[i].Num_el) { vectorElement[i].N2_e = j; break; }
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_1) && j != vectorElement[i].Num_el) { vectorElement[i].N2_e = j; break; }
                    }

                }

                /* Cосед между точками 1 и 3 */
                p1 = vectorElement[i].Num_vert_1;
                p2 = vectorElement[i].Num_vert_3;
                if (vectorPoint[p1].Boundary && vectorPoint[p2].Boundary)
                {

                    vectorElement[i].N3_e = -1;

                    for (int j = 0; j < vectorElement.size(); j++)
                    {

                        if ((p1 == vectorElement[j].Num_bound_vert_1 && p2 == vectorElement[j].Num_bound_vert_2) || (p1 == vectorElement[j].Num_bound_vert_2 && p2 == vectorElement[j].Num_bound_vert_1))
                        {

                            vectorElement[i].Num_bound = vectorElement[j].Num_bound;

                            break;

                        }

                    }

                }
                else
                {

                    for (int j = 0; j < vectorElement.size(); j++)
                    {
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_2 || p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_1) && j != vectorElement[i].Num_el) { vectorElement[i].N3_e = j; break; }
                        if ((p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_2) && j != vectorElement[i].Num_el) { vectorElement[i].N3_e = j; break; }
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_1) && j != vectorElement[i].Num_el) { vectorElement[i].N3_e = j; break; }

                        j = j;

                    }

                }

                //cout << vectorElement[i].Num_el << " \t " << vectorElement[i].N1_e << " \t " << vectorElement[i].N2_e << " \t " << vectorElement[i].N3_e << " \t " << vectorElement[i].Num_bound << endl;

                i = i;

            }

        }

    }

    File_Mesh.close();

}

void Blank()
{

    int Local_count = 0;
    int Local_it = 0;

    ofstream File_Blank("Documents/Blank/Blank_R0.BLN", ios_base::trunc);

    /* Запись заголовка бланкировочного файла */
    for (int i = 1; i <= max_el; i++)
    {

        if (vectorElement[i].Geom_el == 1 && vectorElement[i].Num_bound == 2)
        {

            Local_it++;

            if (Local_it == 1) Local_count = i;

            i = i;

        }

    }

    File_Blank << Local_it + 1 << "\t" << "1" << endl;

    /* Запись коррдинат */
    for (int i = 0; i <= max_el; i++)
    {

        if (vectorElement[i].Geom_el == 1 && vectorElement[i].Num_bound == 2)
        {

            File_Blank << fixed << setprecision(4) << vectorElement[i].Coord_vert_1.x << " \t " << vectorElement[i].Coord_vert_1.y << endl;


        }

    }

    File_Blank << fixed << setprecision(4) << vectorElement[Local_count].Coord_vert_1.x << " \t " << vectorElement[Local_count].Coord_vert_1.y;

}

void Redistricting()
{

    double p = 0.0;

    for (int i = 0; i <= max_el; i++)
    {

        vectorElement[i].t = vectorElement[i].T;

        vectorElement[i].alfa = vectorElement[i].Alfa;

    }

}

void Initial_Conditions()
{

    Blank();

    for (int i = 0; i <= max_el; i++)
    {

        vectorElement[i].t = 0.0;
        vectorElement[i].T = 0.0;

        vectorElement[i].alfa = 0.0;
        vectorElement[i].Alfa = 0.0;

    }

}

double C(double Teta)
{
    double T_dim = (Tw_0 - 273.15) + Teta * (Tw_1 - Tw_0);

    return 297.74 + 6.95 * T_dim - 0.0145 * T_dim * T_dim;

}

void Calculation_Temperature()
{
    /* Условие выбора теплоемкости */
    if (Constant_Heat_Capacity == true)
    {
        for (int i = 1; i < max_el; i++)
        {

            /* Первое граничное условие */
            if (vectorElement[i].Num_bound == 1 && vectorElement[i].Geom_el == 2)
            {

                if (vectorElement[i].N1_e == -1)
                {

                    double dT1 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / vectorElement[i].h_1;
                    double dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);
                    double dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                    vectorElement[i].T = dt / (vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                    i = i;

                }

                if (vectorElement[i].N2_e == -1)
                {

                    double dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                    double dT2 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / vectorElement[i].h_2;
                    double dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                    vectorElement[i].T = dt / (vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                    i = i;

                }

                if (vectorElement[i].N3_e == -1)
                {

                    double dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                    double dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);
                    double dT3 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / vectorElement[i].h_3;

                    vectorElement[i].T = dt / (vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                    if (vectorElement[i].T < 0)
                    {
                        i = i;
                    }

                }

            }

            /* Второе граничное условие */
            if (vectorElement[i].Num_bound == 2 && vectorElement[i].Geom_el == 2)
            {

                if (vectorElement[i].N1_e == -1)
                {

                    double dT1 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / vectorElement[i].h_1;
                    double dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);
                    double dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                    vectorElement[i].T = dt / (vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                    i = i;

                }

                if (vectorElement[i].N2_e == -1)
                {

                    double dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                    double dT2 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / vectorElement[i].h_2;
                    double dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                    vectorElement[i].T = dt / (vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                    i = i;

                }

                if (vectorElement[i].N3_e == -1)
                {

                    double dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                    double dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);
                    double dT3 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / vectorElement[i].h_3;

                    vectorElement[i].T = dt / (vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                    if (vectorElement[i].T < 0)
                    {
                        i = i;
                    }

                }

            }

            /* Вычисление температуры */
            if (vectorElement[i].N1_e != -1 && vectorElement[i].N2_e != -1 && vectorElement[i].N3_e != -1 && vectorElement[i].Geom_el == 2)
            {

                double dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                double dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);
                double dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                vectorElement[i].T = dt / (vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                if (dem_time > 300)
                {
                    i = i;
                }

            }

        }

        double q = 0.0;


    }

    if (Constant_Heat_Capacity == false)
    {
        double dT1 = 0, dT2 = 0, dT3 = 0;

        /* Параметры */
        double n = 0.988;
        double R = 8.314;
        double A_1 = 6.86 * 1e+06;
        double A_2 = 5.06 * 1e+06;
        double E_1 = 9.29 * 1e+04;
        double E_2 = 8.52 * 1e+04;

        double a_1 = C(T1) * rho * L * L * A_1 / lambda;
        double a_2 = C(T1) * rho * L * L * A_2 / lambda;
        double E_1_1 = R * Tw_0 / E_1;
        double E_1_2 = R * (Tw_1 - Tw_0) / E_1;
        double E_2_1 = R * Tw_0 / E_2;
        double E_2_2 = R * (Tw_1 - Tw_0) / E_2;
        double C_scale = C(T1);

        for (int i = 0; i < max_el; i++)
        {
            if (vectorElement[i].Geom_el == 2)
            {
                /* Вычисление потоков */
                {
                    if (vectorElement[i].N1_e == -1 || vectorElement[i].N2_e == -1 || vectorElement[i].N3_e == -1)
                    {
                        if (vectorElement[i].N1_e == -1)
                        {
                            if (vectorElement[i].Num_bound == 1)
                            {

                                dT1 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / vectorElement[i].h_1;

                            }

                            else if (vectorElement[i].Num_bound == 2)
                            {

                                dT1 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / vectorElement[i].h_1;

                            }

                            dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);
                            dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                        }

                        if (vectorElement[i].N2_e == -1)
                        {
                            if (vectorElement[i].Num_bound == 1)
                            {

                                dT2 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / vectorElement[i].h_2;

                            }

                            else if (vectorElement[i].Num_bound == 2)
                            {

                                dT2 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / vectorElement[i].h_2;

                            }

                            dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                            dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                        }

                        if (vectorElement[i].N3_e == -1)
                        {
                            if (vectorElement[i].Num_bound == 1)
                            {

                                dT3 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / vectorElement[i].h_3;

                            }

                            else if (vectorElement[i].Num_bound == 2)
                            {

                                dT3 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / vectorElement[i].h_3;

                            }

                            dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                            dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);

                        }
                    }
                    else
                    {

                        dT1 = (vectorElement[vectorElement[i].N1_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[vectorElement[i].N1_e].h_1);
                        dT2 = (vectorElement[vectorElement[i].N2_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[vectorElement[i].N2_e].h_2);
                        dT3 = (vectorElement[vectorElement[i].N3_e].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[vectorElement[i].N3_e].h_3);

                    }
                }

                if (_time >= 1.4719999999)
                {

                    i = i;

                }

                double CC = C(vectorElement[i].t) / C_scale;

                vectorElement[i].T = dt / (CC * vectorElement[i].Area_el) * (dT1 + dT2 + dT3) + vectorElement[i].t;

                double dalfa = pow((1 - vectorElement[i].alfa), n) * (a_1 * exp(-1.0 / (E_1_1 + E_1_2 * vectorElement[i].T)) + a_2 * pow(vectorElement[i].alfa, 0.892) * exp(-1.0 / (E_2_1 + E_2_2 * vectorElement[i].T)));
                vectorElement[i].Alfa = vectorElement[i].alfa + dt * dalfa;



            }
        }

    }
}

void Development()
{

    E_T = 0.0;

    for (int i = 0; i < vectorElement.size(); i++)
    {

        if (vectorElement[i].Geom_el == 2 && E_T < fabs(vectorElement[i].T - vectorElement[i].t) / vectorElement[i].T)
        {

            E_T = fabs(vectorElement[i].T - vectorElement[i].t) / vectorElement[i].T;

        }
    }

}

double Temp_prof(double xx, double yy, string param)
{

    double dl;
    double temp_1 = 99, temp_2 = 100, temp_3 = 101;

    if (param == "NULL")
    {

        for (int i = 0; i < vectorElement.size(); i++)
        {
            if (vectorElement[i].Geom_el == 2)
            {

                dl = sqrt(pow((vectorElement[i].Coord_center_el.x - xx), 2) + pow((vectorElement[i].Coord_center_el.y - yy), 2));

                if (dl < temp_1)
                {

                    temp_1 = dl;
                    num_el_1 = vectorElement[i].Num_el;

                    i = i;

                }
            }

        }

    }
    else
    {

        for (int i = 0; i < vectorElement.size(); i++)
        {
            if (vectorElement[i].Geom_el == 2)
            {

                dl = sqrt(pow((vectorElement[i].Coord_center_el.x - xx), 2) + pow((vectorElement[i].Coord_center_el.y - yy), 2));

                if (dl < temp_1)
                {

                    temp_1 = dl;
                    num_el_1 = vectorElement[i].Num_el;

                    i = i;

                }
            }

        }

        for (int i = 0; i < vectorElement.size(); i++)
        {

            if (vectorElement[i].Geom_el == 2)
            {

                dl = sqrt(pow((vectorElement[i].Coord_center_el.x - xx), 2) + pow((vectorElement[i].Coord_center_el.y - yy), 2));

                if (dl < temp_2 && vectorElement[i].Num_el != num_el_1)
                {

                    temp_2 = dl;

                    num_el_2 = vectorElement[i].Num_el;

                    i = i;

                }

            }

        }

        for (int i = 0; i < vectorElement.size(); i++)
        {

            if (vectorElement[i].Geom_el == 2)
            {

                dl = sqrt(pow((vectorElement[i].Coord_center_el.x - xx), 2) + pow((vectorElement[i].Coord_center_el.y - yy), 2));

                if (dl < temp_3 && vectorElement[i].Num_el != num_el_1 && vectorElement[i].Num_el != num_el_2)
                {

                    temp_3 = dl;
                    num_el_3 = vectorElement[i].Num_el;

                    i = i;

                }

            }

        }

        double a1, a2, b1, b2, c1, c2;
        double Ax = vectorElement[num_el_1].Coord_center_el.x;
        double Ay = vectorElement[num_el_1].Coord_center_el.y;
        double Ex = xx;
        double Ey = yy;
        double Bx = vectorElement[num_el_2].Coord_center_el.x;
        double By = vectorElement[num_el_2].Coord_center_el.y;
        double Cx = vectorElement[num_el_3].Coord_center_el.x;
        double Cy = vectorElement[num_el_3].Coord_center_el.y;

        a1 = Ay - Ey;
        b1 = Ex - Ax;
        c1 = Ax * Ey - Ex * Ay;
        a2 = By - Cy;
        b2 = Cx - Bx;
        c2 = Bx * Cy - Cx * By;

        double det = a1 * b2 - a2 * b1;
        double Dx = (b1 * c2 - b2 * c1) / det;
        double Dy = (a2 * c1 - a1 * c2) / det;

        double BD = sqrt(pow((Bx - Dx), 2) + pow((By - Dy), 2));
        double BC = sqrt(pow((Bx - Cx), 2) + pow((By - Cy), 2));
        double T_D;

        if (param == "Temp") T_D = vectorElement[num_el_2].T + BD / BC * (vectorElement[num_el_3].T - vectorElement[num_el_2].T);
        if (param == "Alfa") T_D = vectorElement[num_el_2].Alfa + BD / BC * (vectorElement[num_el_3].Alfa - vectorElement[num_el_2].Alfa);

        double AE = sqrt(pow((Ax - xx), 2) + pow((Ay - yy), 2));
        double AD = sqrt(pow((Ax - Dx), 2) + pow((Ay - Dy), 2));

        if (param == "Temp") return vectorElement[num_el_1].T + AE / AD * (T_D - vectorElement[num_el_1].T);
        if (param == "Alfa") return vectorElement[num_el_1].Alfa + AE / AD * (T_D - vectorElement[num_el_1].Alfa);

        double pp = 0;

    }

}

void Write_Save()
{
    if (Iter_Glob == 1 || Iter_Glob % 500 == 0)
    {

        string _path = "Documents/Save/Save_El = " + to_string(max_el);
        CreateDirectoryA(_path.c_str(), NULL);

        ofstream File_Save(_path + "/Save_(El = " + to_string(max_el) + ").DAT", ios_base::trunc);

        File_Save << _time << "\t" << max_el << endl;

        File_Save << Pi << "\t" << Tw_0 << "\t" << Tw_1 << "\t" << T0 << "\t" << T1 << "\t" << num_el_1 << "\t" << num_el_2 << "\t" << num_el_3 << endl;

        File_Save << rho << "\t" << L << "\t" << lambda << "\t" << dt << "\t" << _time << "\t" << final_time << "\t" << dem_time << endl;

        File_Save << xx_1 << "\t" << yy_1 << "\t" << E_T << "\t" << File_Mesh_Name << endl;

        /* Запись структуры точек */
        for (int i = 0; i < max_node; i++)
        {

            File_Save << vectorPoint[i].Num_node << "\t" << vectorPoint[i].Boundary << "\t" << vectorPoint[i].x << "\t" << vectorPoint[i].y << endl;

        }

        /* Запись струтктуры элементов */
        for (int i = 0; i < max_el; i++)
        {

            File_Save << vectorElement[i].Num_el << "\t" << vectorElement[i].Geom_el << "\t" << vectorElement[i].Num_vert_1 << "\t" << vectorElement[i].Num_vert_2 << "\t" << vectorElement[i].Num_vert_3 << endl;
            File_Save << vectorElement[i].Coord_vert_1.x << "\t" << vectorElement[i].Coord_vert_1.y << "\t" << vectorElement[i].Coord_vert_2.x << "\t" << vectorElement[i].Coord_vert_2.y << "\t" << vectorElement[i].Coord_vert_3.x << "\t" << vectorElement[i].Coord_vert_3.y << endl;
            File_Save << vectorElement[i].Num_bound_vert_1 << "\t" << vectorElement[i].Num_bound_vert_2 << "\t" << vectorElement[i].Length_face_el_1 << "\t" << vectorElement[i].Length_face_el_2 << "\t" << vectorElement[i].Length_face_el_3 << "\t" << vectorElement[i].Area_el << endl;
            File_Save << vectorElement[i].Coord_center_el.x << "\t" << vectorElement[i].Coord_center_el.y << "\t" << vectorElement[i].h_1 << "\t" << vectorElement[i].h_2 << "\t" << vectorElement[i].h_3 << endl;
            File_Save << vectorElement[i].N1_e << "\t" << vectorElement[i].N2_e << "\t" << vectorElement[i].N3_e << "\t" << vectorElement[i].t << "\t" << vectorElement[i].T << "\t" << vectorElement[i].alfa << "\t" << vectorElement[i].Alfa << "\t" << vectorElement[i].l << endl;

        }

        File_Save.close();

    }

}

void Write_End()
{

    string _path = "Documents/Figure/El = " + to_string(max_el);
    CreateDirectoryA(_path.c_str(), NULL);

    ofstream Field_T(_path + "/1. Field_T_(El = " + to_string(max_el) + ").DAT");
    ofstream Field_Alfa(_path + "/1. Field_Alfa_(El = " + to_string(max_el) + ").DAT");
    ofstream Profile_T(_path + "/2. Profile_T_(El = " + to_string(max_el) + ").DAT");
    ofstream Profile_Alfa(_path + "/2. Profile_Alfa_(El = " + to_string(max_el) + ").DAT");
    ofstream Order_Acc(_path + "/5. Order_Acc_(El = " + to_string(max_el) + ").DAT");

    Field_T << fixed << setprecision(4) << "Time: " << _time << "\t" << "Mesh (Number of cells): " << max_el << endl;
    Field_Alfa << fixed << setprecision(4) << "Time: " << _time << "\t" << "Mesh (Number of cells): " << max_el << endl;

    double Int_Temp = 0.0, h_average = 0.0;

    /* Запись распределния поля T и Alfa */
    for (int i = 0; i <= max_el; i++)
    {

        if (vectorElement[i].Geom_el == 2)
        {

            Field_T << fixed << setprecision(10) << vectorElement[i].Coord_center_el.x << " \t " << vectorElement[i].Coord_center_el.y << " \t "
                << vectorElement[i].T << " \t " << vectorElement[i].t << " \t " << vectorElement[i].Num_el << endl;
            Field_Alfa << fixed << setprecision(10) << vectorElement[i].Coord_center_el.x << " \t " << vectorElement[i].Coord_center_el.y << " \t "
                << vectorElement[i].Alfa << " \t " << vectorElement[i].alfa << " \t " << vectorElement[i].Num_el << endl;

            Int_Temp = Int_Temp + (vectorElement[i].T * vectorElement[i].Area_el);
            h_average = h_average + (vectorElement[i].Length_face_el_1 + vectorElement[i].Length_face_el_2 + vectorElement[i].Length_face_el_3);

        }

    }

    Order_Acc << fixed << setprecision(5) << "Time: " << _time << "\t" << "Mesh (Number of cells): " << max_el << endl << endl;
    Order_Acc << fixed << setprecision(8) << "Average h" << "\t" << h_average / (3 * max_el) << endl;
    Order_Acc << fixed << setprecision(8) << "Int. Temp" << "\t" << Int_Temp << endl;

    Profile_T << fixed << setprecision(4) << "Time: " << _time << "\t" << "Mesh (Number of cells): " << max_el << endl;

    /* Запись значения T и Alfa в сечении */
    int ii = 0;
    double h = 0.1;
    do
    {

        double angle = Pi / 4.0;

        double x = (0.2 + ii * h) * cos(angle);
        double y = (0.2 + ii * h) * sin(angle);

        if (ii == 0)
        {

            Profile_T << fixed << setprecision(9) << 0.2 + ii * h << "\t" << T0 << endl;

        }
        else Profile_T << fixed << setprecision(9) << 0.2 + ii * h << "\t" << Temp_prof(x, y, "Temp") << endl;

        Profile_Alfa << fixed << setprecision(9) << 0.2 + ii * h << "\t" << Temp_prof(x, y, "Alfa") << endl;

        ii++;

    } while ((0.2 + ii * h) < 1.0);

    Profile_T << fixed << setprecision(9) << 0.2 + ii * h << "\t" << T1 << endl;

    cout << "===========================================================================" << endl;

    cout << "The calculation is OVER: " << endl << File_Mesh_Name << endl;

    cout << "===========================================================================" << endl;

    Field_T.close();
    Field_Alfa.close();
    Profile_T.close();
    Profile_Alfa.close();
    Order_Acc.close();

    Write_Save();

}

void Write()
{

    if (Iter_Glob == 1) cout << fixed << setprecision(4) << "Mesh (Number of elements): " << max_el << endl;

    if (Iter_Glob == 1) cout << "The control element: El.num = "
        << num_el_1 + 1 << endl;

    if (Iter_Glob == 1) cout << "==============================================================================" << endl;

    if (Iter_Glob == 1) cout << " \t" << " \t" << "If everything is correct, then press ENTER" << endl;

    if (Iter_Glob == 1) cout << "==============================================================================" << endl;

    if (Iter_Glob == 1) cin.get();

    if (Iter_Glob == 1 || (Iter_Glob % 1000) == 0)
        cout << fixed << setprecision(4) << "Time: " << _time << " (" << dem_time << " sec) " << "\t" << setprecision(10)
        << "El=" << num_el_1 + 1 << ":(T = " << vectorElement[num_el_1].T << "; Alfa = " << vectorElement[num_el_1].Alfa
        << ") " << "\t" << "Max.Residual = " << E_T << endl;

    string _path = "Documents/Figure/El = " + to_string(max_el);

    if (Iter_Glob == 1)
    {

        CreateDirectoryA(_path.c_str(), NULL);

        ofstream file_E_T(_path + "/3. Err(time)_(El = " + to_string(max_el) + ").DAT", ios_base::trunc);
        ofstream file_Alfa_T(_path + "/4. Alfa(time)_(El = " + to_string(max_el) + ").DAT", ios_base::trunc);

        file_E_T << fixed << setprecision(6) << "Time: " << _time << "\t" << "Mesh (Number of elements): " << max_el << endl;
        file_Alfa_T << fixed << setprecision(6) << "Time: " << _time << "\t" << "Mesh (Number of elements): " << max_el << " Number of control elements: " << num_el_1 + 1 << endl;

    }

    ofstream file_E_T(_path + "/3. Err(time)_(El = " + to_string(max_el) + ").DAT", ios_base::app);
    ofstream file_Alfa_T(_path + "/4. Alfa(time)_(El = " + to_string(max_el) + ").DAT", ios_base::app);

    if (Iter_Glob == 1 || (Iter_Glob % 1000) == 0) file_E_T << fixed << setprecision(6) << _time << fixed << setprecision(10) << "\t" << E_T << endl;
    if (Iter_Glob == 1 || (Iter_Glob % 1000) == 0) file_Alfa_T << fixed << setprecision(6) << dem_time << fixed << setprecision(10) << "\t" << vectorElement[num_el_1].Alfa << endl;

    file_E_T.close();
    file_Alfa_T.close();
}

int main()
{

    Iter_Glob = 0;

    Mesh_Init();
    Arrays_Creation();
    Initial_Conditions();

    Temp_prof(xx_1, yy_1, "NULL");

    do
    {

        Iter_Glob++;

        Redistricting();
        Calculation_Temperature();
        Development();
        Write();
        Write_Save();

        _time += dt;
        dem_time = _time * C(T1) * rho * L * L / lambda;

    } while (E_T >= 0.000000001  /*dem_time <= 500*/);

    Write_End();
    Arrays_Remove();

}
