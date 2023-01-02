#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

#include "Variables.h"

void Arrays_Creation()
{

    T = new double[max_str + 1]; // удалить

}

void Arrays_Remove()
{

    delete[] T; //удалить

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

            cout << fixed << setprecision(4) << "Time: " << _time << endl << "Mesh (Number of nodes): " << max_str << endl;

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

    Point node;

    /* Заполнение структуры точек */
    for (int i = 1; i <= max_str; i++)
    {

        double temp;

        File_Mesh >> node.Num_node >> node.x >> node.y >> temp;

        node.Boundary = false;

        vectorPoint.push_back(node);

        i = i;

    }

    Find_String("$Elements");

    Element el;

    /* Заполнение структуры элементов */
    for (int i = 1; i <= max_str + 1; i++)
    {

        double temp;

        File_Mesh >> el.Num_el >> el.Geom_el;

        if (el.Geom_el == 15)
        {

            File_Mesh >> temp >> temp >> temp >> el.Num_vert_1;
            el.Coord_vert_1 = vectorPoint[(el.Num_vert_1) - 1];
            vectorPoint[(el.Num_vert_1) - 1].Boundary = true;

        }
        if (el.Geom_el == 1)
        {

            File_Mesh >> temp >> temp >> el.Num_bound >> el.Num_vert_1 >> el.Num_vert_2;

            el.Coord_vert_1 = vectorPoint[(el.Num_vert_1) - 1];
            el.Coord_vert_2 = vectorPoint[(el.Num_vert_2) - 1];
            vectorPoint[(el.Num_vert_1) - 1].Boundary = true;
            vectorPoint[(el.Num_vert_2) - 1].Boundary = true;

            el.Num_bound_vert_1 = el.Num_vert_1;
            el.Num_bound_vert_2 = el.Num_vert_2;

            i = i;

        }
        if (el.Geom_el == 2)
        {

            File_Mesh >> temp >> temp >> temp >> el.Num_vert_1 >> el.Num_vert_2 >> el.Num_vert_3;

            el.Coord_vert_1 = vectorPoint[(el.Num_vert_1) - 1];
            el.Coord_vert_2 = vectorPoint[(el.Num_vert_2) - 1];
            el.Coord_vert_3 = vectorPoint[(el.Num_vert_3) - 1];

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

        if (i == 51)
        {

            i = i;

        }

    }

    int p1, p2;

    /* Блок нахождения соседних и граничных элементов */
    {

        for (int i = 0; i < vectorElement.size() - 1; i++)
        {

            if (vectorElement[i].Geom_el == 2)
            {

                /* Cосед между точками 1 и 2 */
                p1 = vectorElement[i].Num_vert_1;
                p2 = vectorElement[i].Num_vert_2;
                if (vectorPoint[p1 - 1].Boundary && vectorPoint[p2 - 1].Boundary)
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
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_2 || p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_1) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N1_e = j + 1; break; }
                        if ((p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_2) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N1_e = j + 1; break; }
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_1) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N1_e = j + 1; break; }
                    }

                }

                /* Cосед между точками 2 и 3 */
                p1 = vectorElement[i].Num_vert_2;
                p2 = vectorElement[i].Num_vert_3;
                if (vectorPoint[p1 - 1].Boundary && vectorPoint[p2 - 1].Boundary)
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
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_2 || p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_1) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N2_e = j + 1; break; }
                        if ((p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_2) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N2_e = j + 1; break; }
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_1) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N2_e = j + 1; break; }
                    }

                }

                /* Cосед между точками 1 и 3 */
                p1 = vectorElement[i].Num_vert_1;
                p2 = vectorElement[i].Num_vert_3;
                if (vectorPoint[p1 - 1].Boundary && vectorPoint[p2 - 1].Boundary)
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
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_2 || p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_1) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N3_e = j + 1; break; }
                        if ((p1 == vectorElement[j].Num_vert_2 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_2) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N3_e = j + 1; break; }
                        if ((p1 == vectorElement[j].Num_vert_1 && p2 == vectorElement[j].Num_vert_3 || p1 == vectorElement[j].Num_vert_3 && p2 == vectorElement[j].Num_vert_1) && (j + 1) != vectorElement[i].Num_el) { vectorElement[i].N3_e = j + 1; break; }

                        j = j;

                    }

                }

            }

            //cout << vectorElement[i].Num_el << "; " << vectorElement[i].N1_e << "; " << vectorElement[i].N2_e << "; " << vectorElement[i].N3_e << endl;

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
    for (int i = 0; i <= max_str; i++)
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
    for (int i = 0; i <= max_str; i++)
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

    for (int i = 0; i <= max_str; i++)
    {

        vectorElement[i].t = vectorElement[i].T;

        if (vectorElement[i].Geom_el == 2)
            i = i;

    }

}

void Initial_Conditions()
{

    Blank();

    for (int i = 0; i <= max_str; i++)
    {

        vectorElement[i].t = 0.0;
        vectorElement[i].T = 0.0;

    }

}

void Boundary_Conditions()
{

    for (int i = 0; i <= max_str; i++)
    {

        /* Первое граничное условие */
        if (vectorElement[i].Num_bound == 1)
        {

            if (vectorElement[i].N1_e == -1)
            {

                double dT1 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / vectorElement[i].h_1;
                double dT2 = (vectorElement[(vectorElement[i].N2_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[(vectorElement[i].N2_e) - 1].h_2);
                double dT3 = (vectorElement[(vectorElement[i].N3_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[(vectorElement[i].N3_e) - 1].h_3);

                vectorElement[i].T = dt / vectorElement[i].Area_el * (dT1 + dT2 + dT3) + vectorElement[i].t;

                i = i;

            }

            if (vectorElement[i].N2_e == -1)
            {

                double dT1 = (vectorElement[(vectorElement[i].N1_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[(vectorElement[i].N1_e) - 1].h_1);
                double dT2 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / vectorElement[i].h_2;
                double dT3 = (vectorElement[(vectorElement[i].N3_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[(vectorElement[i].N3_e) - 1].h_3);

                vectorElement[i].T = dt / vectorElement[i].Area_el * (dT1 + dT2 + dT3) + vectorElement[i].t;

                i = i;

            }

            if (vectorElement[i].N3_e == -1)
            {

                double dT1 = (vectorElement[(vectorElement[i].N1_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[(vectorElement[i].N1_e) - 1].h_1);
                double dT2 = (vectorElement[(vectorElement[i].N2_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[(vectorElement[i].N2_e) - 1].h_2);
                double dT3 = (T1 - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / vectorElement[i].h_3;

                vectorElement[i].T = dt / vectorElement[i].Area_el * (dT1 + dT2 + dT3) + vectorElement[i].t;

                if (vectorElement[i].T < 0)
                {
                    i = i;
                }

            }

        }

        /* Второе граничное условие */
        if (vectorElement[i].Num_bound == 2)
        {

            if (vectorElement[i].N1_e == -1)
            {

                double dT1 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / vectorElement[i].h_1;
                double dT2 = (vectorElement[(vectorElement[i].N2_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[(vectorElement[i].N2_e) - 1].h_2);
                double dT3 = (vectorElement[(vectorElement[i].N3_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[(vectorElement[i].N3_e) - 1].h_3);

                vectorElement[i].T = dt / vectorElement[i].Area_el * (dT1 + dT2 + dT3) + vectorElement[i].t;

                i = i;

            }

            if (vectorElement[i].N2_e == -1)
            {

                double dT1 = (vectorElement[(vectorElement[i].N1_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[(vectorElement[i].N1_e) - 1].h_1);
                double dT2 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / vectorElement[i].h_2;
                double dT3 = (vectorElement[(vectorElement[i].N3_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[(vectorElement[i].N3_e) - 1].h_3);

                vectorElement[i].T = dt / vectorElement[i].Area_el * (dT1 + dT2 + dT3) + vectorElement[i].t;

                i = i;

            }

            if (vectorElement[i].N3_e == -1)
            {

                double dT1 = (vectorElement[(vectorElement[i].N1_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[(vectorElement[i].N1_e) - 1].h_1);
                double dT2 = (vectorElement[(vectorElement[i].N2_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[(vectorElement[i].N2_e) - 1].h_2);
                double dT3 = (T0 - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / vectorElement[i].h_3;

                vectorElement[i].T = dt / vectorElement[i].Area_el * (dT1 + dT2 + dT3) + vectorElement[i].t;

                if (i == 65)
                    i = i;

            }

        }

    }

}

void Calculation_Temperature()
{

    for (int i = 0; i < max_str; i++)
    {

        if (vectorElement[i].N1_e != -1 && vectorElement[i].N2_e != -1 && vectorElement[i].N3_e != -1 && vectorElement[i].Geom_el == 2)
        {

            double dT1 = (vectorElement[(vectorElement[i].N1_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_1 / (vectorElement[i].h_1 + vectorElement[(vectorElement[i].N1_e) - 1].h_1);
            double dT2 = (vectorElement[(vectorElement[i].N2_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_2 / (vectorElement[i].h_2 + vectorElement[(vectorElement[i].N2_e) - 1].h_2);
            double dT3 = (vectorElement[(vectorElement[i].N3_e) - 1].t - vectorElement[i].t) * vectorElement[i].Length_face_el_3 / (vectorElement[i].h_3 + vectorElement[(vectorElement[i].N3_e) - 1].h_3);

            vectorElement[i].T = dt / vectorElement[i].Area_el * (dT1 + dT2 + dT3) + vectorElement[i].t;

            if (vectorElement[i].T < 0)
            {
                i = i;
            }

        }

        if (vectorElement[i].Geom_el == 2)
        {

            //cout << "T=" << vectorElement[i].T << " (EL=" << vectorElement[i].Num_el << "); " << endl;

        }

    }

    double q = 0.0;

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

void Write()
{

    if (Iter_Glob == 1) cout << fixed << setprecision(4) << "Mesh (Number of elements): " << max_str << endl;

    if (Iter_Glob == 1) cout << "The control element and its neighbors: El.num = "
        << num_el_1 << "; El.num = " << num_el_2 << "; El.num = " << num_el_3 << endl;

    if (Iter_Glob == 1) cout << "==============================================================================" << endl;

    if (Iter_Glob == 1) cout << " \t" << " \t" << "If everything is correct, then press ENTER" << endl;

    if (Iter_Glob == 1) cout << "==============================================================================" << endl;

    if (Iter_Glob == 1) cin.get();

    cout << fixed << setprecision(4) << "Time: " << _time << "\t" << setprecision(10)
        << "T (El=" << num_el_1 << ") = " << vectorElement[num_el_1 - 1].T << "\t" << "Max.Residual = " << E_T << endl;

    if (Iter_Glob == 1)
    {

        ofstream file_E_T("Documents/Figure/E_T.DAT", ios_base::trunc);
        file_E_T << fixed << setprecision(6) << "Time: " << _time << "\t" << "Mesh (Number of elements): " << max_str << endl;

    }

    ofstream file_E_T("Documents/Figure/E_T.DAT", ios_base::app);
    file_E_T << fixed << setprecision(6) << _time << fixed << setprecision(10) << "\t" << E_T << endl;

}

double Temp_prof(double xx, double yy)
{

    double dl;
    double temp_1 = 99, temp_2 = 100, temp_3 = 101;

    if (Iter_Glob == 0)
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
        double Ax = vectorElement[num_el_1 - 1].Coord_center_el.x;
        double Ay = vectorElement[num_el_1 - 1].Coord_center_el.y;
        double Ex = xx;
        double Ey = yy;
        double Bx = vectorElement[num_el_2 - 1].Coord_center_el.x;
        double By = vectorElement[num_el_2 - 1].Coord_center_el.y;
        double Cx = vectorElement[num_el_3 - 1].Coord_center_el.x;
        double Cy = vectorElement[num_el_3 - 1].Coord_center_el.y;

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

        double T_D = vectorElement[num_el_2 - 1].T + BD / BC * (vectorElement[num_el_3 - 1].T - vectorElement[num_el_2 - 1].T);

        double AE = sqrt(pow((Ax - xx), 2) + pow((Ay - yy), 2));
        double AD = sqrt(pow((Ax - Dx), 2) + pow((Ay - Dy), 2));
        double T_E;

        return T_E = vectorElement[num_el_1 - 1].T + AE / AD * (T_D - vectorElement[num_el_1 - 1].T);

        double pp = 0;

    }

}

void Write_End()
{

    ofstream Field_T("Documents/Figure/T_Field_(El = " + to_string(max_str) + ").DAT");
    ofstream Profile_T("Documents/Figure/T_Profile_(El = " + to_string(max_str) + ").DAT");

    Field_T << fixed << setprecision(4) << "Time: " << _time << "\t" << "Mesh (Number of cells): " << max_str << endl;

    /* Запись распределния поля температуры */
    for (int i = 0; i < max_str; i++)
    {

        if (vectorElement[i].Geom_el == 2)
            Field_T << fixed << setprecision(10) << vectorElement[i].Coord_center_el.x << " \t " << vectorElement[i].Coord_center_el.y << " \t "
            << vectorElement[i].T << " \t " << vectorElement[i].t << " \t " << vectorElement[i].Num_el << endl;

    }

    Profile_T << fixed << setprecision(4) << "Time: " << _time << "\t" << "Mesh (Number of cells): " << max_str << endl;

    /* Запись значения температуры в сечении */
    int ii = 0;
    double h = 0.1;
    do
    {

        double angle = Pi / 4.0;

        double x = (0.2 + ii * h) * cos(angle);
        double y = (0.2 + ii * h) * sin(angle);

        if (ii == 0)
        {

            Profile_T << fixed << setprecision(9) << 0.2 + ii * h << "\t" << "1.000000000" << endl;

        }
        else Profile_T << fixed << setprecision(9) << 0.2 + ii * h << "\t" << Temp_prof(x, y) << endl;

        ii++;

    } while ((0.2 + ii * h) < 1.0);

    Profile_T << fixed << setprecision(9) << 0.2 + ii * h << "\t" << "0.000000000" << endl;

    cout << "===========================================================================" << endl;

    cout << "The calculation is OVER: " << endl << File_Mesh_Name << endl;

    cout << "===========================================================================" << endl;

}

int main()
{

    Iter_Glob = 0;

    Mesh_Init();
    Arrays_Creation();
    Initial_Conditions();

    Temp_prof(xx_1, yy_1);

    do
    {

        Iter_Glob++;

        Redistricting();
        Boundary_Conditions();
        Calculation_Temperature();
        Development();
        Write();

        _time += dt;

    } while (E_T > 0.000001);

    Write_End();
    Arrays_Remove();

}
