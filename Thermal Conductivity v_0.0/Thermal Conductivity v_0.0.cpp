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
        cout << "ERROR opening the file: "  << File_Mesh_Name << endl << endl << endl;

        exit(1);
    }
    else
    {
        if (max_str != 0) cout << "File opened SUCCESSFULLY: " << File_Mesh_Name << endl;
        if (max_str != 0) cout << fixed << setprecision(4) << "Time: " << _time << "\t" << "Mesh (Number of nodes): " << max_str << endl;
        //return 1;
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
        el.Area_el = 0.25 * sqrt((a + b + c) + (b + c - a) + (a + c - b) + (a + b - c));

        double xc_1, ус_1, xc_2, ус_2, xc_3, ус_3;

        /* Блок нахождения центра элемента */
        {

            xc_1 = 0.5 * (el.Coord_vert_2.x + el.Coord_vert_1.x);
            ус_1 = 0.5 * (el.Coord_vert_2.y + el.Coord_vert_1.y);

            xc_2 = 0.5 * (el.Coord_vert_3.x + el.Coord_vert_2.x);
            ус_2 = 0.5 * (el.Coord_vert_3.y + el.Coord_vert_2.y);

            xc_3 = 0.5 * (el.Coord_vert_1.x + el.Coord_vert_3.x);
            ус_3 = 0.5 * (el.Coord_vert_1.y + el.Coord_vert_3.y);

            /* Уравнение серединного перпендикуляра 1 */
            double a1 = el.Coord_vert_1.x - el.Coord_vert_2.x;
            double b1 = -(el.Coord_vert_2.y - el.Coord_vert_1.y);
            double c1 = -((el.Coord_vert_1.x - el.Coord_vert_2.x) * xc_1 - (el.Coord_vert_2.y - el.Coord_vert_1.y) * ус_1);

            double a2 = el.Coord_vert_2.x - el.Coord_vert_3.x;
            double b2 = -(el.Coord_vert_3.y - el.Coord_vert_2.y);
            double c2 = -((el.Coord_vert_2.x - el.Coord_vert_3.x) * xc_2 - (el.Coord_vert_3.y - el.Coord_vert_2.y) * ус_2);

            el.Coord_center_el.x = (b1 * c2 - b2 * c1) / (b2 * a1 - b1 * a2);
            el.Coord_center_el.y = (a2 * c1 - a1 * c2) / (a1 * b2 - a2 * b1);

        }

        el.h_1 = sqrt(pow((el.Coord_center_el.x - xc_1), 2) + pow((el.Coord_center_el.y - ус_1), 2));
        el.h_2 = sqrt(pow((el.Coord_center_el.x - xc_2), 2) + pow((el.Coord_center_el.y - ус_2), 2));
        el.h_3 = sqrt(pow((el.Coord_center_el.x - xc_3), 2) + pow((el.Coord_center_el.y - ус_3), 2));

        vectorElement.push_back(el);

        i = i;

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

            //if(vectorElement[i].Geom_el == 2 && vectorElement[i].Num_bound == 2) cout << vectorElement[i].Num_el << " \t " << vectorElement[i].N1_e << " \t " << vectorElement[i].N2_e << " \t " << vectorElement[i].N3_e << " \t " << vectorElement[i].Num_bound << endl;

            i = i;

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

        if(vectorElement[i].Geom_el == 2) 
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

        if (vectorElement[i].Num_bound == 1 /*|| vectorElement[i].Num_bound == 2*/)
        {

            vectorElement[i].T = 0.0;
        }

        if (vectorElement[i].Num_bound == 2 /*|| vectorElement[i].Num_bound == 4*/)
        {

            vectorElement[i].T = 1.0;
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

            i = i;

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

void Write()
{

    if (Iter_Glob == 1) cout << fixed << setprecision(4) << "Time: " << _time << "\t" << "Mesh (Number of elements): " << max_str << endl;

    if(Iter_Glob == 1) cin.get();

    cout << fixed << setprecision(4) << "Time: " << _time << "\t" << setprecision(10)
    << "T (El=" << num_el << ") = " << vectorElement[num_el - 1].T << "\t" << "Max.Residual = " << E_T << endl;

    if (Iter_Glob == 1)
    {

        ofstream file_E_T("Documents/Figure/E_T.DAT", ios_base::trunc);
        file_E_T << fixed << setprecision(6) << "Time: " << _time << "\t" << "Mesh (Number of elements): " << max_str << endl;

    }

    ofstream file_E_T("Documents/Figure/E_T.DAT", ios_base::app);
    file_E_T << fixed << setprecision(6) << _time << "\t" << E_T << endl;

}

void Write_End()
{

    ofstream Field_T("Documents/Figure/T_Field.DAT");
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
    for (int i = 0; i < max_str; i++)
    {

        if (vectorElement[i].Geom_el == 2 && 
           (vectorElement[i].Num_el == 1420 || vectorElement[i].Num_el == 181 || 
            vectorElement[i].Num_el == 219 || vectorElement[i].Num_el == 2536 || 
            vectorElement[i].Num_el == 1570 || vectorElement[i].Num_el == 1997 || 
            vectorElement[i].Num_el == 2160 || vectorElement[i].Num_el == 800 ||
            vectorElement[i].Num_el == 1331 || vectorElement[i].Num_el == 1770 ||
            vectorElement[i].Num_el == 884 || vectorElement[i].Num_el == 1501 || 
            vectorElement[i].Num_el == 1465 || vectorElement[i].Num_el == 1914 ||
            vectorElement[i].Num_el == 2707 || vectorElement[i].Num_el == 2385 ||
            vectorElement[i].Num_el == 1023))
        {

            vectorElement[i].l = sqrt(pow((vectorElement[i].Coord_center_el.x), 2) + pow((vectorElement[i].Coord_center_el.y), 2));

            Profile_T << fixed << setprecision(10) << vectorElement[i].l << " \t " << vectorElement[i].T << " \t " << vectorElement[i].Num_el << endl;

        }            

    }

}

void Find_Num_el(double xx, double yy)
{

    double dl;
    double temp = 100;

    for (int i = 0; i < vectorElement.size(); i++)
    {

        if (vectorElement[i].Geom_el == 2)
        {

            dl = sqrt(pow((vectorElement[i].Coord_center_el.x - xx), 2) + pow((vectorElement[i].Coord_center_el.y - yy), 2));

            if (dl < temp)
            {

                temp = dl;
                num_el = vectorElement[i].Num_el;

                i = i;

            }

        }

    }

}

int main()
{

    Mesh_Init();
    Arrays_Creation();
    Initial_Conditions();

    Find_Num_el(0.47919608314923778, 0.14015797430309276);

    Iter_Glob = 0;

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
