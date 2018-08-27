#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
/**
 * @author ionut.brinzan
 * Run: mpirun -np NR ./filtru topologie.in imagini.in statistica.out
 */
using namespace std;
int main(int argc, char* argv[])
{

    MPI_Request request;
    MPI_Status status;
    int nProcesses;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

    //citesc vecinii mei
    ifstream input(argv[1]);
    int i;
    string line;
    for (i = 0; i < rank + 1; i++) {
        getline(input, line);
    }

    //setez vecinii mei
    vector<int> neigh;
    stringstream ss(line);
    string item;
    getline(ss, item, ' ');
    int val;
    while (getline(ss, item, ' ')) {
        sscanf(item.c_str(), "%d", &val);
        neigh.push_back(val);
    }
    input.close();

    int terminate = 0;
    int nrLines = 0;

    //vector cu nr de linii procesate
    long* stats = (long*)malloc(nProcesses * sizeof(long));
    for (i = 0; i < nProcesses; i++) {
        stats[i] = 0;
    }

    if (rank == 0) {

        ifstream input2(argv[2]);
        string line2;
        getline(input2, line2);
        int numberImages;
        int i;
        sscanf(line2.c_str(), "%d", &numberImages);

        //cat timp am imagini de prelucrat
        for (i = 0; i < numberImages; i++) {

            getline(input2, line2);
            istringstream iss(line2);
            iss >> item;
            string filter(item);
            iss >> item;
            string imageIN(item);
            iss >> item;
            string imageOUT(item);
            ifstream input3(imageIN.c_str());
            ofstream output(imageOUT.c_str());

            string line3;
            int width;
            int height;
            int t = 0;

            //citesc antetul
            while (t < 3) {
                getline(input3, line3);
                if (line3[0] != '#') {
                    if (t == 1) {
                        istringstream iss2(line3);
                        iss2 >> item;
                        sscanf(item.c_str(), "%d", &width);
                        iss2 >> item;
                        sscanf(item.c_str(), "%d", &height);
                    }
                    t++;
                }
                output << line3 << "\n";
            }

            //alocari + bordare
            int* data = (int*)malloc((height + 2) * (width + 2) * sizeof(int));
            int* data2 = (int*)malloc(height * (width + 2) * sizeof(int));
            int** pic2 = (int**)malloc(height * sizeof(int*));
            int** pic = (int**)malloc((height + 2) * sizeof(int*));

            int l;
            for (l = 0; l < height + 2; l++) {
                pic[l] = &(data[(width + 2) * l]);
            }
            for (l = 0; l < height; l++) {
                pic2[l] = &(data2[(width + 2) * l]);
            }
            for (l = 0; l < width + 2; l++) {
                pic[0][l] = 0;
                pic[height + 1][l] = 0;
            }
            for (l = 0; l < height + 2; l++) {
                pic[l][0] = 0;
                pic[l][width + 1] = 0;
            }

            long a;
            long b;
            //citesc imaginea propriu zisa
            for (a = 1; a < height + 1; a++) {
                for (b = 1; b < width + 1; b++) {
                    input3 >> pic[a][b];
                }
            }

            width += 2;
            int newPic;
            int p;

            //semnalez ca vreau sa trimit o poza ( tag = 0)
            for (p = 0; p < neigh.size(); p++) {
                MPI_Send(&newPic, 1, MPI_INT, neigh[p], 0, MPI_COMM_WORLD);
            }

            //trimit width si height a bucatii de imagine pe care un copil
            //o sa o primeasca
            //width => tag = 1
            //height => tag = 2
            int toSend = height / neigh.size();
            int last = height - (neigh.size() - 1) * toSend;
            for (p = 0; p < neigh.size(); p++) {
                if (p != neigh.size() - 1) {
                    MPI_Send(&width, 1, MPI_INT, neigh[p], 1, MPI_COMM_WORLD);
                    MPI_Send(&toSend, 1, MPI_INT, neigh[p], 2, MPI_COMM_WORLD);
                }
                else {
                    MPI_Send(&width, 1, MPI_INT, neigh[p], 1, MPI_COMM_WORLD);
                    MPI_Send(&last, 1, MPI_INT, neigh[p], 2, MPI_COMM_WORLD);
                }
            }

            //aflu ce tag pun a poza pentru a indica ce filtru trebuie aplicat
            //de nodurile frunza
            int currentFilter;
            if (filter == "sobel") {
                currentFilter = 3;
            }
            else
                currentFilter = 4;

            //trimit imaginea la copii impartita in bucati
            int line = 1;
            for (p = 0; p < neigh.size(); p++) {
                if (p != neigh.size() - 1) {
                    MPI_Send(&(pic[line - 1][0]), (toSend + 2) * width, MPI_INT, neigh[p], currentFilter, MPI_COMM_WORLD);
                    line += toSend;
                }
                else {
                    MPI_Send(&(pic[line - 1][0]), (last + 2) * width, MPI_INT, neigh[p], currentFilter, MPI_COMM_WORLD);
                }
            }

            //primesc bucatile de imagine cu filtru aplicat
            line = 0;
            for (p = 0; p < neigh.size(); p++) {
                if (p != neigh.size() - 1) {
                    MPI_Recv(&(pic2[line][0]), toSend * width, MPI_INT, neigh[p], 5, MPI_COMM_WORLD, &status);
                    line += toSend;
                }
                else {
                    MPI_Recv(&(pic2[line][0]), last * width, MPI_INT, neigh[p], 5, MPI_COMM_WORLD, &status);
                }
            }

            //scriu in imaginea de output
            for (a = 0; a < height; a++) {
                for (b = 1; b < width - 1; b++) {
                    output << pic2[a][b] << "\n";
                }
            }
            delete[] data;
            delete[] data2;
            delete[] pic;
            delete[] pic2;
            input3.close();
            output.close();
        }
        //am terminat de prelucrat poze => cer statistica
        int p;
        int dummy;
        //semnalizez ca vreau sa primesc statistica(tag == 9)
        for (p = 0; p < neigh.size(); p++) {
            MPI_Send(&dummy, 1, MPI_INT, neigh[p], 9, MPI_COMM_WORLD);
        }
        //primesc statistica
        for (p = 0; p < neigh.size(); p++) {
            long* dummyStat = (long*)calloc(nProcesses, sizeof(long));
            MPI_Recv(dummyStat, nProcesses, MPI_LONG, neigh[p], 8, MPI_COMM_WORLD, &status);
            int i;
            for (int i = 0; i < nProcesses; i++) {
                if (dummyStat[i] != 0) {
                    stats[i] = dummyStat[i];
                }
            }
            delete[] dummyStat;
        }
        //scriu statistica in fisier
        ofstream out(argv[3]);
        for (p = 0; p < nProcesses; p++) {
            out << p << ": " << stats[p] << "\n";
        }
        neigh.clear();
        delete[] stats;
        out.close();
        input2.close();
    }
    else {
        int parent;
        while (terminate == 0) {
            int newPic;
            //astept sa primesc ceva
            MPI_Recv(&newPic, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int tag = status.MPI_TAG;
            //daca am imagine de procesat
            if (tag == 0) {
            	//retin parintele
                parent = status.MPI_SOURCE;
                int p;
                //le semnalez copiilor ca o sa vina o imagine
                for (p = 0; p < neigh.size(); p++) {
                    if (neigh[p] != parent) {
                        MPI_Send(&newPic, 1, MPI_INT, neigh[p], 0, MPI_COMM_WORLD);
                    }
                }
                int width;
                int height;
                //primesc width si height cu tagurile specifice
                MPI_Recv(&width, 1, MPI_INT, parent, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&height, 1, MPI_INT, parent, 2, MPI_COMM_WORLD, &status);
                //daca nu sunt nod frunza
                if (neigh.size() > 1) {
                	//trimit la copii cat o sa aiba bucata de imagine pe care ei o vor procesa
                    int toSend = height / (neigh.size() - 1);
                    int last = height - (neigh.size() - 2) * toSend;
                    for (p = 0; p < neigh.size(); p++) {
                        if (neigh[p] != parent) {
                            if (p == neigh.size() - 1) {
                                MPI_Send(&width, 1, MPI_INT, neigh[p], 1, MPI_COMM_WORLD);
                                MPI_Send(&last, 1, MPI_INT, neigh[p], 2, MPI_COMM_WORLD);
                            }
                            else {
                                MPI_Send(&width, 1, MPI_INT, neigh[p], 1, MPI_COMM_WORLD);
                                MPI_Send(&toSend, 1, MPI_INT, neigh[p], 2, MPI_COMM_WORLD);
                            }
                        }
                    }
                    //aloc memorie pentru ce o sa primesc eu
                    int* data = (int*)malloc((height + 2) * (width) * sizeof(int));
                    int** mypic = (int**)malloc((height + 2) * sizeof(int*));
                    int* data2 = (int*)malloc(height * (width) * sizeof(int));
                    int** mypic2 = (int**)malloc(height * sizeof(int*));
                    int l;
                    for (l = 0; l < height + 2; l++) {
                        mypic[l] = &(data[width * l]);
                    }
                    for (l = 0; l < height; l++) {
                        mypic2[l] = &(data2[width * l]);
                    }
                    //astept sa primesc poza
                    MPI_Recv(&(mypic[0][0]), width * (height + 2), MPI_INT, parent, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    int filter = status.MPI_TAG;
                    int line = 1;
                    //trimit la copii
                    for (p = 0; p < neigh.size(); p++) {
                        if (neigh[p] != parent) {
                            if (p != neigh.size() - 1) {
                                MPI_Send(&(mypic[line - 1][0]), (toSend + 2) * width, MPI_INT, neigh[p], filter, MPI_COMM_WORLD);
                                line += toSend;
                            }
                            else {
                                MPI_Send(&(mypic[line - 1][0]), (last + 2) * width, MPI_INT, neigh[p], filter, MPI_COMM_WORLD);
                            }
                        }
                    }
                    line = 0;
                    //primesc de la copii imaginea procesata
                    for (p = 0; p < neigh.size(); p++) {
                        if (neigh[p] != parent) {
                            if (p != neigh.size() - 1) {
                                MPI_Recv(&(mypic2[line][0]), toSend * width, MPI_INT, neigh[p], 5, MPI_COMM_WORLD, &status);
                                line += toSend;
                            }
                            else {
                                MPI_Recv(&(mypic2[line][0]), last * width, MPI_INT, neigh[p], 5, MPI_COMM_WORLD, &status);
                            }
                        }
                    }
                    //trimit la parinte ce am primit de la toti copiii
                    MPI_Send(&(mypic2[0][0]), width * height, MPI_INT, parent, 5, MPI_COMM_WORLD);
                    delete[] data;
                    delete[] data2;
                    delete[] mypic;
                    delete[] mypic2;
                }
                else {
                	//daca sunt frunza
                	//aloc memorie pentru ce urmeaza sa primesc
                    int* data = (int*)malloc((height + 2) * (width) * sizeof(int));
                    int** mypic = (int**)malloc((height + 2) * sizeof(int*));
                    int* data2 = (int*)malloc((height + 2) * (width) * sizeof(int));
                    int** mypic2 = (int**)malloc((height + 2) * sizeof(int*));
                    int l;
                    for (l = 0; l < height + 2; l++) {
                        mypic[l] = &(data[width * l]);
                    }
                    for (l = 0; l < height + 2; l++) {
                        mypic2[l] = &(data2[width * l]);
                    }
                    //astept sa primesc poza
                    MPI_Recv(&(mypic[0][0]), width * (height + 2), MPI_INT, parent, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    int filter = status.MPI_TAG;
                    int a;
                    int b;
                    int h = height + 1;
                    int w = width - 1;
                    //aplic filtru in functie de tag
                    if (filter == 3) {
                        for (a = 1; a < h; a++) {
                            for (b = 1; b < w; b++) {
                                mypic2[a][b] = 1 * mypic[a - 1][b - 1] - 1 * mypic[a - 1][b + 1];
                                mypic2[a][b] += 2 * mypic[a][b - 1] - 2 * mypic[a][b + 1];
                                mypic2[a][b] += 1 * mypic[a + 1][b - 1] - 1 * mypic[a + 1][b + 1];
                                mypic2[a][b] += 127;
                                if (mypic2[a][b] > 255) {
                                    mypic2[a][b] = 255;
                                }
                                else if (mypic2[a][b] < 0) {
                                    mypic2[a][b] = 0;
                                }
                            }
                        }
                    }
                    else if (filter == 4) {
                        for (a = 1; a < h; a++) {
                            for (b = 1; b < w; b++) {
                                mypic2[a][b] = -(mypic[a - 1][b - 1] + mypic[a - 1][b] + mypic[a - 1][b + 1]);
                                mypic2[a][b] += -1 * mypic[a][b - 1] + 9 * mypic[a][b] - 1 * mypic[a][b + 1];
                                mypic2[a][b] += -(mypic[a + 1][b - 1] + mypic[a + 1][b] + mypic[a + 1][b + 1]);
                                if (mypic2[a][b] > 255) {
                                    mypic2[a][b] = 255;
                                }
                                else if (mypic2[a][b] < 0) {
                                    mypic2[a][b] = 0;
                                }
                            }
                        }
                    }
                    nrLines += height;
                    //trimit la parinte
                    MPI_Send(&(mypic2[1][0]), width * height, MPI_INT, parent, 5, MPI_COMM_WORLD);
                    delete[] data;
                    delete[] data2;
                    delete[] mypic;
                    delete[] mypic2;
                }
            }
            else if (tag == 9) {
            	//daca e tag de terminare
            	//daca nu sunt frunza
                if (neigh.size() > 1) {
                    int dummy;
                    int p;
                    //trimit la copii ca vreau statistica lor
                    for (p = 0; p < neigh.size(); p++) {
                        if (neigh[p] != parent) {
                            MPI_Send(&dummy, 1, MPI_INT, neigh[p], 9, MPI_COMM_WORLD);
                        }
                    }
                    //primesc statistica lor si updatez statistica mea
                    for (p = 0; p < neigh.size(); p++) {
                        long* dummyStat = (long*)calloc(nProcesses, sizeof(long));
                        if (neigh[p] != parent) {
                            MPI_Recv(dummyStat, nProcesses, MPI_LONG, neigh[p], 8, MPI_COMM_WORLD, &status);
                        }
                        int j;
                        for (j = 0; j < nProcesses; j++) {
                            if (dummyStat[j] != 0) {
                                stats[j] = dummyStat[j];
                            }
                        }
                        delete[] dummyStat;
                    }
                    //trimit noua statistica
                    MPI_Send(stats, nProcesses, MPI_LONG, parent, 8, MPI_COMM_WORLD);
                    delete[] stats;
                    terminate = 1;
                }
                else {
                	//daca sunt frunza
                    stats[rank] = nrLines;
                    //trimit statistica mea la parinte
                    MPI_Send(stats, nProcesses, MPI_LONG, parent, 8, MPI_COMM_WORLD);
                    terminate = 1;
                    delete[] stats;
                }
            }
        }
    }

    MPI_Finalize();
    return 0;
}