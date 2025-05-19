#include <fstream>
#include <sstream>
#include <omp.h>
#include <iostream>
#include <string>
#include <argumentum/argparse.h>
#include <vector>
#include <cmath>
#include <chrono>
#include "progressbar.hpp"
#include "misc.h"
#include "MagneticSystem.h"

using namespace std;
using namespace argumentum;


static double qmin;
static double qmax;
static double qstep;
static double angle;
static double normPosVal;
static double normMVal;

std::string filename,output,savefilename;
bool enablePG = false; // флаг для включения прогресс-бара
bool simpleFile = false; // сохранять в упрощенном формате файла

// функция возвращает вектор S перпендикулярный вектору Q
inline Vect recip(Vect & qNorm, Vect & s){
    double sc = scalar(s,qNorm);
    return {s.x - qNorm.x * sc, s.y - qNorm.y * sc};
}

int main(int argc, char* argv[])
{
    int num_threads = 0;
    int max_threads = 0;
    auto start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            max_threads = omp_get_max_threads();
        }
    }
    
    auto parser = argumentum::argument_parser{};
    auto params = parser.params();

    parser.config().program("msf")
        .description("Program to calculate the magnetic structure factor");
    params.add_parameter(filename,"-f","--filename").nargs(1).required().metavar("FILE.mfsys")
        .help("Path to text file with structure of the system. \
            Format is the mfsys file. txt or dat files are also is supported.");
    params.add_parameter(output,"-o","--output").absent("").nargs(1).metavar("FILE.txt")
        .help("New file where to save the result. By default it is the old file with '.dat' extention.");
    params.add_parameter(simpleFile,"","--simple")
        .help("Save only Real part as 2D matrix to the file. Needed to reduce the resulting file size. This file further may be directly readed by np.loadtxt() and plotted by plt.imshow()");
    params.add_parameter(savefilename,"-s","--save").absent("").nargs(1).metavar("FILE.mfsys")
        .help("If set, save scaled and rotated system to .mfsys file.");
    params.add_parameter(qmin,"","--qmin").nargs(1)
        .help("q minimal value");
    params.add_parameter(qmax,"","--qmax").nargs(1)
        .help("q maximal value");
    params.add_parameter(qstep,"","--qstep").nargs(1)
        .help("size of pixel (dQ)");
    params.add_parameter(angle,"-r","--rotate").absent(0).nargs(1)
        .help("Define the angle in degrees if needed to rotate the system");
    params.add_parameter(normPosVal,"-p","--normalisePositions").absent(1).nargs(1)
        .help("Divide (normalise) all space soordinates by this value. If set 0, use the minimal distance between spins as a value.");
    params.add_parameter(normMVal,"-m","--normaliseMoments").absent(1).nargs(1)
        .help("Divide (normalise) all M vector by this value. If set 0, use minimal length of M vector as a value.");
    params.add_parameter(enablePG,"-b","--bar")
        .help("Enables progress bar. Disabled by default.");

    auto res = parser.parse_args( argc, argv, 1 );

    if ( !res )
      return 1;

    if (output == ""){
        output = filename + ".dat";
    }

    const size_t rowsCount = (qmax-qmin)/qstep + 1;
    const size_t pixelsCount = rowsCount*rowsCount;

    MagneticSystem sys(filename);
    cout << "OpenMP threads: " << num_threads <<"; max: "<< max_threads << endl;
    cout<<"file: "<<filename<<" size: "<<sys.N()<<endl;
    cout<<"output: "<<output<<endl;
    
    if (angle!=0){
        sys.rotate(angle);
        cout<<"rotation by "<<angle<<" degrees applied"<<endl;
    }
    if (normPosVal != 1){
        normPosVal = sys.normPos(normPosVal);
        cout<<"positions normalised by "<<normPosVal<<endl;
    }
    if (normMVal != 1){
        normMVal = sys.normM(normMVal);
        cout<<"m vectors normalised by "<<normMVal<<endl;
    }
    if (savefilename != ""){
        sys.save(savefilename);
        cout<<"system saved to: "<<savefilename<<endl;
    }

    cout<<"file contains "<<pixelsCount<<" lines ("<<rowsCount<<"*"<<rowsCount<<")"<<endl;

    const size_t N = sys.N();

    vector <saveElement> saveArray(pixelsCount);

    vector <Vect> distances(N*N);
    for(std::size_t i = 0; i < N; ++i){
        for(std::size_t j = 0; j < N; ++j){
            distances[i * N + j] = {sys.parts[i].p.x - sys.parts[j].p.x, sys.parts[i].p.y - sys.parts[j].p.y};
        }
    }

    progressbar *pg;
    if (enablePG)
        pg = new progressbar(pixelsCount, enablePG);
    

    #pragma omp parallel for
    for (int i=0; i<pixelsCount; i++){
        if (enablePG){
            #pragma omp critical
                pg->update();
        }
        
        int row = i / rowsCount;
        int col = i % rowsCount;
        Vect r_temp;
        Vect q {qmin + col*qstep, qmin + row*qstep};
        if (q.x == 0. && q.y == 0.){
            saveArray[i] = {0,0,0,0};
            continue;
        } 
        Vect qNorm = normalise(q);

        vector <Vect> recipArr(N);
        for (int k=0; k<N; k++){
            recipArr[k] = recip(qNorm,sys.parts[k].m);
        }

        saveArray[i] = {q.x,q.y,0,0};
        double fourierRoot;
        for(size_t k=0; k<N; k++){
            for(size_t m=0; m<N; m++){
                fourierRoot = scalar(recipArr[k],recipArr[m]);
                saveArray[i].im += fourierRoot * sin(scalar( q, distances[k * N + m] ));
                saveArray[i].re += fourierRoot * cos(scalar( q, distances[k * N + m] ));
            }
        }
        saveArray[i].im /= N;
        saveArray[i].re /= N;
    }

    auto end = std::chrono::high_resolution_clock::now();
    // Вычисляем длительность
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << endl << "Time: " << duration.count()/1000. << " seconds" << std::endl;

    ofstream f(output.c_str());

    f<<"# ";
    for (int i=0; i<argc; i++){
        f<<argv[i]<<" ";
    }
    f<<endl;

    f<< "# Time: " << duration.count()/1000. << "s.; "<<num_threads <<" threads of "<< max_threads << std::endl;
    
    if (simpleFile){
        for (int i=0; i<pixelsCount; i++){
            if (i>0 && i%rowsCount==0)
                f<<endl;
            f<<saveArray[i].re<<"\t";
        }
    } else {
        f<<"#legend:"<<endl<<"#N\tqx\tqy\tre\tim"<<endl;
        for (int i=0; i<pixelsCount; i++){
            f<<i<<"\t"<<saveArray[i].qx<<"\t"<<saveArray[i].qy<<"\t"<<saveArray[i].re<<"\t"<<saveArray[i].im<<endl;
        }
    }

    f.close();
    if (enablePG)
        delete pg;

    return 0;
}