#include "Vect.h"
#include "PartArray.h"
#include <fstream>
#include <sstream>
#include <omp.h>
#include <iostream>
#include <string>
#include <argumentum/argparse.h>
#include <vector>
#include <cmath>
#include "ProgressBar.h"

using namespace std;
using namespace argumentum;


static double qmin;
static double qmax;
static double qstep;
static double angle;
static double normPosVal;
static double normMVal;

/**
 * @brief Rotate a point counterclockwise by a given angle around a given origin.
    The angle should be given in degrees.
 * 
 * @param sys system
 * @param angleDegree 
 * @param origin 
 */
void rotate(PartArray &sys,double angleDegree, Vect point = {0,0,0}){
    double angleRadians = angleDegree * (M_PI / 180.);
    double anSin = sin(angleRadians);
    double anCos = cos(angleRadians);

    double tmpx, tmpy;
    for (auto & p : sys.parts){
        tmpx = /*p->m.x + */anCos * (point.x - p->m.x) - anSin * (point.y - p->m.y);
        tmpy = /*p->m.y + */anSin * (point.x - p->m.x) + anCos * (point.y - p->m.y);
        p->m.x = tmpx;
        p->m.y = tmpy;

        tmpx = /*p->pos.x + */anCos * (point.x - p->pos.x) - anSin * (point.y - p->pos.y);
        tmpy = /*p->pos.y + */anSin * (point.x - p->pos.x) + anCos * (point.y - p->pos.y);
        p->pos.x = tmpx;
        p->pos.y = tmpy;
    }
    return;
}

double normPos(PartArray &sys,double val){
    if (val == 0){
        bool fst = true;
        for (auto p1 : sys.parts){
            for (auto p2 : sys.parts){
                if (p1 != p2){
                    double sp = p1->pos.space(p2->pos);
                    if (fst){
                        val = sp;
                        fst = false;
                    } else {
                        if (val>sp)
                            val = sp;
                    }
                }
            }
        }
    }

    for (auto & p : sys.parts){
        p->pos.x /= val;
        p->pos.y /= val;
    }
    return val;
}

double normM(PartArray &sys,double val){
    if (val == 0){
        bool fst = true;
        for (auto p : sys.parts){
            double sp = p->m.length();
            if (fst){
                val = sp;
                fst = false;
            } else {
                if (val>sp)
                    val = sp;
            }
        }
    }

    for (auto & p : sys.parts){
        p->m.x /= val;
        p->m.y /= val;
    }
    return val;
}

// функция возвращает вектор S перпендикулярный вектору Q
inline Vect recip(Vect & qNorm, Vect & s){
    return s - (qNorm * s.scalar(qNorm));
}

inline double scalar(const Vect & a, const Vect & b) {
	return (a.x * b.x) + (a.y * b.y);
}

PartArray readTxt(string filename){
    ifstream f(filename);
    string buf;
    PartArray sys;
    while (getline(f,buf))
    {
        istringstream ss(buf);
        Part *p = new Part();
        ss>>(p->pos.x);
        ss>>(p->pos.y);
        ss>>(p->m.x);
        ss>>(p->m.y);
        sys.insert(p);
    }
    
    f.close();
    return sys;
}

int main(int argc, char* argv[])
{

    std::string filename,output,savefilename;
    
    auto parser = argumentum::argument_parser{};
    auto params = parser.params();

    parser.config().program("msf")
        .description("Program to calculate the magnetic structure factor");
    params.add_parameter(filename,"-f","--filename").nargs(1).required().metavar("FILE.mfsys")
        .help("Path to text file with structure of the system. \
            Format is the mfsys file. Or txt file is supported.");
    params.add_parameter(output,"-o","--output").absent("").nargs(1).metavar("FILE.txt")
        .help("New file where to save the result. By default it is the old file with '.dat' extention.");
    params.add_parameter(savefilename,"-s","--save").absent("").nargs(1).metavar("FILE.mfsys")
        .help("If set, save scaled and rotated system to .mfsys file.");
    params.add_parameter(qmin,"","--qmin").nargs(1).absent(-5)
        .help("q minimal value");
    params.add_parameter(qmax,"","--qmax").nargs(1).absent(5)
        .help("q maximal value");
    params.add_parameter(qstep,"","--qstep").nargs(1).absent(0.1)
        .help("size of pixel (dQ)");
    params.add_parameter(angle,"-r","--rotate").absent(0).nargs(1)
        .help("Define the angle in degrees if needed to rotate the system");
    params.add_parameter(normPosVal,"-p","--normalisePositions").absent(1).nargs(1)
        .help("Divide (normalise) all space soordinates by this value. If set 0, use the minimal distance between spins as a value.");
    params.add_parameter(normMVal,"-m","--normaliseMoments").absent(1).nargs(1)
        .help("Divide (normalise) all M vector by this value. If set 0, use minimal length of M vector as a value.");

    auto res = parser.parse_args( argc, argv, 1 );

    if ( !res )
      return 1;

    if (output == ""){
        output = filename + ".dat";
    }

    const int totsteps = (qmax-qmin)/qstep + 1;

    PartArray sys;
    if (filename.rfind(".txt") == string::npos)
        sys.load(filename);
    else
        sys = readTxt(filename);
    cout<<"file: "<<filename<<" size: "<<sys.size()<<endl;
    cout<<"output: "<<output<<endl;
    
    if (angle!=0){
        rotate(sys,angle);
        cout<<"rotation by "<<angle<<" degrees applied"<<endl;
    }
    if (normPosVal != 1){
        normPosVal = normPos(sys,normPosVal);
        cout<<"positions normalised by "<<normPosVal<<endl;
    }
    if (normMVal != 1){
        normMVal = normM(sys,normMVal);
        cout<<"m vectors normalised by "<<normMVal<<endl;
    }
    if (savefilename != ""){
        sys.save(savefilename);
        cout<<"system saved to: "<<savefilename<<endl;
    }
    ofstream f(output.c_str());
    cout<<"file contains "<<totsteps*totsteps<<" lines ("<<totsteps<<"*"<<totsteps<<")"<<endl;

    struct saveElement
    {
        double qx;
        double qy;
        double re;
        double im;
    };

    vector <saveElement> saveArray(totsteps*totsteps);

    vector <Vect> distances(sys.size()*sys.size());
    for(auto &parta : sys.parts){
        for(auto &partb : sys.parts){
            distances[parta->Id() * sys.size() + partb->Id()] = parta->pos - partb->pos;
        }
    }

    ProgressBar pg;
    pg.start(totsteps*totsteps);
    

    #pragma omp parallel for
    for (int i=0; i<totsteps*totsteps; i++){
        pg.update(i+1);
        int row = i / totsteps;
        int col = i % totsteps;
        Vect r_temp;
        Vect q( qmin + row*qstep, qmin + col*qstep, 0);
        if (q.x == 0. && q.y == 0.){
            saveArray[i] = {0,0,0,0};
            continue;
        } 
        Vect qNorm = q.normalize();

        vector <Vect> recipArr(sys.size());
        for(auto &part : sys.parts){
            recipArr[part->Id()] = recip(qNorm,part->m);
        }

        saveArray[i] = {q.x,q.y,0,0};
        double fourierRoot;
        for(auto &parta : sys.parts){
            for(auto &partb : sys.parts){
                fourierRoot = scalar(recipArr[parta->Id()],recipArr[partb->Id()]);
                saveArray[i].im += fourierRoot * sin(scalar( q, distances[parta->Id() * sys.size() + partb->Id()] ));
                saveArray[i].re += fourierRoot * cos(scalar( q, distances[parta->Id() * sys.size() + partb->Id()] ));
            }
        }
        saveArray[i].im /= sys.size();
        saveArray[i].re /= sys.size();
    }

    f<<"#legend:"<<endl<<"#N\tqx\tqy\tre\tim"<<endl;
    for (int i=0; i<totsteps*totsteps; i++){
        f<<i<<"\t"<<saveArray[i].qx<<"\t"<<saveArray[i].qy<<"\t"<<saveArray[i].re<<"\t"<<saveArray[i].im<<endl;
    }

    f.close();

    return 0;
}