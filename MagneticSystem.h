#ifndef MAGNETICSYSTEM_H
#define MAGNETICSYSTEM_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <optional>
#include <numeric>
#include "stringfunctions.h"
#include "misc.h"

using namespace std;

class MagneticSystem
{
private:
    string _fileContent;
    int _fileVersion;
    void load_txt();
    void load_v1();
    void load_v2();
    void readFileToString(string fileName);

public:
    MagneticSystem(string fileName);
    ~MagneticSystem();

    inline size_t N() const { return parts.size(); }

    /**
     * @brief Возвращает версию файла с магнитной системой
     * 
     * @param file путь до файла
     * @return int 
     * -1 - если формат не опознан
     * 0 - если это текстовый файл с частицами
     * 1 - если это .mfsys версии 1 (https://github.com/uyras/partsEngine/wiki/Формат-файла-магнитной-системы#версия-1)
     * 1 - если это .mfsys версии 2 (https://github.com/uyras/partsEngine/wiki/Формат-файла-магнитной-системы#версия-2)
     */
    static int fileVersion(std::string file);

    void save(string filename) const;

    /**
     * @brief Rotate a point counterclockwise by a given angle around a given origin.
        The angle should be given in degrees.
    * 
    * @param angleDegree 
    * @param origin 
    */
    void rotate(double angleDegree, Vect point = {0,0});

    double normPos(double val);

    double normM(double val);

    vector< Part >   parts;
};


#endif //MAGNETICSYSTEM_H