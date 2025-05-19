#ifndef MISC_H
#define MISC_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>

using namespace std;

//------------------структуры данных--------------------//

struct saveElement
{
    double qx;
    double qy;
    double re;
    double im;
};

struct Vect{
    double x;
    double y;
};

struct Part{
    Vect p;
    Vect m;
};

inline double distance(const Vect &a, const Vect &b)
{
    Vect d = {a.x-b.x,a.y-b.y};
    return sqrt(
        d.x * d.x +
        d.y * d.y
        );
}

inline double distance_2(const Vect &a, const Vect &b)
{
    Vect d = {a.x-b.x,a.y-b.y};
    return
        d.x * d.x +
        d.y * d.y;
}

inline double length(const Vect &d){
    return sqrt(
        d.x * d.x +
        d.y * d.y
        );
}

inline Vect normalise(const Vect &q){
    double l = length(q);
    return {q.x / l, q.y / l};
}

inline double scalar(const Vect &a, const Vect &b)
{
    return (a.x * b.x) + (a.y * b.y);
}


#endif //MISC_H