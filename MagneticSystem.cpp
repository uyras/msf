#include "MagneticSystem.h"

MagneticSystem::MagneticSystem(string fileName)
{
    _fileVersion = fileVersion(fileName); 
    readFileToString(fileName); //кэшируем содержимое файла во внутреннюю строку

    if (this->_fileVersion == 0) this->load_txt();
    if (this->_fileVersion == 1) this->load_v1();
    if (this->_fileVersion == 2) this->load_v2();

}

MagneticSystem::~MagneticSystem()
{
}

void MagneticSystem::load_txt(){
    stringstream f(_fileContent);
    string buf;
    this->parts.clear(); //удаляем все частицы
    while (getline(f,buf))
    {
        istringstream ss(buf);
        Part part;
        ss>>(part.p.x);
        ss>>(part.p.y);
        ss>>(part.m.x);
        ss>>(part.m.y);
        parts.push_back(part);
    }
    return;
}

void MagneticSystem::load_v1()
{
    stringstream f(_fileContent);

    this->parts.clear(); //удаляем все частицы

    //сначала сохраняем xyz
    double dummy;
    f >> dummy;
    f >> dummy;
    f >> dummy;

    int i=0;

    //пропускаем строку с заголовками
    char c[256];
    f.getline(c,256,'\n');
    f.getline(c,256,'\n');

    //затем читаем все магнитные моменты системы и положения точек
    double radius = 0;
    string shape;
    while (!f.eof()) {
        Part temp;
        if (!(f >> temp.p.x).good()) break; //если не получилось считать - значит конец файла
        f >> temp.p.y;
        f >> dummy;
        f >> temp.m.x;
        f >> temp.m.y;
        f >> dummy;
        f >> dummy; //w
        f >> dummy; //h
        //f >> temp.sector; для MPI реализации, @todo потом перегрузить
        f >> dummy; //r

        f >> shape;

        parts.push_back(temp);
        i++;
    }
}

void MagneticSystem::load_v2()
{
    this->parts.clear(); //удаляем все частицы

    double dummy;

    stringstream f(_fileContent);
  
    f.seekg(0);
    string section = "[parts]";
    std::string str;
    while (!f.eof() && str != section){
        std::getline(f,str);
        rtrim(str);
    }
    if (str!=section)
        throw(string("section [parts] not found in file")); //todo сделать нормальные классы для исключений


    while (!f.eof()){
        getline(f,str);
        trim(str);
        if (str.empty()) continue;
        if (str[0]=='[' && str[str.length()-1]==']') break;
        stringstream helper(str);

        Part temp;
        helper >> dummy; //id

        helper >> temp.p.x;
        helper >> temp.p.y;
        helper >> dummy;
        helper >> temp.m.x;
        helper >> temp.m.y;
        helper >> dummy;
        helper >> dummy; //state

        parts.push_back(temp);
    }
}

void MagneticSystem::readFileToString(string fileName)
{
    ifstream inFile;
    inFile.open(fileName);
    if (!inFile.good()) throw(string("Error reading file ") + fileName);
    std::stringstream strStream;
    strStream << inFile.rdbuf(); //read the file
    _fileContent.clear();
    _fileContent = strStream.str(); //str holds the content of the file
    inFile.close();
}

void MagneticSystem::save(string filename) const
{
    ofstream f;
    f.open(filename, ios_base::out|ios_base::trunc);
    if (f.fail())
        throw(string("saveHelper: file "+filename+" is unwritable or not found"));
    
    f<<"[header]"<<endl;
    f<<"version=2"<<endl;
    f<<"dimensions=3"<<endl;
    f<<"type=standart"<<endl;
    f<<"size="+std::to_string(this->N())<<endl;
    f<<"state="+string(this->N(),'0')<<endl;
    f<<"interactionrange=1"<<endl;
    f<<"sizescale=1"<<endl;
    f<<"magnetizationscale=1"<<endl;
    f<<"[parts]"<<endl;

    for (size_t i=0; i<this->N(); i++) {
        f << i << "\t";
        f << parts[i].p.x << "\t";
        f << parts[i].p.y << "\t";
        f << 0 << "\t";
        f << parts[i].m.x << "\t";
        f << parts[i].m.y << "\t";
        f << 0 << "\t";
        f << "0";
        f<<endl;
    }

    f.close();
}

void MagneticSystem::rotate(double angleDegree, Vect point)
{
    double angleRadians = angleDegree * (M_PI / 180.);
    double anSin = sin(angleRadians);
    double anCos = cos(angleRadians);

    double tmpx, tmpy;
    for (auto & p : parts){
        tmpx = /*p->m.x + */anCos * (point.x - p.m.x) - anSin * (point.y - p.m.y);
        tmpy = /*p->m.y + */anSin * (point.x - p.m.x) + anCos * (point.y - p.m.y);
        p.m.x = tmpx;
        p.m.y = tmpy;

        tmpx = /*p->pos.x + */anCos * (point.x - p.p.x) - anSin * (point.y - p.p.y);
        tmpy = /*p->pos.y + */anSin * (point.x - p.p.x) + anCos * (point.y - p.p.y);
        p.p.x = tmpx;
        p.p.y = tmpy;
    }
    return;
}

double MagneticSystem::normPos(double val)
{
    if (val == 0){
        bool fst = true;
        for (int i=0; i<this->N(); i++){
            for (int j=0; j<this->N(); j++){
                if (i != j){
                    double sp = distance(this->parts[i].p, this->parts[j].p);
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

    for (auto & part : this->parts){
        part.p.x /= val;
        part.p.y /= val;
    }
    return val;
}

double MagneticSystem::normM(double val)
{
    if (val == 0){
        bool fst = true;
        for (auto p : this->parts){
            double sp = length(p.m);
            if (fst){
                val = sp;
                fst = false;
            } else {
                if (val>sp)
                    val = sp;
            }
        }
    }

    for (auto & p : this->parts){
        p.m.x /= val;
        p.m.y /= val;
    }
    return val;
}

int MagneticSystem::fileVersion(std::string file)
{
    if (ends_with(file,".csv") || ends_with(file,".txt") || ends_with(file,".dat")){
        return 0;
    } else if (ends_with(file,".mfsys")){
        std::ifstream f(file);
        if (f.good()) {
            std::string s;
            std::getline(f,s);
            rtrim(s);
            if (s=="[header]"){
                f.close();
                return 2;
            } else {
                std::getline(f,s); //read 2 line
                std::getline(f,s); //read 3 line
                std::getline(f,s); //read 4 line
                if (s=="x\ty\tz\tMx\tMy\tMz\tr"){
                    f.close();
                    return 1;
                } else {
                    f.close();
                    return -1;
                }
            }
            f.close();
        } else {
            throw(std::string("file "+file+" not found"));
        }
    } else {
        throw(std::string("Wrong input file extention. Only mfsys, txt, dat and csv files are supported!"));
    }
    return 0;
}
