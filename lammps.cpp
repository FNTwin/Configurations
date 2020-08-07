#include <iostream>
#include <array>
#include <math.h>
#include <random>
#include <fstream>

constexpr int idN=100;

constexpr float Lx=2.f;
constexpr float Ly=2.f;
constexpr float Lz=2.f;

//Particle spec
constexpr int Ntype=1;
constexpr int N=300000;
constexpr float masses=1.0;
using vec=std::array<std::array<float,3>,N>;
vec rpos;

constexpr float V=Lx*Ly*Lz;
constexpr float rho=N/V;
float Ncheck=floor(V*rho);

void InitConf(const int Num, vec& r_pos) noexcept {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<> distr(-10, 10); // define the range

    for(int n=0; n<Num; ++n){
        for(int i=0; i<3; ++i) {
            r_pos[n][i] = distr(eng);
        }
    }// generate numbers
}
void ulozLAMMPS(const int Num,const float L_x,const float L_y, const float L_z, const int type, const float m, const vec& r_pos) noexcept {
    std::ofstream file;
    file.open( "configuration.conf");
    float Lxp=L_x/2, Lyp=L_y/2, Lzp=L_z/2;
    file<<std::fixed<<std::setprecision(6);
    file<< "\t" << -Lxp << "\t" << Lxp << " " << " xlo xhi\n";
    file<< "\t" << -Lyp << "\t" << Lyp << " " << " ylo yhi\n";
    file<< "\t" << -Lzp << "\t" << Lzp << " " << " zlo zhi\n" <<std::endl;
    file<< "\t\t" <<Num<< " atoms\n" <<std::endl;
    file<<"\t\t\t" <<type<< " atom types\n" <<std::endl;
    file<<"Masses\n" <<std::endl;
    file<<"\t\t\t" <<type<<"\t"<<std::fixed<<std::setprecision(6)<<m<<"\n" <<std::endl;
    file<<"Atoms\n" <<std::endl;
    for(int i=1; i<=Num;i++){
        file<<"\t\t\t"<<i<<"\t" <<type<<"\t"<<r_pos[i-1][0]<<"\t" <<r_pos[i-1][1]<<"\t" <<r_pos[i-1][2]<<"\n";
    }
    file.close();
}
int main() {
    if (N == Ncheck){
        std::cout<<"rho = " << rho<<" \n";
        std::cout<<"-----------------------------------------------\n";
        std::cout<<"|----------Generating configurations----------|\n";
        InitConf(N,rpos);
        std::cout<<"|-----------Configuration generated-----------|\n";
        std::cout<<"|------------Generating .conf file------------|\n";
        ulozLAMMPS(N,Lx,Ly,Lz,Ntype,masses, rpos);
        std::cout<<"|-------------.conf file generated------------|\n";
        std::cout<<"-----------------------------------------------\n";
    }
    else{
        std::cout<< "System specification seems to be wrong";
    }
    return 0;

}
