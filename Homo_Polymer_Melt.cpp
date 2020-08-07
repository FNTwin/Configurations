#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <random>
#include <iomanip>
//------------Parameters------------------
constexpr int idN=100; //simulation frames
constexpr int N_type=1, B_type=1; //Number of particles and bonds
constexpr int N_particles=3000; //Number of CG-DPD beads
constexpr int N_chained=10; // Number of beads in a chain
constexpr int calc_N_chains(const int num_p,const int num_c) {
	return num_p/num_c;
}
constexpr int N_chains=calc_N_chains(N_particles,N_chained); // Number of chains in the system

//--------------Property----------------
constexpr float Lx=10.f,Ly=10.f,Lz=10.f; //Dimension of the box
constexpr float mass=1.f; //Masses of the particles
constexpr float calc_density(const float volume, const int num_p){
    return (float)num_p/volume;
}
constexpr float vol=Lx*Ly*Lz; //Volume
constexpr float rho=calc_density(vol,N_particles); //Density of the system

//-------------Positions----------------
using Vec_rcom=std::array<std::array<float,3>,N_chains>;
using Vec_chain=std::array<std::array<float,3>,N_chained>;
using Vec_system=std::array<std::array<float,3>,N_particles>;
using Vec=std::array<float,3>;

Vec_rcom rcom;
Vec_chain rchain;
Vec_system rsystem;

//----------Flags-----------------
Vec_system flags;

//---------Functions declarations----------
void InitCOMS(int N, Vec_rcom& rpos) noexcept ;
void InitChain(int num_chains,Vec_chain& rpos) noexcept ;
void InitConf(int N, int N_bc, int N_c,
        Vec_rcom& r_com, Vec_chain& r_chain, Vec_system r_sys,
        float L_x,float L_y,float L_z, Vec_system flag) noexcept ;
void writefile(int N,int N_bc, int N_c,
               float L_x,float L_y, float L_z, float m,
               int type_N, int type_B,
               Vec_system& flag, Vec_system &r_system) noexcept ;

//----------------------------------------
//----------------Main-------------------
//----------------------------------------

int main(){
    std::cout<<"rho : "<<rho<<std::endl;
    std::cout<<"------------------------------------------------------------------\n";
    std::cout<<"|-------------------Initial Position of coms---------------------|\n";
    InitCOMS(N_particles,rcom);
    std::cout<<"|--------------------------Generated-----------------------------|\n";
    std::cout<<"|--------------------Initial Chains Generated--------------------|\n";
    InitChain(N_chains,rchain);
    std::cout<<"|--------------------------Generated-----------------------------|\n";
    std::cout<<"|-----------------Initial Configuration Generated----------------|\n";
    InitConf(N_particles,N_chained,N_chains,rcom,rchain,rsystem,Lx,Ly,Lz,flags);
    std::cout<<"|--------------------------Generated-----------------------------|\n";
    std::cout<<"|---------------------Writing HomoMelt.conf----------------------|\n";
    writefile(N_particles,N_chained,N_chains,Lx,Ly,Lz,mass,N_type,B_type,flags,rsystem);
    std::cout<<"------------------------------------------------------------------\n";
    return 0;
}

//---------------Functions Implementation----------------------------------
void InitCOMS(int N, Vec_rcom& rpos) noexcept {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<> distr(-10, 10); // define the range

    for(int n=0; n<N; ++n){
        for(int i=0; i<3; ++i) {
            rpos[n][i] = distr(eng);
        }
    }// generate numbers
}

void InitChain(int num_chains, Vec_chain& rpos) noexcept {
    //Inizialize variable
    Vec gamma,r_com;
    float r_cutoff=3.f;
    rpos[0][0]=1.f;
    rpos[0][1]=1.f;
    rpos[0][2]=1.f;

    r_com[0]+=rpos[0][0];
    r_com[1]+=rpos[0][1];
    r_com[2]+=rpos[0][2];

    //Inizialize random
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<> distr(-1, 1);

    //Generate positions
    for (int i = 1; i < num_chains ; ++i) {
        rpos[i][0] = rpos[i-1][0] + (2. * distr(eng) - 1) * r_cutoff;
        rpos[i][1] = rpos[i-1][1] + (2. * distr(eng) - 1) * r_cutoff;
        rpos[i][2] = rpos[i-1][2] + (2. * distr(eng) - 1) * r_cutoff;

        r_com[0]+=rpos[i][0];
        r_com[1]+=rpos[i][1];
        r_com[2]+=rpos[i][2];
    }

    //Normalize and rescale
    for (auto& x: r_com) {
        x/=num_chains;
    }
    for (int j = 0; j < num_chains; ++j) {
        for (int k = 0; k < 3; ++k) {
            rpos[j][k] -= r_com[k];
        }
    }
}

void InitConf(int N, int N_bc, int N_c,
        Vec_rcom& r_com, Vec_chain& r_chain, Vec_system r_sys,
        float L_x,float L_y,float L_z, Vec_system flag) noexcept {
    float Lxp=L_x/2.f , Lyp=L_y/2.f, Lzp=L_z/2.f;
    int iter=0;

    //Copy rchain to rcom with Periodic Boundary Conditions
    for (int i = 0; i < N_c ; ++i) {
        for (int j = 0; j < N_bc; ++j) {
            r_sys[iter][0] = r_com[i][0] + r_chain[j][0];
            r_sys[iter][1] = r_com[i][1] + r_chain[j][1];
            r_sys[iter][2] = r_com[i][2] + r_chain[j][2];

            //Shift system
            r_sys[iter][0] -= Lxp;
            r_sys[iter][1] -= Lyp;
            r_sys[iter][2] -= Lzp;

            //Apply PBC and store flags
            flag[iter][0] = r_sys[iter][0] > Lxp ? 1 : -1;
            flag[iter][1] = r_sys[iter][1] > Lyp ? 1 : -1;
            flag[iter][2] = r_sys[iter][2] > Lzp ? 1 : -1;

            //Normalize
            r_sys[iter][0] -= L_x *std::roundf(r_sys[iter][0] /L_x);
            r_sys[iter][1] -= L_y *std::roundf(r_sys[iter][1] /L_y);
            r_sys[iter][2] -= L_z *std::roundf(r_sys[iter][2] /L_z);

            iter++;
        }
    }
}

void writefile(int N,int N_bc, int N_c,
        float L_x,float L_y, float L_z, float m,
        int type_N, int type_B,
        Vec_system& flag, Vec_system& r_system) noexcept {

    float Lxp=L_x/2.f , Lyp=L_y/2.f, Lzp=L_z/2.f;
    std::ofstream file;
    file.open( "HomoMelt.conf");
    file<<std::fixed<<std::setprecision(6);
    //Box writing
    file<< "\t" << -Lxp << "\t" << Lxp << " " << " xlo xhi\n";
    file<< "\t" << -Lyp << "\t" << Lyp << " " << " ylo yhi\n";
    file<< "\t" << -Lzp << "\t" << Lzp << " " << " zlo zhi\n" <<std::endl;
    //Atoms & bonds
    file<< "\t\t" <<N<< " atoms\n" <<std::endl;
    file<< "\t\t" <<N_c*(N_bc-1)<<" bonds\n"<<std::endl;
    file<<"\t\t\t" <<type_N<< " atom types\n" <<std::endl;
    file<<"\t\t\t" <<type_B<< " bond types\n" <<std::endl;
    //If you have more than 1 atoms you should do a loop for writing the masses
    file<<"Masses\n" <<std::endl;
    file<<"\t\t\t" <<type_N<<"\t"<<std::fixed<<std::setprecision(6)<<m<<"\n" <<std::endl;
    //Positions of Atoms
    file<<"Atoms\n" <<std::endl;
    int molID=0,iter=0;
    for(int i=0; i<N_c;i++){
        for (int j = 0; j < N_bc; ++j) {
            file<<"\t\t\t"<<iter<<"\t" <<type_N<<"\t"<<molID<<"\t";
            file<<r_system[iter][0]<<"\t" <<r_system[iter][1]<<"\t" <<r_system[iter][2]<<"\t";
            file<<flag[iter][0]<<"\t"<<flag[iter][1]<<"\t"<<flag[iter][2]<<"\n";
            iter++;
        }
        molID++;
    }
    //Bonds
    std::cout<<"\nBonds\n"<<std::endl;
    iter=0;
    int iterID=1;
    for(int i=0; i<N_c;i++){
        for (int j = 0; j < N_bc -1; ++j) {
            file<<"\t\t\t"<<iterID<<"\t"<<1<<"\t"<<iter<<"\t"<<iter+1<<"\n";
            iter++;
            iterID++;
        }
        iter++;
    }
    file.close();

}

