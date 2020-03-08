#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <mpi.h>


//--------------------------------------------
struct _spin{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    
    _spin operator+=(const _spin& sp){
        this->x += sp.x;
        this->y += sp.y;
        this->z += sp.z;
        
        return *this;
    }

    _spin operator/=(const int& stat){
        this->x /= stat;
        this->y /= stat;
        this->z /= stat;
        
        return *this;
    }
};

struct local_data_type{
    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;
    double ex = 0.0;
    double ey = 0.0;
    double ez = 0.0;
};

struct output_data_type{
    double M = 0.0;
    double E = 0.0;
};
//--------------------------------------------



//--------------------------------------------
typedef struct _spin spin;
using Gen_rand = std::mt19937;
using int_distr = std::uniform_int_distribution<>;
using real_distr = std::uniform_real_distribution<>;
//--------------------------------------------

spin neighbourcount(const spin* const* const* const sp, const int& i, const int& j, const int& k, const int& L){
    spin res;
    res.x = 0.0;
    res.y = 0.0;
    res.z = 0.0;
    
    if(i == 0)
        res += sp[L-1][j][k];
    else
        res += sp[i-1][j][k];
    if(i == (L-1))
        res += sp[0][j][k];
    else
        res += sp[i+1][j][k];
    if(j == 0)
        res += sp[i][L-1][k];
    else
        res += sp[i][j-1][k];
    if(j == (L-1))
        res += sp[i][0][k];
    else
        res += sp[i][j+1][k];
    if(k == 0)
        res += sp[i][j][L-1];
    else
        res += sp[i][j][k-1];
    if(k == (L-1))
        res += sp[i][j][0];
    else
        res += sp[i][j][k+1];
    
    return res;
}

spin get_magnetization(const spin* const* const* const sp, const int& L){
    spin m;
    m.x = 0.0;
    m.y = 0.0;
    m.z = 0.0;
    
    for (int i = 0; i<L; ++i) {
        for (int j = 0; j<L; ++j) {
            for (int k = 0; k<L; ++k) {
                m += sp[i][j][k];
            }
        }
    }
    
    const int spin_amount = L*L*L;
    
    m /= spin_amount;
    
    return m;
}

spin get_energy(const spin* const* const* const sp, const int& L){
    spin e;
    e.x = 0.0;
    e.y = 0.0;
    e.z = 0.0;
    
    for (int i = 0; i<L; ++i) {
        for (int j = 0; j<L; ++j) {
            for (int k = 0; k<L; ++k) {
                const spin neigh = neighbourcount(sp, i, j, k, L);
                
                e.x += sp[i][j][k].x*neigh.x;
                e.y += sp[i][j][k].y*neigh.y;
                e.z += sp[i][j][k].z*neigh.z;
            }
        }
    }
    
    const int spin_amount = L*L*L;
    
    e /= (-2.0)*spin_amount;
    
    return e;
}

local_data_type get_local_data(const spin* const* const* const sp, const int& L){
    local_data_type local_data;
    
    local_data.mx = get_magnetization(sp, L).x;
    local_data.my = get_magnetization(sp, L).y;
    local_data.mz = get_magnetization(sp, L).z;
    local_data.ex = get_energy(sp, L).x;
    local_data.ey = get_energy(sp, L).y;
    local_data.ez = get_energy(sp, L).z;
    
    return local_data;
}

void output_data_initialize(output_data_type& output_data){
    
    output_data.M = 0.0;
    output_data.E = 0.0;
}

void output_data_update(output_data_type& output_data, const local_data_type& local_data){
    
    output_data.M = std::sqrt(local_data.mx*local_data.mx + local_data.my*local_data.my + local_data.mz*local_data.mz);
    output_data.E = local_data.ex + local_data.ey + local_data.ez;
}

void output_data_normalize(output_data_type& output_data, const int& stat){
    
    output_data.M /= stat;
    output_data.E /= stat;
}

void spin_lattice_update(spin***&, const int&);

void allocate_spin_lattice(spin***& sp, const int& L){
    
    sp = new spin**[L];
    for (int i = 0; i<L; ++i) {
        sp[i] = new spin*[L];
        for (int j = 0; j<L; ++j) {
            sp[i][j] = new spin[L];
        }
    }
    
    spin_lattice_update(sp, L);
}

void spin_lattice_update(spin***& sp, const int& L){
    
    for (int i = 0; i<L; ++i) {
        for (int j = 0; j<L; ++j) {
            for (int k = 0; k<L; ++k) {
                sp[i][j][k].x = 1.0;
                sp[i][j][k].y = 0.0;
                sp[i][j][k].z = 0.0;
            }
        }
    }
}

void deallocate_spin_lattice(spin***& sp, const int& L){
    
    for (int i = 0; i<L; ++i) {
        for (int j = 0; j<L; ++j) {
            delete[] sp[i][j];
        }
        delete[] sp[i];
    }
    
    delete sp;
}

int get_linear_size(const int& L_step){
    
    switch(L_step){
        case 0: {return 16;} break;
    }
    return 1;
}

int main(int argc, char** argv){
    
    int T_step_amount, L_step_amount, conf_amount, mcs_average;
    int mcs_thermalization_noinit, mcs_thermalization_init;
    double T_end, T_initial, dT;
    
    const double Pi = 3.1415926535;
    
    FILE* INPUT_FILE;
    
    INPUT_FILE = fopen("configuration.cfg", "r");
    
    fscanf(INPUT_FILE, "T_step_amount:  %d\n", &T_step_amount);
    fscanf(INPUT_FILE, "L_step_amount:  %d\n", &L_step_amount);
    fscanf(INPUT_FILE, "conf_amount:  %d\n", &conf_amount);
    fscanf(INPUT_FILE, "mcs_average:  %d\n", &mcs_average);
    fscanf(INPUT_FILE, "mcs_thermalization_noinit:  %d\n", &mcs_thermalization_noinit);
    fscanf(INPUT_FILE, "mcs_thermalization_init:  %d\n", &mcs_thermalization_init);
    fscanf(INPUT_FILE, "T_initial:  %lf\n", &T_initial);
    fscanf(INPUT_FILE, "T_end:  %lf\n", &T_end);
    fclose(INPUT_FILE);
    
    std::cout << "T_step_amount - " << T_step_amount << std::endl;
    std::cout << "T_step_amount - " << L_step_amount << std::endl;
    std::cout << "conf_amount - " << conf_amount << std::endl;
    std::cout << "mcs_average - " << mcs_average << std::endl;
    std::cout << "mcs_thermalization_noinit - " << mcs_thermalization_noinit << std::endl;
    std::cout << "mcs_thermalization_init - " << mcs_thermalization_init << std::endl;
    std::cout << "T_initial - " << T_initial << std::endl;
    std::cout << "T_end - " << T_end << std::endl;
    
    
    dT = (T_end - T_initial)/T_step_amount;
    
    Gen_rand gen_rand(std::random_device{}());
    
    output_data_type** output_data = new output_data_type*[L_step_amount];
    for (int L_step = 0; L_step<L_step_amount; ++L_step) {
        output_data[L_step] = new output_data_type[T_step_amount];
        for (int T_step = 0; T_step<T_step_amount; ++T_step) {
            output_data_initialize(output_data[L_step][T_step]);
        }
    }
    
    for (int L_step = 0; L_step<L_step_amount; ++L_step) {
        const int L = get_linear_size(L_step);
        
        spin*** sp;
        spin sp_new;
        
        allocate_spin_lattice(sp, L);
        for (int conf = 0; conf<conf_amount; ++conf) {
            std::cout << "Current conf - " << conf << std::endl;
            
            spin_lattice_update(sp, L);
            for (int T_step = 0; T_step<T_step_amount; ++T_step) {
                const double T = (T_initial + dT*T_step);
                
                const int mcs_therm = (T_step == 0 ? mcs_thermalization_init : mcs_thermalization_noinit);
                
                const int mcs_amount = mcs_therm + mcs_average;
                
                for (int mcs_step = 0; mcs_step<mcs_amount; ++mcs_step) {
                    const int spin_Amount = L*L*L;
                    
                    int_distr get_node{0, L-1};
                    real_distr get_ph{0.0, 2*Pi};
                    real_distr get_ks{0.0, Pi};
                    real_distr get_Ran{0.0, 1.0};
                    
                    auto get_node_index = [&gen_rand, &get_node](){
                        return get_node(gen_rand);
                    };
                    
                    auto get_phi_angle = [&gen_rand, &get_ph](){
                        return get_ph(gen_rand);
                    };
                    
                    auto get_ksi_angle = [&gen_rand, &get_ks](){
                        return get_ks(gen_rand);
                    };
                    
                    auto get_r = [&gen_rand, &get_Ran](){
                        return get_Ran(gen_rand);
                    };
                    
                    if(mcs_step > mcs_therm)
                        output_data_update(output_data[L_step][T_step], get_local_data(sp, L));
                    
                    for (int spin_step = 0; spin_step<spin_Amount; ++spin_step) {
                        const int i = get_node_index();
                        const int j = get_node_index();
                        const int k = get_node_index();
                        
                        const double phi = get_phi_angle();
                        const double ksi = get_ksi_angle();
                        
                        sp_new.x = std::sin(ksi)*std::cos(phi);
                        sp_new.y = std::sin(ksi)*std::sin(phi);
                        sp_new.z = std::cos(ksi);
                        
                        const spin neigh = neighbourcount(sp, i, j, k, L);
                        
                        const double dE = (sp[i][j][k].x - sp_new.x)*neigh.x +
                                          (sp[i][j][k].y - sp_new.y)*neigh.y +
                                          (sp[i][j][k].z - sp_new.z)*neigh.z;
                        
                        if(dE <= 0.0 || std::exp(-dE/T) >= get_r()){
                            sp[i][j][k].x = sp_new.x;
                            sp[i][j][k].y = sp_new.y;
                            sp[i][j][k].z = sp_new.z;
                        }
                    }
                }
            
            }
        }
        deallocate_spin_lattice(sp, L);
    }
    
    for (int L_step = 0; L_step<L_step_amount; ++L_step) {
        for (int T_step = 0; T_step<T_step_amount; ++T_step) {
            output_data_normalize(output_data[L_step][T_step], conf_amount*mcs_average);
        }
    }
    
    const int file_name_length = 100;
    char file_name[file_name_length];
    sprintf(file_name, "3d_Heis_equil_L16.dat");
    FILE* OUTPUT_FILE = fopen(file_name, "w+");
    
    for (int T_step = 0; T_step<T_step_amount; ++T_step) {
        const double T = (T_initial + dT*T_step);
        
        fprintf(OUTPUT_FILE, "%10.9lf\t", T);
        
        for (int L_step = 0; L_step<L_step_amount; ++L_step) {
            const double m = output_data[L_step][T_step].M;
            
            fprintf(OUTPUT_FILE, "\t%10.9lf", m);
        }
        
        fprintf(OUTPUT_FILE, "\t");
        
        for (int L_step = 0; L_step<L_step_amount; ++L_step) {
            const double e = output_data[L_step][T_step].E;
            
            fprintf(OUTPUT_FILE, "\t%10.9lf", e);
        }
        
        fprintf(OUTPUT_FILE, "\n");
    }
    
    fflush(OUTPUT_FILE);
    fclose(OUTPUT_FILE);
    
    
    for (int L_step = 0; L_step<L_step_amount; ++L_step) {
        delete[] output_data[L_step];
    }
    delete[] output_data;
    
    return 0;
}




































