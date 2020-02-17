#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

//------------------------------------------------------------
using spin = int;
using Gen_rand = std::mt19937;
using int_distr = std::uniform_int_distribution<>;
using real_distr = std::uniform_real_distribution<>;

struct local_data_type{
    double m = 0.0;
    double e = 0.0;
};

struct output_data_type{
    double M = 0.0;
    double M2 = 0.0;
    double E = 0.0;
    double E2 = 0.0;
};

//------------------------------------------------------------

spin neighbourcount(const spin* const* const sp, const int& i, const int& j, const int& L){
    spin res = 0;
    
    if(i == 0)
        res += sp[L-1][j];
    else
        res += sp[i-1][j];
    if(i == (L-1))
        res += sp[0][j];
    else
        res += sp[i+1][j];
    if(j == 0)
        res += sp[i][L-1];
    else
        res += sp[i][j-1];
    if(j == (L-1))
        res += sp[i][0];
    else
        res += sp[i][j+1];
        
    return res;
}

double get_magnetization(const spin* const* const sp, const int& L){
    double m = 0.0;
    
    for (int i = 0; i<L; ++i) {
        for (int j = 0; j<L; ++j) {
            m += sp[i][j];
        }
    }
    
    const int spin_amount = L*L;
    
    m /= spin_amount;
    
    return std::abs(m);
}

double get_energy(const spin* const* const sp, const int& L){
    double e = 0.0;
    
    for (int i = 0; i<L; ++i) {
        for (int j = 0; j<L; ++j) {
            e += double(sp[i][j])*neighbourcount(sp, i, j, L);
        }
    }
    
    const int spin_amount = L*L;
    
    e /= (-2.0)*spin_amount;
    
    return e;
}

local_data_type get_local_data(const spin* const* const sp, const int& L){
    local_data_type local_data;
    
    local_data.m = get_magnetization(sp, L);
    local_data.e = get_energy(sp, L);
    
    return local_data;
}

void output_data_initialize(output_data_type& output_data){
    output_data.M = 0.0;
    output_data.M2 = 0.0;
    output_data.E = 0.0;
    output_data.E2 = 0.0;
}

void output_data_update(output_data_type& output_data, const local_data_type& local_data){
    output_data.M = local_data.m;
    output_data.M2 = local_data.m*local_data.m;
    output_data.E = local_data.e;
    output_data.E2 = local_data.e*local_data.e;
}

void output_data_normalize(output_data_type& output_data, const int& stat){
    output_data.M /= stat;
    output_data.M2 /= stat;
    output_data.E /= stat;
    output_data.E2 /= stat;
}

void spin_lattice_update(spin**&, const int&);

void allocate_spin_lattice(spin**& sp, const int& L){
    sp = new spin*[L];
    for (int i = 0; i<L; ++i) {
        sp[i] = new spin[L];
    }
    
    spin_lattice_update(sp, L);
}

void spin_lattice_update(spin**& sp, const int& L){
    for (int i = 0; i<L; ++i) {
        for (int j = 0; j<L; ++j) {
            sp[i][j] = 1;
        }
    }
}

int get_linear_size(const int& L_step){
    switch(L_step){
        case 0: {return 16;} break;
    }
    return 1;
}

int main(){
    
    FILE* FILE_INPUT;
    
    int T_step_amount, conf_amount, L_step_amount;
    int mcs_thermalization_noinit, mcs_average, mcs_thermalization_init;
    double T_start, T_end, dT;
    
    FILE_INPUT = fopen("conf.cfg", "r");
    
    fscanf(FILE_INPUT, "T_step_amount:  %d\n", &T_step_amount);
    fscanf(FILE_INPUT, "L_step_amount: %d\n", &L_step_amount);
    fscanf(FILE_INPUT, "conf_amount: %d\n", &conf_amount);
    fscanf(FILE_INPUT, "mcs_thermalization_noinit: %d\n", &mcs_thermalization_noinit);
    fscanf(FILE_INPUT, "mcs_thermalization_init: %d\n", &mcs_thermalization_init);
    fscanf(FILE_INPUT, "mcs_average: %d\n", &mcs_average);
    fscanf(FILE_INPUT, "T_start: %lf\n", &T_start);
    fscanf(FILE_INPUT, "T_end: %lf\n", &T_end);
    
    dT = (T_end - T_start)/T_step_amount;
    
    Gen_rand gen_rand(std::random_device{}());
    
    output_data_type** output_data = new output_data_type*[T_step_amount];
    for (int T_step = 0; T_step<T_step_amount; ++T_step) {
        output_data[T_step] = new output_data_type[L_step_amount];
        for (int L_step = 0; L_step<L_step_amount; ++L_step) {
            output_data_initialize(output_data[T_step][L_step]);
        }
    }
    
    
    
    
    for (int L_step = 0; L_step<L_step_amount; ++L_step) {
        const int L = get_linear_size(L_step);
        
        spin** sp;
        
        allocate_spin_lattice(sp, L);
        for (int conf = 0; conf < conf_amount; ++conf) {
            std::cout << "Current conf - " << conf << std::endl;
            
            for (int T_step = 0; T_step<T_step_amount; ++T_step) {
                const int T = (T_start + dT*T_step);
                
                spin_lattice_update(sp, L);
                
                const double w[] = {std::exp(-2.0/T), std::exp(-4.0/T),
                                    std::exp(-6.0/T), std::exp(-8.0/T)};
                
                const int mcs_therm = (T_step == 0 ? mcs_thermalization_init : mcs_thermalization_noinit);
                
                
                const int mcs_amount = mcs_average + mcs_therm;
                for (int mcs_step = 0; mcs_step < mcs_amount; ++mcs_step) {
                    const int spin_amount = L*L;
                    
                    int_distr get_node_i{0, L-1};
                    real_distr get_random{0.0, 1.0};
                    
                    auto get_node_index = [&get_node_i, &gen_rand](){
                        return get_node_i(gen_rand);
                    };
                    
                    auto get_r = [&get_random, &gen_rand](){
                        return get_random(gen_rand);
                    };
                    
                    output_data_update(output_data[T_step][L_step], get_local_data(sp, L));
                    
                    for (int spin_step = 0; spin_step<spin_amount; ++spin_step) {
                        const int i = get_node_index();
                        const int j = get_node_index();
                        
                        const spin dE = 2*sp[i][j]*neighbourcount(sp, i, j, L);
                        
                        if(dE <= 0 || w[dE/2 - 1] >= get_r())
                            sp[i][j] = -sp[i][j];
                    }
                }
            }
        }
        
    }
    
    for (int T_step = 0; T_step<T_step_amount; ++T_step) {
        for (int L_step = 0; L_step<L_step_amount; ++L_step) {
            output_data_normalize(output_data[T_step][L_step], mcs_average*conf_amount);
        }
    }
    
    
    const int file_name_length = 100;
    char fname[file_name_length];
    sprintf(fname, "2d_Ising_equil_L16.dat");
    FILE* FILE_OUTPUT = fopen(fname, "w+");
    
    for (int T_step = 0; T_step<T_step_amount; ++T_step) {
        const double T = (T_start + dT*T_step);
        
        fprintf(FILE_OUTPUT, "%10.9lf\t", T);
        
        for (int L_step = 0; L_step<L_step_amount; ++L_step) {
            const double m = output_data[T_step][L_step].M;
            
            fprintf(FILE_OUTPUT, "\t%10.9lf", m);
        }
        
        fprintf(FILE_OUTPUT, "\t");
        
        for (int L_step = 0; L_step<L_step_amount; ++L_step) {
            const double e = output_data[T_step][L_step].E;
            
            fprintf(FILE_OUTPUT, "\t%10.9lf", e);
        }
        
        fprintf(FILE_OUTPUT, "\n");
    }
    
    fflush(FILE_OUTPUT);
    fclose(FILE_OUTPUT);
    
    return 0;
}

































