
#include<bits/stdc++.h>
#include "de.h"
#include "time.h"

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;
int kind_function;
int g_function_number;
int g_problem_size;
unsigned int g_max_num_evaluations;
void (*gl_func)(double *, double *,int,int,int);
Fitness g_optimize_fitness;
double Fx[8][11]= {
 {0,0,0,0,0,0,0,0,0,0,0 },
 {0,100,1100,700,1900,1700,1600,2100,2200,2400,2500}, 
 {0,0,0,0,0,0,0,0,0,0,0 },
 {0,0,0,0,0,0,0,0,0,0,0 }, 
 {0,100,1100,700,1900,1700,1600,2100,2200,2400,2500},
 {0,100,1100,700,1900,1700,1600,2100,2200,2400,2500},
 {0,0,0,0,0,0,0,0,0,0,0 },
 {0,100,1100,700,1900,1700,1600,2100,2200,2400,2500}
 };
//number of runs
int num_runs = 30;

//
int g_pop_size;
double g_arc_rate;
int g_memory_size;
double g_p_best_rate;

//for output file
void print_results(string file_name);
double records[30][16];
string code;
int run_id;

int main(int argc, char **argv) {
  // rand_seed
  cout << scientific << setprecision(8);
  string code2;
  double seed[1000];
  fstream seed_file;
  seed_file.open("input_data/Rand_Seeds.txt",ios::in);
  for(int i=0;i<1000;i++) {
  	seed_file>>seed[i];
  }

    //dimension size. please select from 10, 20
  g_problem_size = 20;
  //available number of fitness evaluations
  //g_max_num_evaluations
  if(g_problem_size==10) g_max_num_evaluations=200000;
  if(g_problem_size==20)  g_max_num_evaluations=1000000;
  //MLS-LSHADE parameters
  g_pop_size = (int)round(g_problem_size*18);
  g_memory_size = 6*(g_problem_size/10)*(g_problem_size/10);
  g_arc_rate = 2.6;
  g_p_best_rate = 0.11;
  
  
  for(int kindOfFunction=7;kindOfFunction<8;kindOfFunction++) {
  	switch (kindOfFunction) {
        case 0:
            code = "000";
            code2="Basic";
            gl_func=&cec21_basic_func;
            break;

        case 1:
            code = "100";
            gl_func=&cec21_bias_func;
            code2="Bias";
            break;

        case 2:
            code = "010";
          gl_func=&cec21_shift_func;
          code2="Shift";
            break;

        case 3:
            code = "001";
           gl_func=&cec21_rot_func;
           code2="Rotation";
            break;

        case 4:
            code = "110";
            gl_func=&cec21_bias_shift_func;
            code2="BiasAndShift";
            break;

        case 5:
            code = "101";
            gl_func=&cec21_bias_rot_func;
            code2="BiasAndRotation";
            break;

        case 6:
            code = "011";
            gl_func=&cec21_shift_rot_func;
            code2="ShiftAndRotation";
            break;

        case 7:
            code = "111";
            code2="BiasShiftRotation";
           gl_func=&cec21_bias_shift_rot_func;
            break;
        default:
            cout << "Problem ID out of range." << endl;
            exit(0);
    }
  	kind_function=kindOfFunction;
  	stringstream file_name2;
    file_name2<< "tables/MLS-LSHADE_(Results_for_" <<g_problem_size<<"D_"<< code2 <<")"<<".txt";
    ofstream fout2;
    fout2.open((file_name2.str()).c_str(),ios::out );
    for (int i =0;i <10; i++) {
        g_function_number = i + 1;
        g_optimize_fitness=Fx[kind_function][g_function_number];
        cout<<"kind_function= "<<kind_function<<"........function_number= "<<g_function_number<<endl;
        cout << "\n-------------------------------------------------------" << endl;
        cout << "Function = " << g_function_number << ", Dimension size = " << g_problem_size << "\n" << endl;
        vector<Fitness> bsf_fitness_array(num_runs);
        Fitness mean_bsf_fitness = 0;
        Fitness std_bsf_fitness = 0;
        Fitness best_bsf_fitness=0;
        Fitness worst_bsf_fitness=0;
        Fitness median_bsf_fitness=0;
        for (int j = 0; j < num_runs; j++) {
            //init seed 
            run_id=j;
           	int seed_ind = (g_problem_size / 10 * g_function_number * num_runs + j) - num_runs;
            seed_ind %= 1000;
            srand(seed[seed_ind]);
                
            searchAlgorithm *alg = new MLS_LSHADE();
            bsf_fitness_array[j] = alg->run();
            cout << j+1 << "th run, " << "error value = " << bsf_fitness_array[j] << endl;
        }     
        
        // print results to file
        stringstream file_name;
        file_name << "results/MLS-LSHADE_(" << code << ")_" << g_function_number << "_" << g_problem_size << ".txt";
        print_results(file_name.str());
        sort(bsf_fitness_array.begin(),bsf_fitness_array.end());
        best_bsf_fitness=bsf_fitness_array[0];
        worst_bsf_fitness=bsf_fitness_array[num_runs-1];
        median_bsf_fitness=(bsf_fitness_array[14]+bsf_fitness_array[15])/2.0;
        
        for (int j = 0; j < num_runs; j++) mean_bsf_fitness += bsf_fitness_array[j];
        mean_bsf_fitness /= num_runs;
       
        for (int j = 0; j < num_runs; j++) std_bsf_fitness += pow((mean_bsf_fitness - bsf_fitness_array[j]), 2.0);
        std_bsf_fitness /= num_runs;
        std_bsf_fitness = sqrt(std_bsf_fitness);

        cout  << "\nmean = " << mean_bsf_fitness << ", std = " << std_bsf_fitness << endl;
        fout2<<best_bsf_fitness<<'\t'<<worst_bsf_fitness<<'\t'<<median_bsf_fitness<<'\t'<<mean_bsf_fitness<<'\t'<<std_bsf_fitness<<endl;
        bsf_fitness_array.clear();    
    }
  }
    return 0;
}
void print_results(string file_name) {
    ofstream fout;
    fout.open(file_name.c_str(),ios::out );
    if (fout.is_open()) {
        fout << scientific << setprecision(8);
        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < num_runs; j++) {
                fout << records[j][i] << " ";
            }
            fout << endl;
        }
    }
    else {
        cout << "Error! cannot onpen output file " << file_name << "." << endl;
    }

    fout.close();
}
