#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <type_traits>
#include <typeinfo>

#include <random>
#include <algorithm>

#include "bitmap_image.hpp"

#define DIM 130
#define SQ_SIZE 5

#define BLACK 0,0,0
#define WHITE 255,255,255
#define RED 255,0,0
#define GREEN 0,255,0
#define BLUE 0,0,255
#define YELLOW 255,255,0
#define GRAY 128,128,128

using namespace std;

enum State {Normal, Proliferating, Complex, Dead, Necrotic};

typedef struct Cell_t{
    enum State state = Normal;
    double ECM_value = 0;
    int ttl          = 2000;
}Cell;


typedef struct params_t {
    double r_prolif;
    double r_binding;
    double r_escape;
    double r_lysis;
    double r_decay;
    int carrying_capacity;
    double ECM_treshold;
    double ECM_degrade_rate;
}params;

random_device d;
default_random_engine rng(d());
mt19937 gen(d());
uniform_real_distribution<> dist(0, 1);

uniform_int_distribution<> time_to_die(1,8);
uniform_int_distribution<> mutation(0,1);



void draw(int x,int y,bitmap_image *image, int r, int g, int b){
    for (int i=0; i<SQ_SIZE; i++){
        for (int j=0; j<SQ_SIZE; j++){
            image->set_pixel(SQ_SIZE*x+j,SQ_SIZE*y+i,r,g,b);
        }
    }
}

void draw_field(Cell field[DIM][DIM],bitmap_image *image){
    for(int y = 0; y<DIM; y++){
        for(int x = 0; x<DIM; x++){
            if(field[y][x].state == Dead) draw(x,y,image,RED);
            //if(field[y][x].state == Normal) draw(x,y,image,BLACK);
            if(field[y][x].state == Proliferating) draw(x, y, image, GREEN);
            if(field[y][x].state == Complex) draw(x,y,image,YELLOW);
            if(field[y][x].state == Necrotic) draw(x,y,image,GRAY);
        }
    }
}

void field_setup(Cell field[DIM][DIM]){
    int a = time_to_die(gen);
    int b = mutation(gen);
    field[DIM/2][DIM/2].state=Proliferating;
    field[DIM/2][DIM/2].ttl=20+a+b;

    a = time_to_die(gen);
    b = mutation(gen);
    field[DIM/2][DIM/2-1].state=Proliferating;
    field[DIM/2][DIM/2-1].ttl = 20+a+b;

    a = time_to_die(gen);
    b = mutation(gen);
    field[DIM/2][DIM/2+1].state=Proliferating;
    field[DIM/2][DIM/2+1].ttl=20+a+b;

    a = time_to_die(gen);
    b = mutation(gen);
    field[DIM/2-1][DIM/2].state=Proliferating;
    field[DIM/2-1][DIM/2].ttl=20+a+b;

    a = time_to_die(gen);
    b = mutation(gen);
    field[DIM/2+1][DIM/2].state=Proliferating;
    field[DIM/2+1][DIM/2].ttl=20+a+b;

    uniform_real_distribution<> ecm_strength(0.8,1.2);

    for(int i =0; i<DIM; i++){
        for(int j = 0; j<DIM; j++){
            field[i][j].ECM_value=ecm_strength(gen);
        }
    }


}

//returns the number of new proliferating cells
int proliferating_process(int y, int x, Cell field[DIM][DIM], double r, params test_params){

    //If cancer cell is cabable of proliferating, check if there are normal cells in range to invade
    if (r<=test_params.r_prolif){
        vector<pair<int, int>> normal_cells_in_range;

        if(field[y][x-1].state==Normal){
            normal_cells_in_range.emplace_back(pair<int, int>{y, x-1});
        }
        if(field[y][x+1].state==Normal){
            normal_cells_in_range.emplace_back(pair<int, int>{y, x+1});
        }
        if(field[y-1][x].state==Normal){
            normal_cells_in_range.emplace_back(pair<int, int>{y-1, x});
        }
        if(field[y+1][x].state==Normal){
            normal_cells_in_range.emplace_back(pair<int, int>{y+1, x});
        }

        //No normal cells in range
        if (normal_cells_in_range.empty()){
            return 0;
        }
            //Choose one normal cell randomly
        else{
            uniform_int_distribution<> int_dist(0, normal_cells_in_range.size()-1);
            int number = int_dist(gen);
            pair<int,int> chosen_cell = normal_cells_in_range[number];

            //If cell can place the daughter cell (the ECM is weak enough)
            if(field[chosen_cell.first][chosen_cell.second].ECM_value<test_params.ECM_treshold) {

                int number = int_dist(gen);
                int a = time_to_die(gen);
                int b = mutation(gen);

                field[chosen_cell.first][chosen_cell.second].state = Proliferating;
                field[chosen_cell.first][chosen_cell.second].ttl = 20+a+b;
                return 1;
            }
            else{
                return 0;
            }
        }
    }

    //If the cancer cell cannot proliferate, it can be consumed by the immune system and form a complex
    if (r>=1-test_params.r_binding){
        field[y][x].state=Complex;
        return -1;
    }

        //Otherwise nothing happens and cancer cell remains on its place.
    else{
        return 0;
    }

}

//returns the number of new proliferating cells
int complex_process(int y, int x, Cell field[DIM][DIM], double r, params test_params){

    //If the cancer cell escapes the immune complex, it will become proliferating again
    if (r<=test_params.r_escape){
        field[y][x].state=Proliferating;
        return 1;
    }
    //If it does not escape, the cell is killed by the immune system (a programmed death called "lysis")
    if(r>=1-test_params.r_lysis){
        field[y][x].state=Dead;
        return 0;
    }
    else{
        return 0;
    }
}


void dead_process(int y, int x, Cell field[DIM][DIM], double r, params test_params){
    //The cell can decay and let a normal cell form
    if(r<=test_params.r_decay){
        field[y][x].state=Normal;
        return;
    }

    //Otherwise it can remain dead
    else{
        return;
    }
}


int process_cells(Cell field[DIM][DIM], params test_params, int* number_of_proliferating_cells){

    // Find all cells that can undergo a transformation (Proliferating, complex, dead)
    vector<pair<int, int>> cancer_cells;
    for (int y = 0; y < DIM; y++) {
        for (int x = 0; x < DIM; x++) {
            if (field[y][x].state > 0 ) {
                cancer_cells.emplace_back(pair<int, int>{y, x});
            }
        }
    }

    //shuffle their order so that they are each chosen randomly
    shuffle(cancer_cells.begin(), cancer_cells.end(), gen);


    for (const auto& ij: cancer_cells){
        int y = ij.first;
        int x = ij.second;
        double r = dist(rng);

        switch(field[y][x].state) {
            case Proliferating:
                *number_of_proliferating_cells += proliferating_process(y, x, field, r, test_params);
                break;
            case Complex:
                *number_of_proliferating_cells += complex_process(y, x, field, r, test_params);
                break;
            case Dead:
                dead_process(y, x, field, r, test_params);
                break;
            default:
                break;
        }
    }

    return 0;

}

void degrade_matrix(Cell field[DIM][DIM], params test_params){
    for (int i = 1; i< DIM-1; i++){
        for(int j = 1;j<DIM-1;j++){
            int proliferating_neighbors = 0;
            if(field[i][j].state==Normal){
                if (field[i][j+1].state==Proliferating) proliferating_neighbors++;
                if (field[i][j-1].state==Proliferating) proliferating_neighbors++;
                if (field[i+1][j].state==Proliferating) proliferating_neighbors++;
                if (field[i-1][j].state==Proliferating) proliferating_neighbors++;

                field[i][j].ECM_value -= test_params.ECM_degrade_rate*(double)proliferating_neighbors*field[i][j].ECM_value;
            }
        }
    }
}

void age_cells(Cell field[DIM][DIM], params test_params, int* number_of_number_of_proliferating_cells){
    for (int i = 1; i< DIM-1; i++){
        for(int j = 1;j<DIM-1;j++){

            if(field[i][j].state == Proliferating){
                if(field[i][j].ttl == 0){
                    field[i][j].state=Necrotic;
                    *number_of_number_of_proliferating_cells=(*number_of_number_of_proliferating_cells)-1;
                }
                else {
                    field[i][j].ttl--;
                }
            }
        }
    }
}

void apply_immunotherapy(params* test_params,params test_steps){
    test_params->r_prolif  = test_steps.r_prolif;
    test_params->r_binding  = test_steps.r_binding;
    test_params->r_escape  = test_steps.r_escape;
    test_params->r_lysis   = test_steps.r_lysis;

}

int count_complex(Cell field[DIM][DIM]){
    int count = 0;
    for (int i = 1; i< DIM-1; i++){
        for(int j = 1;j<DIM-1;j++){
            if(field[i][j].state == Complex) {
                count++;
            }
        }
    }
    return count;
}


void no_cure(params test_params){

    Cell field[DIM][DIM] = {Cell()};
    field_setup(field);

    int number_of_proliferating_cells = 5;
    double original_prolif_probability = test_params.r_prolif;

    int complex_cells;



    for (int i = 0; i<801; i++) {
        test_params.r_prolif = original_prolif_probability * (1.0 - ((double) number_of_proliferating_cells) / (double) test_params.carrying_capacity);

        //do all the processes of a time step
        process_cells(field,test_params,&number_of_proliferating_cells);

        //degrade the extra-cellular matrix
        degrade_matrix(field,test_params);

        //make all the cells age
        age_cells(field,test_params,&number_of_proliferating_cells);



        if(1) {
            complex_cells = count_complex(field);
            printf("%d;%d\n",number_of_proliferating_cells,complex_cells);
            //bitmap_image image(DIM * SQ_SIZE, DIM * SQ_SIZE);
            //draw_field(field, &image);
            //char name[40] = {0};
            //sprintf(name, "%d_obrazok.bmp", i);
            //image.save_image(name);
        }
    }
}


void cure(params test_params, params therapy_params){

    Cell field[DIM][DIM] = {Cell()};
    field_setup(field);

    int number_of_proliferating_cells = 5;
    double original_prolif_probability = test_params.r_prolif;
    params before_therapy;
    int complex_cells;
    //100 time steps of the simulation
    for (int i = 0; i<801; i++) {

        //calculate the new rate of proliferating
        if(200<i && i<=400) {
            if(i == 201){
                before_therapy = test_params;
            }
            apply_immunotherapy(&test_params, therapy_params);
        }
        else {
            if (i == 401) {
                test_params = before_therapy;
            } else {
                test_params.r_prolif = original_prolif_probability * (1.0 - ((double) number_of_proliferating_cells) /
                                                                            (double) test_params.carrying_capacity);
            }
        }


        //do all the processes of a time step
        process_cells(field,test_params,&number_of_proliferating_cells);

        //degrade the extra-cellular matrix
        degrade_matrix(field,test_params);

        //make all the cells age
        age_cells(field,test_params,&number_of_proliferating_cells);



        if(1) {
            complex_cells = count_complex(field);
            printf("%d;%d\n",number_of_proliferating_cells,complex_cells);
            //bitmap_image image(DIM * SQ_SIZE, DIM * SQ_SIZE);
            //draw_field(field, &image);
            //char name[80] = {0};
            //sprintf(name, "%d_obrazok.bmp", i);
            //image.save_image(name);
        }
    }
}


int main(int argc, char *argv[]) {
    params test_params;
    params therapy_params;

    //NO CURE
    test_params.r_prolif=0.85;
    test_params.r_binding=0.1;
    test_params.r_escape=0.5;
    test_params.r_lysis=0.35;
    test_params.r_decay=0.35;
    test_params.carrying_capacity=550;
    test_params.ECM_treshold = 0.5;
    test_params.ECM_degrade_rate=0.1;

    no_cure(test_params);
    /*----------------------------*/

    //WEAK CURE
    test_params.r_prolif=0.85;
    test_params.r_binding=0.1;
    test_params.r_escape=0.5;
    test_params.r_lysis=0.35;
    test_params.r_decay=0.35;
    test_params.carrying_capacity=550;
    test_params.ECM_treshold = 0.5;
    test_params.ECM_degrade_rate=0.1;

    therapy_params.r_lysis = 0.35;
    therapy_params.r_escape = 0.5;
    therapy_params.r_prolif = 0.2;
    therapy_params.r_binding = 0.2;

    cure(test_params, therapy_params);
    /*----------------------------*/


    //MODERATE CURE
    test_params.r_prolif=0.85;
    test_params.r_binding=0.1;
    test_params.r_escape=0.5;
    test_params.r_lysis=0.35;
    test_params.r_decay=0.35;
    test_params.carrying_capacity=550;
    test_params.ECM_treshold = 0.5;
    test_params.ECM_degrade_rate=0.1;

    therapy_params.r_lysis = 0.85;
    therapy_params.r_escape = 0.1;
    therapy_params.r_prolif = 0.58;
    therapy_params.r_binding = 0.38;

    cure(test_params, therapy_params);
    /*----------------------------*/



    //STRONG CURE
    test_params.r_prolif=0.85;
    test_params.r_binding=0.1;
    test_params.r_escape=0.5;
    test_params.r_lysis=0.35;
    test_params.r_decay=0.35;
    test_params.carrying_capacity=550;
    test_params.ECM_treshold = 0.5;
    test_params.ECM_degrade_rate=0.1;

    therapy_params.r_lysis = 0.85;
    therapy_params.r_escape = 0.1;
    therapy_params.r_prolif = 0.45;
    therapy_params.r_binding = 0.50;

    cure(test_params, therapy_params);
    return 0;
}
