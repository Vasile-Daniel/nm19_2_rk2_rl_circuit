#include <stdio.h>
#include <math.h>
#define R         10.0    // resistance R [Ohms]
#define L         0.01     // inductance L [H]
#define E         5.0      // constant input voltage E [V] 
#define DT        0.00001    // time step
#define NUM_STEPS 500    // number of steps

// Function to calculate the derivative (rate of change of current through the inductor)
double derivative(double i) {
    double tau = L / R;
    return (E - R * i) / L; // differential equation of the RL circuit
}

// Runge-Kutta method of second order for solving the differential equation
void runge_kutta(double t[], double i[]) {
    for (int j = 0; j < NUM_STEPS - 1; j++) {
        double k1 = DT * derivative(i[j]);
        double k2 = DT * derivative(i[j] + k1 / 2.0);
        i[j + 1] = i[j] + k2; // update current for the next step
        t[j + 1] = t[j] + DT; // update time for the next step
    }
}

int main() {
    double t[NUM_STEPS]; // array for time
    double i[NUM_STEPS]; // array for current

    t[0] = 0.0;  // initial time
    i[0] = 0.0;  // initial current through the inductor

    printf("Applying the Runge-Kutta method to the RL circuit with constant input voltage E = %.2f V\n", E);

    FILE *fp = NULL; 
    fp = fopen("rl_rk2_output_data.txt", "w");

    if(fp == NULL){
        printf("Error opening the file");
        return 1;
    }

    runge_kutta(t, i);

    // Displaying the values of time and current
    for (int j = 0; j < NUM_STEPS; j++) {
        printf("Time: %.2f, Current: %.4f A\n", t[j], i[j]);
        fprintf(fp, "%lf\t %lf\n", t[j], i[j]);
    }

    // ================= GNU PLOT =========================== // 

    FILE *gnuplotPipe = NULL; 
    gnuplotPipe = popen("gnuplot -persist", "w");

    if (gnuplotPipe == NULL){
        printf("Error opening Gnuplot");
        return 1;
    }

    fprintf(gnuplotPipe, "set terminal png\n");
    fprintf(gnuplotPipe, "set output 'rl_rk2_plot.png'\n");
    fprintf(gnuplotPipe, "set grid\n");
    fprintf(gnuplotPipe, "set title 'Behaviour of current through the inductor'\n");
    fprintf(gnuplotPipe, "set xlabel 't [s]'\n");
    fprintf(gnuplotPipe, "set ylabel 'i_L [A]'\n");
    fprintf(gnuplotPipe, "set label 'RK2 method' at 3, 3\n");

    fprintf(gnuplotPipe, "plot 'rl_rk2_output_data.txt' using 1:2 with lines linewidth 2 linecolor rgb 'blue'\n");

    fflush(gnuplotPipe);
    fprintf(gnuplotPipe, "exit\n");
    pclose(gnuplotPipe);

    return 0;
}
