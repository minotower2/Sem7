#include <math.h>
#include <stdio.h>

#define ALPHA (-0.2)
#define H 0.0001
#define EPSILON 1e-12
#define EPS 1e-7


double f_x(double t, double  x, double  v) {
    return v;
}

double f_y(double t, double  x, double  v) {
    return cos(t) - (1 + ALPHA * x * x) * x;
}


void rk_step(double end_time, double x, double v, double h, double half, double sixth,
             double * res_x, double * res_v, int * err) {
    double t = 0;
    double k1_x = 0, k1_y = 0;
    double k2_x = 0, k2_y = 0;
    double k3_x = 0, k3_y = 0;
    double k4_x = 0, k4_y = 0;
    
    while(t < end_time) {
        if (fabs(x) >= 1e2 || fabs(v) >= 1e2) {
            *err = 1;
            return;
        }
        if (t + h > end_time) {
            h = end_time - t;
            half = h/2;
            sixth = h/6;
        }

        k1_x = f_x(t, x, v);
        k1_y = f_y(t, x, v);

        k2_x = f_x(t + half, x + half * k1_x, v + half * k1_y);
        k2_y = f_y(t + half, x + half * k1_x, v + half * k1_y);

        k3_x = f_x(t + half, x + half * k2_x, v + half * k2_y);
        k3_y = f_y(t + half, x + half * k2_x, v + half * k2_y);

        k4_x = f_x(t + h, x + h * k3_x, v + h * k3_y);
        k4_y = f_y(t + h, x + h * k3_x, v + h * k3_y);

        x += sixth * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
        v += sixth * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);

        t += h;
    }

    *res_x = x;
    *res_v = v;
    
}


int main(int argc, char **argv) {
    double x_init = 0, v_init = 0;
    int err = 0;
    double half = H/2;
    double sixth = H/6;
    double end_time = 2 * M_PI;
    double x = 0, v = 0;
    double delta_x = 1, delta_v = 1;
    double f_1 = 0, f_2 = 0;
    double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
    double temp_x = 0, temp_v = 0;
    double det = 0;
    
    if (argc != 3 ||
    sscanf(argv[1], "%lf", &x_init) != 1 ||
    sscanf(argv[2], "%lf", &v_init) != 1
    ) {
        printf("Usage %s: x_init v_init", argv[0]);
        return 1;        
    }


    while (fabs(delta_x) + fabs(delta_v) > EPSILON) {   
        rk_step(end_time, x_init, v_init, H, half, sixth, &x, &v, &err);
        if (err == 1) {
            printf("Bad initial values\n");
            return 1;
        }
        f_1 = x - x_init;
        f_2 = v - v_init;

        rk_step(end_time, x_init + EPS, v_init, H, half, sixth, &temp_x, &temp_v, &err);
        if (err == 1) {
            printf("Bad initial values\n");
            return 1;
        }
        J11 = (temp_x - x) / EPS - 1.0;
        J21 = (temp_v - v) / EPS;
        
        rk_step(end_time, x_init, v_init + EPS, H, half, sixth, &temp_x, &temp_v, &err);
        if (err == 1) {
            printf("Bad initial values\n");
            return 1;
        }
        J12 = (temp_x - x) / EPS;
        J22 = (temp_v - v) / EPS - 1.0;

        det = (J11 * J22 - J12 * J21);
        if (fabs(det) < 1e-14) {
            printf("Jacobian equals 0\n");
            return 2;
        }

        delta_x = (J12 * f_2 - J22 * f_1) / det;
        delta_v = (J21 * f_1 - J11 * f_2) / det;

        x_init += delta_x;
        v_init += delta_v;
        
    }

    printf("Solution: \n");
    
    printf("x_0: %.12f\n", x_init);
    printf("v_0: %.12f\n", v_init);

    printf("The difference: \n");
    printf("x: %.12f\n", f_1);
    printf("v: %.12f\n", f_2);
    

    return 0;
}
