#include <math.h>
#include <stdio.h>

#define ALPHA (0.2)


double f_x(double t, double  x, double  v) {
  return v;
}

double f_y(double t, double  x, double  v) {
  return cos(t) - (1 + ALPHA * x * x) * x;
}


int main(int argc, char ** argv) {
  FILE *fp;
  double x = 0;
  double v = 0;
  double end_time = 0;
  double t = 0;
  double h = 0.0001;
  double half = h/2;
  double sixth = h/6;
  double k1_x, k2_x, k3_x, k4_x;
  double k1_y, k2_y, k3_y, k4_y;
  if (argc != 4 ||
      (sscanf(argv[1], "%lf", &x)) != 1 ||
      (sscanf(argv[2], "%lf", &v)) != 1 ||
      (sscanf(argv[3], "%lf", &end_time) != 1 ||
      end_time < 0)
   ) {
    printf("Usage %s x(0) x'(0) end_time \n", argv[0]);
    return 1;
  }

  if (fabs(end_time) < 1e-14) end_time = 2*M_PI;

  fp = fopen("output.txt", "w");
  if (!fp) {
    printf("Couldn't open output.txt");
    return 2;
  }


  fprintf(fp, "%.12f %.12f\n", t, x);

  while(t < end_time) {

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

    fprintf(fp, "%.12f %.12f\n", t, x);
    t += h;
  }


  fclose (fp);
  return 0;
}
