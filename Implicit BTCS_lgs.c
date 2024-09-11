#include <stdio.h>
#include <math.h>

int main() 
{
    int nodes;
    int iterationCount = 0;
    double H = 1.0, Tstep = 0.01, Err = 0.0;
    double Tcurr = 0.0, Tmax, Tstepzone;
    double Vin = 1.0, Reyn = 100.0;

    printf("Enter the number of nodes: \n");
    scanf("%d", &nodes);
    printf("Enter the time to converge: \n");
    scanf("%lf", &Tmax);
    Tstepzone = Tcurr/Tstep;

    double dy = H / (nodes - 1);
    double Gamma = (1.0 / Reyn);
    double Vis = Gamma * (Tstep / (dy * dy));

    double New_V[nodes], Old_V;

    for (int i = 0; i < nodes; i++) 
       {
        New_V[i] = 0.0;
       }
    New_V[nodes - 1] = 1.0;

    double di[nodes - 2], pi[nodes - 2], qi[nodes - 2];
    double Nj = Vis, Pj = -(1 + 2 * Vis), Sj = Vis;
    qi[nodes - 1] = 1.0;

    FILE *file1;
    file1 = fopen("Error_btcs_lgauss_seidel.txt", "w");
    fprintf(file1, "%.6f\t%.6f\n", Tcurr, log10(Err));
    printf("%.6f\t%d\t%.6f\n", Tcurr, iterationCount, Err);

    while (Tcurr < Tmax) 
    {
        do
          {
            Err = 0.0;
            for (int i = 1; i < (nodes - 1); i++) 
             {
                di[i] = -New_V[i];
                pi[0] = -Nj / Pj;
                qi[0] = di[i] / Pj;
                pi[i] = -Nj / (Pj + Sj * pi[i - 1]);
                qi[i] = (di[i] - Sj * qi[i - 1]) / (Pj + Sj * pi[i - 1]);
            }

            for (int i = (nodes - 2); i > 0; i--)
            {
                Old_V = New_V[i];
                New_V[i] = pi[i] * New_V[i + 1] + qi[i];
                Err += pow((New_V[i] - Old_V), 2.0);
            }

            Err = sqrt(Err / (nodes - 1));
            fprintf(file1, "%.6f\t%.6f\n", Tcurr, log10(Err));
            printf("%.6f\t%d\t%.6f\n", Tcurr, iterationCount, Err);
            iterationCount++;
        } while (Err > 1e-6);
        Tcurr = Tcurr + Tstep;
    }
    fclose(file1);

    FILE *file2;
    file2 = fopen("velocity_lgauss_seidel.plt", "w");
    fprintf(file2, "VARIABLES = \"u\", \"Y\"\n");
    fprintf(file2, "ZONE T = \"BTCS\"\n\n");
    for (int j = 0; j < nodes; j++) {
        fprintf(file2, "%.6f\t%.6f\n",New_V[j], j * dy);
    }
    fclose(file2);

    return 0;
}
