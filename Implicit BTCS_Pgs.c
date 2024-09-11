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
    New_V[nodes - 1] = Vin;

    FILE *file1;
    file1 = fopen("Error_BTCS_Pgs.txt", "w");
    fprintf(file1,"Current time \t log10(ERROR)\n");
    fprintf(file1, "%.6f\t%.6f\n", Tcurr, log10(Err));
   
    while (Tcurr < Tmax) 
     {
        do 
         {
            Err = 0.0;
            for (int i = 1; i < (nodes - 1); i++) 
               {
                Old_V = New_V[i];
                New_V[i] = (1.0 / (1 + 2 * Vis)) * (Vis * New_V[i + 1] + New_V[i] + Vis * New_V[i - 1]);
                Err += pow((New_V[i] - Old_V), 2.0);
               }
            Err = sqrt(Err / (nodes - 1));
            fprintf(file1, "%.6f\t%.6f\n", Tcurr, log10(Err));
            printf("%.6f\t%d\t%.6f\n", Tcurr, iterationCount, Err);
            iterationCount++;
         } 
        while (Err> 1e-6);
        Tcurr = Tcurr + Tstep;
     }
    fclose(file1);

    FILE *file2;
    file2 = fopen("velocity_BTCS_Pgs.plt", "w");
    fprintf(file2, "VARIABLES = \"New_V\", \"Y\"\n");
    fprintf(file2, "ZONE T = \"BTCS\"\n\n");
    for (int i = 0; i < nodes; i++)
       {
        fprintf(file2, "%.6f\t%.6f\n", New_V[i], i * dy);
       }
    fclose(file2);

    return 0;
}
