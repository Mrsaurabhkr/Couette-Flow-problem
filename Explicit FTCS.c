#include <stdio.h>
#include <math.h>

void main() 
{
    int nodes;
    int iterationCount = 0;
    double H = 1.0, Tstep = 0.005, Err = 0.0;
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

    double New_V[nodes], Old_V[nodes];

    for (int i = 0; i < nodes; i++) 
       {
        New_V[i] = 0.0;
       }
    New_V[nodes - 1] = Vin;

    FILE *file1;
    file1 = fopen("Error_FTCS.txt", "w");
    fprintf(file1, "%.6f\t%.6f\n", Tcurr, log10(Err));
    printf("%.6f\t%d\t%.6f\n", Tcurr, iterationCount, Err);

    while (Tcurr < Tmax) 
     {
        do 
         {
            for (int i = 0; i < nodes; i++) 
               {
                Old_V[i] = New_V[i];
               }
            Err = 0.0;
            for (int j = 1; j < (nodes - 1); j++)
              {
                New_V[j] = Old_V[j] + Vis * (Old_V[j + 1] - 2 * Old_V[j] + Old_V[j - 1]);
                Err = Err +  pow((New_V[j] - Old_V[j]), 2.0);
               }
            Err =Err + sqrt(Err /(double)(nodes - 1));
            fprintf(file1, "%.6f\t%.6f\n", Tcurr, log10(Err));
            printf("%.6f\t%d\t%.6f\n", Tcurr, iterationCount, Err);
            iterationCount++;
         } 
        while (Err> 1e-6);
        Tcurr = Tcurr + Tstep;
     }
    fclose(file1);

    FILE *file2;
    file2 = fopen("velocity_FTCS.plt", "w");
    fprintf(file2, "VARIABLES = \"New_V\", \"Y\"\n");
    fprintf(file2, "ZONE T = \"FTCS\"\n\n");
    for (int j = 0; j < nodes; j++)
       {
        fprintf(file2, "%.6f\t%.6f\n", New_V[j], j * dy);
       }
    fclose(file2);

    
}
