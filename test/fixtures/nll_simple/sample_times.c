/* Independent fixture queries through NonLinLoc's public grid/projection API. */
#include "GridLib.h"
int main(int argc, char **argv) {
    if (argc != 2) return 2;
    SetConstants();
    if (get_transform(0, "SIMPLE 64.0 -20.0 27.0") < 0) return 4;
    FILE *grid_file = NULL, *header_file = NULL;
    GridDesc grid = {0}; SourceDesc source = {0};
    if (OpenGrid3dFile(argv[1], &grid_file, &header_file, &grid, "time", &source, 0) < 0) return 3;
    double points[][3] = {{64.003,-19.998,0.85},{63.991,-20.015,1.7},
        {64.011,-20.023,2.1},{63.994,-19.976,2.85},{64.005,-20.01,0.4}};
    puts("latitude,longitude,depth_km,model_x_km,model_y_km,time_s");
    for (int i=0; i<5; i++) {
        double x,y; latlon2rect(0,points[i][0],points[i][1],&x,&y);
        double t = ReadAbsInterpGrid3d(grid_file,&grid,x,y,points[i][2],0);
        if (!(t >= 0.0 && t < 100.0)) return 5;
        printf("%.12f,%.12f,%.12f,%.15f,%.15f,%.15f\n",points[i][0],points[i][1],points[i][2],x,y,t);
    }
    CloseGrid3dFile(&grid,&grid_file,&header_file);
    return 0;
}
