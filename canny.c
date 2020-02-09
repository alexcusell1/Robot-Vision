#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


void     print(double** matrix, int rows, int cols);
double** getMemory(int rows, int cols);
void     freeMemory(double** mem, int rows, int cols);
void     printMatrix(char* title, double** matrix, int rows, int cols, int padding);
void     printpicMatrix(char* title, int** matrix, int rows, int cols, int padding);
void     writeFileInt(int** pic, char* fileName, int numRows, int numCols);
void     writeFileDouble(double** pic, char* fileName, int numRows, int numCols);

int row, col;
int i, j;
int lthreshold, hthreshold;
char * input;
char temp[50];

int main(int argc, char **argv) 
{
    int      numRows = 256, numCols = 256;
    int      gausWidth, gausCenter;
    int      picMaxValue = 255;
    int      row, col;
    int      toDo = 1;
    int**    pic;
    double** gaus;
    double** maskx;
    double** masky;
    double** xgaus;
    double** ygaus;
    double** pic_mag_gaus;
    int**    peaks;
    int**    final_pic;
    double   sigma;
    double   maxMagVal = 0;
    double   percent;
    double   val = 0;

    //this takes the input file from the user
    FILE* imgFile = fopen(argv[1], "rb");
    //this takes the input value for sigma from the user
    sigma = atof(argv[2]);
    //this takes the input value for percentage from the user
    percent = atof(argv[3]) / 100;

    fgets(temp, 50, imgFile);
    fgets(temp, 50, imgFile);
    fgets(temp, 50, imgFile);

    if ( !( (temp[0]=='2') && (temp[1]=='5') && (temp[2]=='5')))
    {
        fgets(temp, 50, imgFile);
    }

    //this decides the dimensions of the matrices
    gausWidth = 1 + (6 * sigma);
    gausCenter = (gausWidth / 2);

    pic = calloc(numRows, sizeof(int*));
    final_pic = calloc(numRows, sizeof(int*));
    peaks = calloc(numRows, sizeof(int*));

    for(i = 0; i < numRows; i++) 
    {
        pic[i] = calloc(numCols, sizeof(int));
        final_pic[i] = calloc(numCols, sizeof(int*));
        peaks[i] = calloc(numCols, sizeof(int));
    }

    gaus = getMemory(gausWidth, gausWidth);
    maskx = getMemory(gausWidth, gausWidth);
    masky = getMemory(gausWidth, gausWidth);
    xgaus = getMemory(numRows, numCols);
    ygaus = getMemory(numRows, numCols);
    pic_mag_gaus = getMemory(numRows, numCols);

    for(i = 0; i < numRows; i++) 
    {
        for(j = 0; j < numCols; j++) 
        {
            pic[i][j] = getc(imgFile);
            pic[i][j] = pic[i][j] & 0377;
        }
    }
    fclose(imgFile);

    for(row = -gausCenter; row <= gausCenter; row++) 
    {
        for(col = -gausCenter; col <= gausCenter; col++) 
        {
            val = exp(-((col * col) + (row * row)) / (2 * sigma * sigma));
            gaus[row + gausCenter][col + gausCenter] = val;
            maskx[row + gausCenter][col + gausCenter] = col * val;
            masky[row + gausCenter][col + gausCenter] = row * val;
        }
    }

    //convolution
    for(row = gausCenter; row < numRows - gausCenter; row++) 
    {
        for(col = gausCenter; col < numCols - gausCenter; col++) 
        {
            double verticalSum = 0, horizontalSum = 0;
            int row_c, col_c;
            for(row_c = -gausCenter; row_c <= gausCenter; row_c++) 
            {
                for(col_c = -gausCenter; col_c <= gausCenter; col_c++) 
                {
                    verticalSum += (((double)(pic[row + row_c][col + col_c])) * masky[row_c + gausCenter] [col_c + gausCenter]);
                    horizontalSum += (((double)(pic[row + row_c][col + col_c])) * maskx[row_c + gausCenter] [col_c + gausCenter]);
                }
            }

            xgaus[row][col] = horizontalSum;
            ygaus[row][col] = verticalSum;
        }
    }

    //this finds the magnitude
    for(row = 0; row < numRows; row++) 
    {
        for(col = 0; col < numCols; col++) 
        {
            pic_mag_gaus[row][col] = sqrt((xgaus[row][col] * xgaus[row][col]) + (ygaus[row][col] * ygaus[row][col]));
            if(maxMagVal < pic_mag_gaus[row][col]) 
            {
                maxMagVal = pic_mag_gaus[row][col];
            }
        }
    }

    //scale
    for(row = 0; row < numRows; row++) 
    {
        for(col = 0; col < numCols; col++) 
        {
            pic_mag_gaus[row][col] = (pic_mag_gaus[row][col] / maxMagVal) * 255;
        }
    }

    for(row = gausCenter; row < numRows - gausCenter; row++) 
    {
        for(col = gausCenter; col < numCols - gausCenter; col++) 
        {

            if(xgaus[row][col] == 0.0) 
            {
               xgaus[row][col] = .00001;
            }

            double slope = ygaus[row][col] / xgaus[row][col];

            if((-.4142 < slope) && (slope <= .4142)) 
            {
                if((pic_mag_gaus[row][col - 1] < pic_mag_gaus[row][col]) && (pic_mag_gaus[row][col + 1] < pic_mag_gaus[row][col])) 
                {
                    peaks[row][col] = 255;
                }
            } 

            else if((.4142 < slope) && (slope <= 2.4142)) 
            {
                if((pic_mag_gaus[row - 1][col - 1] < pic_mag_gaus[row][col]) && (pic_mag_gaus[row + 1][col + 1] < pic_mag_gaus[row][col])) 
                {
                    peaks[row][col] = 255;
                }
            } 

            else if((-2.4142 < slope) && (slope <= -.4142)) 
            {
                if((pic_mag_gaus[row + 1][col - 1] < pic_mag_gaus[row][col]) && (pic_mag_gaus[row - 1][col + 1] < pic_mag_gaus[row][col])) 
                {
                    peaks[row][col] = 255;
                }
            } 

            else 
            {
                if((pic_mag_gaus[row - 1][col] < pic_mag_gaus[row][col]) && (pic_mag_gaus[row + 1][col] < pic_mag_gaus[row][col])) 
                {
                    peaks[row][col] = 255;
                }
            }
        }
    }

    writeFileInt(peaks, "cannypeaks.pgm", numRows, numCols);

    //this makes the histogram
    int* histogram = calloc(256, sizeof(int)); // represent values 0 -> 255, zeros out indexes

    for(row = 0; row < numRows; row++) 
    {
        for(col = 0; col < numCols; col++) 
        {
            histogram[(int)(pic_mag_gaus[row][col])]++;
        }
    }

    //mark edges
    int cutOff = percent * (numRows) * (numCols);
    int area = 0;
    for(row = 255; (0 < row) && (area <= cutOff); row--) 
    {
        area += histogram[row];
    }

    //this assigns the thresholds
    hthreshold = row;
    lthreshold = (int)(.35 * hthreshold);
    printf("The high threshold = %d\n", hthreshold);
    printf("The low threshold = %d\n", lthreshold);

    //double thresholding
    for(row = 0; row < numRows; row++) 
    {
        for(col = 0; col < numCols; col++) 
        {
            if(peaks[row][col] == 255) 
            {
                if(hthreshold <= pic_mag_gaus[row][col]) 
                {
                    peaks[row][col] = 0;
                    final_pic[row][col] = 255;
                } 

                else if(pic_mag_gaus[row][col] < lthreshold) 
                {
                    peaks[row][col] = 0;
                    final_pic[row][col] = 0;
                }
            }
        }
    }

    while(toDo == 1) 
    {
        toDo = 0;
        for(row = 0; row < numRows; row++) 
        {
            for(col = 0; col < numCols; col++) 
            {
                if(peaks[row][col] == 255) 
                {
                    for(i = -1; i <= 1; i++) 
                    {
                        for(j = -1; j <= 1; j++) 
                        {
                            if(peaks[row + i][col + j] == 255) 
                            {
                                peaks[row][col] = 0;
                                final_pic[row][col] = 255;
                                toDo = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    writeFileInt(final_pic, "cannyfinal.pgm", numRows, numCols);
    writeFileDouble(pic_mag_gaus, "cannymag.pgm", numRows, numCols);

    //FREE EVERYTHING!!!
    for(i = 0; i < numRows; i++) 
    {
        free(pic[i]);
        free(final_pic[i]);
        free(peaks[i]);
    }

    free(pic);
    free(histogram);
    free(final_pic);
    free(peaks);
    freeMemory(gaus, gausWidth, gausWidth);
    freeMemory(maskx, gausWidth, gausWidth);
    freeMemory(masky, gausWidth, gausWidth);
    freeMemory(xgaus, numRows, numCols);
    freeMemory(ygaus, numRows, numCols);
    freeMemory(pic_mag_gaus, numRows, numCols);

    return 0;
}

double** getMemory(int rows, int cols) 
{
    double** tempMem = calloc(rows, sizeof(double*));
    for(i = 0; i < rows; i++) 
    {
        tempMem[i] = calloc(cols, sizeof(double));
    }
    return tempMem;
}

void freeMemory(double** mem, int rows, int cols) 
{
    for(i = 0; i < rows; i++) 
    {
        free(mem[i]);
    }
    free(mem);
}

void printMatrix(char* title, double** matrix, int rows, int cols, int padding) 
{
    printf("\n\n*************| %s |*************\n\n", title);
    for(i = padding; i < rows - padding; i++) 
    {
        for(j = padding; j < cols - padding; j++) 
        {
            printf("| %10lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void printpicMatrix(char* title, int** matrix, int rows, int cols, int padding) 
{
    printf("\n\n*************| %s |*************\n\n", title);
    for(i = padding; i < rows - padding; i++) 
    {
        for(j = padding; j < cols - padding; j++) 
        {
            printf("| %10d ", matrix[i][j]);
        }
        printf("\n");
    }
}

void writeFileInt(int** pic, char* fileName, int numRows, int numCols) 
{
    FILE* imgFile = fopen(fileName, "wb");
    int row, col;

    fprintf(imgFile, "%s\n","P5");
    fprintf(imgFile, "%d %d\n%d\n", numRows, numCols, 255);

    for(row = 0; row < numRows; row++) 
    {
        for(col = 0; col < numCols; col++) 
        {
            fprintf(imgFile, "%c", (char) (pic[row][col]));
        }
    }
    fclose(imgFile);
}

void writeFileDouble(double** pic, char* fileName, int numRows, int numCols) 
{
    FILE* imgFile = fopen(fileName, "wb");

    fprintf(imgFile, "%s\n","P5");
    fprintf(imgFile, "%d %d\n%d\n", numRows, numCols, 255);
    
    for(row = 0; row < numRows; row++) 
    {
        for(col = 0; col < numCols; col++) 
        {
            fprintf(imgFile, "%c", (char)((int) (pic[row][col])));
        }
    }
    fclose(imgFile);
}