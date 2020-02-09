#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
double ival[256][256],maxival;

int main(int argc, char **argv)
{
  	int i, j, p, q, mr, sum1, sum2;
	double lthreshold, hthreshold;
	FILE * fo1, * fo2, * fo3, * imgFile, * fopen();
	char * input;
	char temp[50];

	argc--;
	argv++;
	input = * argv;
	imgFile = fopen(input, "rb");


	if (imgFile == NULL) 
	{
		printf("Cannot open %s to read\n", input);
  		exit(0);
  	}

  	argc--;
	argv++;
	input = * argv;
	

	sprintf(temp, "%s%s", input, "mag.pgm");
	fo1 = fopen(temp, "wb");
	printf("output: %s\n", temp);


	if (fo1 == NULL) 
	{
		printf("Cannot open %s to write\n", temp);
  		exit(0);
  	}

  	sprintf(temp,"%s%s",input,"1.pgm");
	fo2=fopen(temp,"wb");
	printf("output: %s\n", temp);
	if(fo2 == NULL) {
		printf("Cannot open %s to write\n", temp);
		exit(0);
	}

  	sprintf(temp, "%s%s", input, "2.pgm");
	fo3 = fopen(temp, "wb");
	printf("output: %s\n", temp);

	if (fo3 == NULL) 
	{
  		printf("Cannot open %s to write\n", temp);
  		exit(0);
	}

	argc--;
	argv++;
	input = * argv;
	lthreshold = atof(input);

	argc--;
	argv++;
	input = * argv;
	hthreshold = atof(input);

	/* Step1: read in image file */
	for (i = 0; i < 256; i++) 
	{
	    for (j = 0; j < 256; j++) 
	    {
	      	pic[i][j] = getc(imgFile);
	      	pic[i][j] &= 0377;
	    }
	}

	/* Step2: Convolve image with maskx and masky 
	to get fx, fy at each pixel*/
	mr = 1;
	for (i = mr; i < 256 - mr; i++) 
	{
	    for (j = mr; j < 256 - mr; j++) 
	    {
	      	sum1 = 0;
	      	sum2 = 0;
	      	for (p = -mr; p <= mr; p++) 
	      	{
	        	for (q = -mr; q <= mr; q++) 
	        	{
	          		sum1 += pic[i + p][j + q] * maskx[p + mr][q + mr];
	          		sum2 += pic[i + p][j + q] * masky[p + mr][q + mr];
	        	}
	      	}
	      outpicx[i][j] = sum1;
	      outpicy[i][j] = sum2;
	    }
	}

	/* Step3: find magnitude */
	maxival = 0;
	for (i = mr; i < 256 - mr; i++) 
	{
	    for (j = mr; j < 256 - mr; j++) 
	    {
	      	ival[i][j] = sqrt((double)((outpicx[i][j] * outpicx[i][j]) + (outpicy[i][j] * outpicy[i][j])));

	      	if (ival[i][j] > maxival)
	      	{
	        	maxival = ival[i][j];
	      	}
		}
	}

	/* Step4: output result images */
	fprintf(fo1, "P5\n");
	fprintf(fo2, "P5\n");
	fprintf(fo3, "P5\n");
	fprintf(fo1, "256 256\n");
	fprintf(fo2, "256 256\n");
	fprintf(fo3, "256 256\n");
	fprintf(fo1, "255\n");
	fprintf(fo2, "255\n");
	fprintf(fo3, "255\n");

	/* output Magnitude */
	for (i = 0; i < 256; i++) 
	{
  		for (j = 0; j < 256; j++) 
  		{
    		ival[i][j] = (ival[i][j] / maxival) * 255;
    		fprintf(fo1, "%c", (char)((int)(ival[i][j])));
  		}
	}

	/* output with low threshold */
	for (i = 0; i < 256; i++) 
	{
  		for (j = 0; j < 256; j++) 
  		{
    		if (ival[i][j] > lthreshold)
    		{
    			fprintf(fo2, "%c", (char)(255));
    		}
      	
    		else
    		{
      			fprintf(fo2, "%c", (char)(0));
    		}
  		}
	}
	
	/* output with high threshold */
	for (i = 0; i < 256; i++) 
	{
  		for (j = 0; j < 256; j++) 
  		{
			if (ival[i][j] > hthreshold)
			{
  				fprintf(fo3, "%c", (char)(255));
			}

			else
			{
  				fprintf(fo3, "%c", (char)(0));
			}
  		}
	}

  fclose(fo1);
  fclose(fo2);
  fclose(fo3);

}
