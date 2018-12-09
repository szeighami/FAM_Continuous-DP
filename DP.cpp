#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstring>

using namespace std;

ofstream myfile;
ofstream gains;
ifstream gainsRead;
ifstream fileRead;
double* maxInData;

int memoryUsage = 0;

#ifndef WIN32
#include <sys/resource.h>
#include <sys/times.h>
#endif

struct points{
    int index;
    points* next;
};


#ifndef WIN32
void calculateExecutionTime(struct rusage *myTimeStart, struct rusage *myTimeEnd, double *userTime, double *sysTime)
{
        (*userTime) =
                ((double) (myTimeEnd->ru_utime.tv_sec  - myTimeStart->ru_utime.tv_sec)) +
                ((double) (myTimeEnd->ru_utime.tv_usec - myTimeStart->ru_utime.tv_usec)) * 1e-6;
        (*sysTime) =
                ((double) (myTimeEnd->ru_stime.tv_sec  - myTimeStart->ru_stime.tv_sec)) +
                ((double) (myTimeEnd->ru_stime.tv_usec - myTimeStart->ru_stime.tv_usec)) * 1e-6;

}
#endif



double*** opt;
double*** arr;
double** data;
double* bestPointLowerLimit;
double* bestPointUpperLimit;
int n;

double integrateExtra(int pointInS, int pointInD, double upperLimit, double w1)
{
    double coef = (data[pointInS][0]*data[pointInD][1]-data[pointInD][0]*data[pointInS][1])/(data[pointInD][1]*data[pointInD][1]);
    double extra = 0;
    double a = 2*(w1*w1)*log(w1*(data[pointInD][0]+upperLimit*data[pointInD][1])) - (w1*w1);
    double b = -2*((data[pointInD][1]*data[pointInD][1])-(data[pointInD][0]*data[pointInD][0])*(w1*w1))*log(data[pointInD][0]*w1+data[pointInD][1]) + ((data[pointInD][0]*w1)*(2*data[pointInD][1]-(data[pointInD][0]*w1)));
    double c = (data[pointInS][1]/data[pointInD][1])*(0.5*upperLimit*w1*w1-w1);
    extra = 0.25*(a - (1/(data[pointInD][0]*data[pointInD][0]))*b)+c;
    return extra;

}
double integrate(int pointInS, int pointInD, double lowerLimit, double upperLimit)
{
    if (lowerLimit >= upperLimit)
        return 0;
    double arr = 0;

    
    if (upperLimit <= 1)
    {
        double coef = (data[pointInS][0]*data[pointInD][1]-data[pointInD][0]*data[pointInS][1])/(data[pointInD][1]*data[pointInD][1]);
        arr = coef * log((data[pointInD][0]+upperLimit*data[pointInD][1])/(data[pointInD][0]+lowerLimit*data[pointInD][1]));
        arr += (data[pointInS][1]/data[pointInD][1]) * (upperLimit - lowerLimit);
        arr = arr*0.5;
    }
    if (lowerLimit >= 1)
    {
        double coef = (data[pointInS][1]*data[pointInD][0]-data[pointInD][1]*data[pointInS][0])/(data[pointInD][0]*data[pointInD][0]);
        arr = coef * log(((data[pointInD][0]/lowerLimit)+data[pointInD][1])/((data[pointInD][0]/upperLimit)+data[pointInD][1]));
        arr += (data[pointInS][0]/data[pointInD][0]) * (1/lowerLimit - 1/upperLimit);
        arr = arr*0.5;

    }
     
    if (upperLimit >1 && lowerLimit < 1)
    {
        double coef = (data[pointInS][0]*data[pointInD][1]-data[pointInD][0]*data[pointInS][1])/(data[pointInD][1]*data[pointInD][1]);
        arr = coef * log((data[pointInD][0]+1*data[pointInD][1])/(data[pointInD][0]+lowerLimit*data[pointInD][1]));
        arr += (data[pointInS][1]/data[pointInD][1]) * (1 - lowerLimit);
        arr = arr*0.5;
        
        coef = (data[pointInS][1]*data[pointInD][0]-data[pointInD][1]*data[pointInS][0])/(data[pointInD][0]*data[pointInD][0]);
        double arr2 = coef * log(((data[pointInD][0]/1)+data[pointInD][1])/((data[pointInD][0]/upperLimit)+data[pointInD][1]));
        arr2 += (data[pointInS][0]/data[pointInD][0]) * (1/1 - 1/upperLimit);
        arr2 = arr2*0.5;
        arr = arr+arr2;

    }

    
    double probIntegra = 0;
    probIntegra = 0.5*(upperLimit-lowerLimit);

    if (upperLimit > 1)
    {
        probIntegra -= 0.5*(upperLimit)-1+(1/(2*upperLimit));
    }

    if (lowerLimit > 1)
    {
        probIntegra += 0.5*(lowerLimit)-1+(1/(2*lowerLimit));
    }
    return probIntegra-arr;
}

double getAngle(int i, int j, bool& upperHalf)
{
    if (data[i][0] == data[j][0])
    {
        upperHalf = data[i][1]>data[j][1];
        return 0;
    }
    
    if (data[i][1] == data[j][1])
    {
        upperHalf = data[i][0]<data[j][0];
        return 999999;
    }
    
    upperHalf = data[i][1]>data[j][1];
    return (data[j][0]-data[i][0])/(data[i][1]-data[j][1]);

}

double calcARR(int i, int l, int u)
{
    bool upperHalf;
    double lowerLimit;
    double upperLimit;
    if (l == n)
    {
        lowerLimit = 0.000001;
    }
    else
    {
        lowerLimit = getAngle(i, l, upperHalf);
        if (!upperHalf)
            return 1;
    }
    if (u == n)
    {
        upperLimit = 999999;
    }
    else
    {
        upperLimit = getAngle(i, u, upperHalf);
        if (upperHalf)
            return 1;
    }
    double arr = 0;
    for (int p = 0; p < n; p++)
    {

        double lowerLimitForP = lowerLimit;
        if (bestPointLowerLimit[p] > lowerLimit)
            lowerLimitForP = bestPointLowerLimit[p];

        double upperLimitForP = upperLimit;
        if (bestPointUpperLimit[p] < upperLimit)
            upperLimitForP = bestPointUpperLimit[p];

        if (lowerLimitForP >= upperLimitForP)
            continue;
        double integral = integrate(i, p, lowerLimitForP, upperLimitForP);
        
        arr += integral;
    }
    return arr;
}

double getOpt(int r, int i, int z)
{
    if (opt[r][i][z] != -1)
        return opt[r][i][z];

    if (r == 0)
    {
        opt[r][i][z] = arr[i][z][n];
        return opt[r][i][z];
    }
    bool upperHalf;
    double angle_l;
    if (z == n)
        angle_l = 0;
    else
        angle_l = getAngle(i, z, upperHalf);
 
    double min = 1;
    for (int j = 0; j < n; j++)
    {
        if (j == i)
            continue;
        double angle = getAngle(i, j, upperHalf);
        if (angle < 0 && !upperHalf)
        {
            opt[r][i][z] = 1;
            return opt[r][i][z];
        }
        if (upperHalf || angle < angle_l)
            continue;
    
        double arrVal;
        arrVal = arr[i][z][j] + getOpt(r-1, j, i);

        if (arrVal < min)
            min = arrVal;
    }

    opt[r][i][z] = min;
    return opt[r][i][z];

}


int main(int argc, char** argv)
{
#ifndef WIN32
	double  userTime, sysTime;
	struct rusage myTime_start, myTime_end;
	double  userTime2, sysTime2;
	struct rusage myTime_start2, myTime_end2;
    double  userTime3, sysTime3;
	struct rusage myTime_start3, myTime_end3;

#endif

    int k = 0, d = 2;
	k = atoi(argv[2]);
	n = atoi(argv[3]);
    string filename = argv[1];
	ifstream inFile;
	inFile.open(filename.c_str());

	memoryUsage = 0;
	data = new double*[n];
	memoryUsage += sizeof(double) * n;


#ifndef WIN32
	// start time
	getrusage(RUSAGE_SELF,&myTime_start);
#endif

	for(int i = 0; i < n; i++)
	{
		data[i] = new double[d];
		memoryUsage += sizeof(double) * d;
		for(int j = 0; j < d; j++)
		{
			inFile >> data[i][j];
		}
        
	}
	inFile.close();

       
	FILE * fp;
    fp = fopen ("result.txt", "a");
	fprintf(fp, "value of k = %d\n", k);

    bool upperHalf;
    bestPointLowerLimit = new double[n];
    bestPointUpperLimit = new double[n];
    int* skylineIndex = new int[n];
    int skylineCount = 0;
    for (int i = 0; i < n; i++)
    {
        double max = 0;
        double min = 999999;
        for (int j = 0; j < n; j++)
        {
            if (i==j)
                continue;
            double angle = getAngle(i, j, upperHalf);
            if (upperHalf)
            {
                if (angle > max)
                    max = angle;
            }
            else
            {
                if (angle < min)
                    min = angle;
                if (angle <= 0)
                {
                    min = 0;
                }
            }
        }
        bestPointLowerLimit[i] = max;
        bestPointUpperLimit[i] = min;
        if (min > max)
        {
            skylineIndex[skylineCount++] = i;
        }
    }

    for (int i = 0; i < skylineCount; i++)
    {
        data[i][0] = data[skylineIndex[i]][0];
        data[i][1] = data[skylineIndex[i]][1];
        bestPointLowerLimit[i] = bestPointLowerLimit[skylineIndex[i]];
        bestPointUpperLimit[i] = bestPointUpperLimit[skylineIndex[i]];
    }

    n = skylineCount;
    arr = new double**[n];
    for (int i = 0; i < n; i++)
    {
        arr[i] = new double*[n+1];
        for (int l = 0; l < n+1; l++)
        {
            arr[i][l] = new double[n+1];
            for (int u = 0; u < n+1; u++)
            {
                arr[i][l][u] = calcARR(i, l, u);
            }
        }
    }

    opt = new double**[k];
    for (int r = 0; r < k; r++)
    {
        opt[r] = new double*[n];
        for (int i = 0; i < n; i++)
        {
            opt[r][i] = new double[n+1];
            for (int j = 0; j < n+1; j++)
            {
                opt[r][i][j] = -1;
            }
        }
    }

    double min = 1;
    for (int i = 0; i < n; i++)
    {
        if (getOpt(k-1, i, n) < min)
            min = getOpt(k-1, i, n);
    }
    cout << "ARR: " << min << endl;
#ifndef WIN32
	// end time
	getrusage(RUSAGE_SELF,&myTime_end);
	// output execution time
	calculateExecutionTime(&myTime_start, &myTime_end, &userTime, &sysTime);

	fprintf(fp, "User time : %f seconds\n", userTime);
	fprintf(fp, "System time : %f seconds\n", sysTime);
	fprintf(fp, "arr : %f\n\n", min);
#endif


    return 0;
}


