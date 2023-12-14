#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

void Gauss(vector<vector<double>> A, vector<double> B, vector<double>& x, int n) {
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[maxRow], A[i]);
        swap(B[maxRow], B[i]);

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            B[k] -= factor * B[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        x.push_back(B[i]);
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    for (int i = 0; i < n; i++) {
        cout << "a[" << i << "] = " << x[i] << endl;
    }
}

void gragh(vector<double> X, vector<double> a, int m, int n) 
{
    ofstream dataFile("data.txt");
   
    for (int i = 0; i < n; i++)
    {
        double Y = 0;
        for (int q = 0; q <= m; q++)
        {
            Y += a[q] * pow(X[i], q);
        }
        dataFile << X[i] << " " << Y << endl;
    }
    dataFile.close();
    ofstream scriptFile("script.plt");
    scriptFile << "set term png\n";
    scriptFile << "set output 'graph.png'\n";
    scriptFile << "set multiplot\n";
    scriptFile << "plot 'data.txt' with lines\n";
    scriptFile << "plot 'oldData.txt' with lines\n";
    scriptFile << "unset multiplot";
    scriptFile.close();
    
}

int main()
{
    vector<double> X = { 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    vector<double> Y = { 0.957, 0.969, 0.976, 0.978, 0.975, 0.968, 0.954, 0.939, 0.918, 0.894};
    double n = 10;
    int m = 2;
    vector<double> powerX;
    for (int i = 1; i <= 2 * m; i++)
    {
        double sum = 0;
        for (int q = 0; q < n; q++)
        {
            sum += pow(X[q], i);
        }
        powerX.push_back(sum);
    }
    vector<vector<double>> sum_X = { {n,0,0},{0,0,0},{0,0,0} };
    for (int i = 0; i <= m; i++)
    {
        for (int q = 0; q <= m; q++)
        {
            if (i + q >= 1)
            {
                sum_X[i][q] = powerX[i + q - 1];
            }
        }
    }
    for (int i = 0; i < m + 1; i++)
    {
        for (int q = 0; q < m + 1; q++)
        {
            cout << sum_X[i][q] << " ";
        }
        cout << endl;
    }
    vector<double> praw;
    for (int i = 0; i < m + 1; i++)
    {
        double sum = 0;
        for (int q = 0; q < n; q++)
        {
            sum += Y[q] * pow(X[q], i);
        }
        praw.push_back(sum);
    }
    vector<double> a;
    Gauss(sum_X, praw, a, m + 1);
    double deviation = 0;
    for (int i = 0; i < n; i++)
    {
        double sum = Y[i];
        for (int q = 0; q <= m; q++)
        {
            sum -= a[q] * pow(X[i], q);
        }
        deviation += sum * sum / (n - m - 1);
    }
    cout << "Deviation - " << sqrt(deviation) << endl;
    gragh(X, a, m, n);
}
