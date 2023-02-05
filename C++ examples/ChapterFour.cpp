#include <fstream> // Header for file output
#include <vector> // Header for vectors
#include <cmath> // Header for sqrt(), fabs()
#include <ctime> // Header for calculating the runtime
#include <iostream> // Header for output on screen

int main()
{
// Specify simulation parameters
double Theta = -15.0;// Dimensionless applied potential
double Tmax = 1.0;// Dimensionless duration of the potential pulse
double omegaX = 1.05;
double h = 1.e-4;
double omegaT = 1.015;
double deltaT = Tmax / 4000.0;
double Ts = Tmax / 2.0;
clock_t runtime;
runtime = clock();
//Calculate other parameters
// Expanding spatial grid
double maxX = 6.0*sqrt(Tmax);
std::vector<double> X;
X.push_back(0.0);
while(X.back() < maxX)
{
X.push_back(X.back() + h);
h *= omegaX;
}
int n = X.size();// number of spacesteps
//Expanding time grid
std::vector<double> T;
T.push_back(0.0);
while(T.back() < Tmax) // where maxX = 6*sqrt(Tmax)
{
T.push_back(T.back() + deltaT);
if(T.back() > Ts) deltaT *= omegaT;
}
int m = T.size();// number of timesteps
T[m-1] = Tmax;
// Create containers
std::vector<double> alpha(n-1, 0.0);
std::vector<double> beta(n-1, 0.0);
std::vector<double> gamma(n-1, 0.0);
std::vector<double> g_mod(n-2, 0.0);
std::vector<double> delta(n-1, 1.0);
std::vector<double> d_mod(n-1, 0.0);
std::vector<double> C(n, 1.0); // concentration profile
// Calculate Thomas coefficients for T < Ts
double delT = T[1]-T[0];
std::vector<double> delX(n, 0.0);
delX[0] = X[1] - X[0];
for(int i=1; i<n-1; i++)
{
delX[i] = X[i+1] - X[i];
alpha[i] = -(2.0* delT) / ( delX[i-1] * (delX[i-1] + delX[i]) );
gamma[i] = -(2.0* delT) / ( delX[i] * (delX[i-1] + delX[i]) );
beta[i] = 1.0 - alpha[i] - gamma[i];
}
// Modified gamma coefficients for T < Ts
g_mod[0] = 0.0; // boundary condition
for(int i=1; i<n-2; i++) g_mod[i] = gamma[i]
/ (beta[i] - g_mod[i-1] * alpha[i]);
// Open file to output chronoamperogram and concentration profile
std::ofstream Chrono("Chrono_Output.txt");
std::ofstream Profile("Profile_Output.txt");
// BEGIN SIMULATION
for(int k=1; k<m; k++)
{
// Calculate Thomas coefficients for T > Ts
if(T[k] > Ts){
for(int i=1; i<n-1; i++)
{
delT = T[k] - T[k-1];
alpha[i] = -(2.0* delT) / ( delX[i-1] * (delX[i-1] + delX[i]) );
gamma[i] = -(2.0* delT) / ( delX[i] * (delX[i-1] + delX[i]) );
beta[i] = 1.0 - alpha[i] - gamma[i];
}
// Modified gamma coefficients for T > Ts
g_mod[0] = 0.0; // boundary condition
for(int i=1; i<n-2; i++) g_mod[i] = gamma[i]
/ (beta[i] - g_mod[i-1] * alpha[i]);
}
// Forward sweep - create modified deltas
d_mod[0] = 1.0 / (1.0 + exp(-Theta));
for(int i=1; i<n-1; i++) {
d_mod[i] = ( delta[i] - d_mod[i-1]*alpha[i] )
/ ( beta[i] - g_mod[i-1] * alpha[i] );
}
// Back Substitution
C[n-1] = 1.0;
for(int i=n-2; i>=0; i--) {
C[i] = d_mod[i] - g_mod[i] * C[i+1];
delta[i] = C[i];
}
//Output current (and error)
double flux = -(C[1] - C[0])/(X[1] - X[0]);
double reference = - 1.0 / sqrt(M_PI*T[k]);
double QJ = (flux - reference) / reference *100.0;
Chrono << T[k] << "\t" << flux << "\t" << QJ << "\n";
}
// END SIMULATION
runtime = (float)(clock() - runtime)/CLOCKS_PER_SEC*1e3;
std::cout << "Runtime = " << runtime << " ms" << std::endl;
//Output the concentration profile at the end of the pulse
for(int i=0; i<n; i++) Profile << X[i] << "\t" << C[i] << "\n";
}