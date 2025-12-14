#ifndef COCG_H
#define COCG_H

#include <complex>
#include <vector>

class SolverCOCG {
   private:
    int n = 0;
    std::vector<int> gi;
    std::vector<int> gj;

    std::vector<std::complex<double>> di;
    std::vector<std::complex<double>> gg;
    std::vector<std::complex<double>> rp;
    std::vector<std::complex<double>> r;
    std::vector<std::complex<double>> x0;
    std::vector<std::complex<double>> z;
    std::vector<std::complex<double>> p;
    std::vector<std::complex<double>> s;

   public:
    void init(std::vector<int>& gi_s, std::vector<int>& gj_s,
              std::vector<std::complex<double>>& di_s,
              std::vector<std::complex<double>>& gg_s, int n_s);

    void solve(std::vector<std::complex<double>>& solution, double eps,
               std::vector<std::complex<double>>& rp_s, int maximum_iterations);
};

#endif