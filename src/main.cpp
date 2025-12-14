#include <stdint.h>

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "../include/slae.h"
using namespace std;
int main() {
    const std::string folder_index = "1";

    ComplexSlae complex_slae;

    complex_slae.inputSLAEData(folder_index);
    complex_slae.solveSLAE();

    // FILE* fp;
    // string folder = "data/complex/1";
    // int _n, _maxiter;
    // size_t n, maxiter;
    // double eps;

    // fp = fopen(string(folder + "/kuslau").c_str(), "r");
    // fscanf(fp, "%d%lf%d", &_n, &eps, &_maxiter);
    // fclose(fp);
    // maxiter = (size_t)_maxiter;
    // n = (size_t)(_n / 2);

    // uint32_t a;
    // double b, c;

    // vector<size_t> idi;
    // idi.resize(n + 1);
    // fp = fopen((folder + "/idi").c_str(), "rb");
    // for (size_t i = 0; i < idi.size(); i++) {
    //     size_t rdb = fread(&a, sizeof(a), 1, fp);
    //     if (rdb == 0) throw;
    //     idi[i] = --a;
    // }
    // fclose(fp);

    // complex<double>* di = new complex<double>[n];
    // fp = fopen((folder + "/di").c_str(), "rb");
    // for (size_t i = 0, j = 0; i < idi[n]; i++, j++) {
    //     size_t rdb = fread(&b, sizeof(b), 1, fp);
    //     if (rdb == 0) throw;
    //     if (idi[j + 1] - idi[j] == 2) {
    //         rdb = fread(&c, sizeof(c), 1, fp);
    //         if (rdb == 0) throw;
    //         di[j].imag(c);
    //         i++;
    //     }
    //     di[j].real(b);
    // }
    // fclose(fp);

    // size_t ig_size = n + 1;
    // size_t* ig = new size_t[ig_size];
    // fp = fopen((folder + "/ig").c_str(), "rb");
    // for (size_t i = 0; i < ig_size; i++) {
    //     size_t rdb = fread(&a, sizeof(a), 1, fp);
    //     if (rdb == 0) throw;
    //     ig[i] = --a;
    // }
    // fclose(fp);

    // size_t jg_size = ig[n];
    // size_t* jg = new size_t[jg_size];
    // fp = fopen((folder + "/jg").c_str(), "rb");
    // for (size_t i = 0; i < jg_size; i++) {
    //     size_t rdb = fread(&a, sizeof(a), 1, fp);
    //     if (rdb == 0) throw;
    //     jg[i] = --a;
    // }
    // fclose(fp);

    // vector<size_t> ijg;
    // ijg.resize(ig[n] + 1);
    // fp = fopen((folder + "/ijg").c_str(), "rb");
    // for (size_t i = 0; i < ijg.size(); i++) {
    //     size_t rdb = fread(&a, sizeof(a), 1, fp);
    //     if (rdb == 0) throw;
    //     ijg[i] = --a;
    // }
    // fclose(fp);

    // size_t gg_size = ig[n];
    // complex<double>* gg = new complex<double>[gg_size];
    // fp = fopen((folder + "/gg").c_str(), "rb");
    // for (size_t i = 0, j = 0; i < ijg[ijg.size() - 1]; i++, j++) {
    //     size_t rdb = fread(&b, sizeof(b), 1, fp);
    //     if (rdb == 0) throw;
    //     if (ijg[j + 1] - ijg[j] == 2) {
    //         rdb = fread(&c, sizeof(c), 1, fp);
    //         if (rdb == 0) throw;
    //         gg[j].imag(c);
    //         i++;
    //     }
    //     gg[j].real(b);
    // }
    // fclose(fp);

    // complex<double>* rp = new complex<double>[n];
    // fp = fopen((folder + "/pr").c_str(), "rb");
    // for (size_t i = 0; i < n; i++) {
    //     size_t rdb = fread(&b, sizeof(b), 1, fp);
    //     if (rdb == 0) throw;
    //     rp[i].real(b);
    //     rdb = fread(&c, sizeof(c), 1, fp);
    //     if (rdb == 0) throw;
    //     rp[i].imag(c);
    // }
    // fclose(fp);

    // ijg.clear();
    // idi.clear();

    // complex<double>* x0 = new complex<double>[n];

    // std::cout << "Read sizes:\n"
    //           << "  ig:  " << 2 * ig_size << '\n'
    //           << "  idi: " << 2 * (n + 1) << '\n'
    //           << "  di:  " << 2 * n << "\n"
    //           << "  jg:  " << 2 * jg_size << '\n'
    //           << "  ijg: " << 2 * (ig[n] + 1) << '\n'
    //           << "  gg:  " << 2 * gg_size << "\n"
    //           << "  rp:  " << 2 * n << '\n'
    //           << "  x0:  " << 2 * n << "\n";
}