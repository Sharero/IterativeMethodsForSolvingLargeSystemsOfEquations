import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

out_dir = Path("output/")

def read_data(path):
    iterations = []
    residuals = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or ';' not in line:
                continue
            iter_str, res_str = line.split(';')
            iterations.append(int(iter_str.strip()))
            residuals.append(float(res_str.strip()))
    return np.array(iterations), np.array(residuals)

def draw_plot_residual_vs_iterations(slae_number, method_name):
    file1 = out_dir / "omp1" / f"{method_name}_Di_{slae_number}.txt"
    file2 = out_dir / "omp1" / f"{method_name}_Di_Smooth_{slae_number}.txt"

    it1, res1 = read_data(file1)
    it2, res2 = read_data(file2)

    plt.figure(figsize=(8,5))
    plt.plot(it1, res1, linestyle='-', label=f'{method_name}')
    plt.plot(it2, res2, linestyle='-', label=f'{method_name}_Smooth')
    plt.yscale('log')
    plt.xlabel('Количество итераций')
    plt.ylabel('Норма относительной невязки')
    plt.title(f'Зависимость невязки от числа итераций для СЛАУ {slae_number}')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # for i in range(1,5):
    #     draw_plot_residual_vs_iterations(i, "COCG")

    for i in range(1,5):
        draw_plot_residual_vs_iterations(i, "COCR")

    