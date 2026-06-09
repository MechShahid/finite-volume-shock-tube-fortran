"""Plot numerical and exact solutions for the 1D Sod shock tube problem."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

GAMMA = 1.4
X0 = 0.5
T_FINAL = 0.25
LEFT = (1.0, 0.0, 1.0)      # rho, u, p
RIGHT = (0.125, 0.0, 0.1)   # rho, u, p
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)


def pressure_function(p, state):
    rho_k, u_k, p_k = state
    a_k = np.sqrt(GAMMA * p_k / rho_k)
    if p > p_k:  # shock
        A = 2.0 / ((GAMMA + 1.0) * rho_k)
        B = (GAMMA - 1.0) / (GAMMA + 1.0) * p_k
        f = (p - p_k) * np.sqrt(A / (p + B))
        df = np.sqrt(A / (p + B)) * (1.0 - 0.5 * (p - p_k) / (p + B))
    else:  # rarefaction
        pr = p / p_k
        f = 2.0 * a_k / (GAMMA - 1.0) * (pr ** ((GAMMA - 1.0) / (2.0 * GAMMA)) - 1.0)
        df = (1.0 / (rho_k * a_k)) * pr ** (-(GAMMA + 1.0) / (2.0 * GAMMA))
    return f, df


def star_pressure_velocity():
    rho_l, u_l, p_l = LEFT
    rho_r, u_r, p_r = RIGHT
    a_l = np.sqrt(GAMMA * p_l / rho_l)
    a_r = np.sqrt(GAMMA * p_r / rho_r)
    p_guess = max(1e-8, 0.5 * (p_l + p_r) - 0.125 * (u_r - u_l) * (rho_l + rho_r) * (a_l + a_r))
    p = p_guess
    for _ in range(50):
        f_l, df_l = pressure_function(p, LEFT)
        f_r, df_r = pressure_function(p, RIGHT)
        p_new = p - (f_l + f_r + u_r - u_l) / (df_l + df_r)
        if p_new < 0:
            p_new = 0.5 * p
        if abs(p_new - p) / (p_new + p + 1e-12) < 1e-12:
            p = p_new
            break
        p = p_new
    f_l, _ = pressure_function(p, LEFT)
    f_r, _ = pressure_function(p, RIGHT)
    u_star = 0.5 * (u_l + u_r) + 0.5 * (f_r - f_l)
    return p, u_star


def exact_solution(x):
    p_star, u_star = star_pressure_velocity()
    rho_l, u_l, p_l = LEFT
    rho_r, u_r, p_r = RIGHT
    a_l = np.sqrt(GAMMA * p_l / rho_l)
    a_r = np.sqrt(GAMMA * p_r / rho_r)
    xi = (x - X0) / T_FINAL
    rho = np.zeros_like(x)
    u = np.zeros_like(x)
    p = np.zeros_like(x)

    for j, s in enumerate(xi):
        if s <= u_star:  # left side of contact
            if p_star > p_l:  # left shock
                s_l = u_l - a_l * np.sqrt((GAMMA + 1.0) / (2.0 * GAMMA) * p_star / p_l + (GAMMA - 1.0) / (2.0 * GAMMA))
                if s <= s_l:
                    rho[j], u[j], p[j] = rho_l, u_l, p_l
                else:
                    rho[j] = rho_l * ((p_star / p_l + (GAMMA - 1.0) / (GAMMA + 1.0)) / ((GAMMA - 1.0) / (GAMMA + 1.0) * p_star / p_l + 1.0))
                    u[j], p[j] = u_star, p_star
            else:  # left rarefaction
                a_star_l = a_l * (p_star / p_l) ** ((GAMMA - 1.0) / (2.0 * GAMMA))
                s_head = u_l - a_l
                s_tail = u_star - a_star_l
                if s <= s_head:
                    rho[j], u[j], p[j] = rho_l, u_l, p_l
                elif s > s_tail:
                    rho[j] = rho_l * (p_star / p_l) ** (1.0 / GAMMA)
                    u[j], p[j] = u_star, p_star
                else:
                    u[j] = 2.0 / (GAMMA + 1.0) * (a_l + 0.5 * (GAMMA - 1.0) * u_l + s)
                    a = 2.0 / (GAMMA + 1.0) * (a_l + 0.5 * (GAMMA - 1.0) * (u_l - s))
                    rho[j] = rho_l * (a / a_l) ** (2.0 / (GAMMA - 1.0))
                    p[j] = p_l * (a / a_l) ** (2.0 * GAMMA / (GAMMA - 1.0))
        else:  # right side of contact
            if p_star > p_r:  # right shock
                s_r = u_r + a_r * np.sqrt((GAMMA + 1.0) / (2.0 * GAMMA) * p_star / p_r + (GAMMA - 1.0) / (2.0 * GAMMA))
                if s >= s_r:
                    rho[j], u[j], p[j] = rho_r, u_r, p_r
                else:
                    rho[j] = rho_r * ((p_star / p_r + (GAMMA - 1.0) / (GAMMA + 1.0)) / ((GAMMA - 1.0) / (GAMMA + 1.0) * p_star / p_r + 1.0))
                    u[j], p[j] = u_star, p_star
            else:  # right rarefaction
                a_star_r = a_r * (p_star / p_r) ** ((GAMMA - 1.0) / (2.0 * GAMMA))
                s_head = u_r + a_r
                s_tail = u_star + a_star_r
                if s >= s_head:
                    rho[j], u[j], p[j] = rho_r, u_r, p_r
                elif s <= s_tail:
                    rho[j] = rho_r * (p_star / p_r) ** (1.0 / GAMMA)
                    u[j], p[j] = u_star, p_star
                else:
                    u[j] = 2.0 / (GAMMA + 1.0) * (-a_r + 0.5 * (GAMMA - 1.0) * u_r + s)
                    a = 2.0 / (GAMMA + 1.0) * (a_r - 0.5 * (GAMMA - 1.0) * (u_r - s))
                    rho[j] = rho_r * (a / a_r) ** (2.0 / (GAMMA - 1.0))
                    p[j] = p_r * (a / a_r) ** (2.0 * GAMMA / (GAMMA - 1.0))
    return rho, u, p


def load_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, comments="#", skip_header=4)


def plot_for_n(nx):
    files = {
        "Godunov": RESULTS / f"godunov_n{nx}.csv",
        "MUSCL": RESULTS / f"muscl_n{nx}.csv",
        "Hybrid": RESULTS / f"hybrid_n{nx}.csv",
    }
    x_exact = np.linspace(0, 1, 2000)
    rho_e, u_e, p_e = exact_solution(x_exact)
    variables = [("rho", "Density", rho_e), ("u", "Velocity", u_e), ("p", "Pressure", p_e)]

    for key, ylabel, exact_values in variables:
        plt.figure(figsize=(7, 4.5))
        for label, path in files.items():
            data = load_csv(path)
            plt.plot(data["x"], data[key], label=label)
        plt.plot(x_exact, exact_values, label="Exact", linestyle="--")
        plt.xlabel("x")
        plt.ylabel(ylabel)
        plt.title(f"Sod shock tube: {ylabel}, nx={nx}, t={T_FINAL}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(FIGURES / f"{key}_comparison_n{nx}.png", dpi=200)
        plt.close()


if __name__ == "__main__":
    plot_for_n(201)
    plot_for_n(1001)
    print(f"Figures written to: {FIGURES}")
