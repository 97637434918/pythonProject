import numpy as np
import matplotlib.pyplot as plt

from nl import Bilinear

plt.style.use("ggplot")
plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
plt.rcParams['axes.unicode_minus'] = False


def draw_response(title, ta, a, t, u):
    plt.figure(title, (12, 6))
    plt.subplot(2, 1, 1)
    plt.plot(ta, a, label=r"输入地震波加速度时程")
    plt.grid(True)
    plt.legend()
    plt.xlim(0, t[-1])
    plt.subplot(2, 1, 2)
    plt.plot(t, u, label=r"SDOF体系位移响应时程")
    plt.xlabel(r"时间 (s)")
    plt.grid(True)
    plt.legend()
    plt.xlim(0, t[-1])
    plt.show()


def draw_hysteretic_curve(u, f):
    plt.figure("Hysteretic Curve", (12, 9))
    plt.plot(u, f)
    plt.grid(True)
    plt.xlabel(r"变形")
    plt.ylabel(r"力")
    plt.show()


def newton_raphson(bl: Bilinear, kh, ph, tol=1e-8, maxiter=50):
    x = ph / (kh + bl.E)
    bl.set_strain(x)
    bl.calc_stress_tangent()
    error = kh * x + bl.sig - ph
    iter = 0
    while abs(error) > tol and iter < maxiter:
        x = x - error / (kh + bl.E)
        bl.set_strain(x)
        bl.calc_stress_tangent()
        error = kh * x + bl.sig - ph
        iter += 1
    if iter == maxiter: print(error)
    bl.confirm()
    return x


def solve_sdof_eqwave_nmk_bl(omg, zeta, ag, dt, bl: Bilinear):
    n = len(ag)
    omg2 = omg * omg

    gma = 0.5
    bta = 0.25

    u0 = 0.0
    v0 = 0.0

    c1 = 1.0 / bta / dt / dt
    c2 = 1.0 / bta / dt
    c3 = gma / bta / dt
    c4 = 1.0 - gma / bta
    c5 = 1.0 - 0.5 * gma / bta
    c6 = 0.5 / bta - 1.0

    c = 2.0 * zeta * omg
    c7 = c1 + c3 * c
    c8 = c2 - c4 * c
    c9 = c6 - dt * c5 * c

    u = np.zeros(n)
    v = np.zeros(n)
    a = np.zeros(n)
    f = np.zeros(n)

    s = np.zeros(n)

    bl.set_strain(u0)
    bl.calc_stress_tangent()
    bl.confirm()

    a0 = -2.0 * zeta * omg * v0 - bl.sig
    kh = c7
    u[0] = u0
    v[0] = v0
    a[0] = -ag[0] - c * v[0] - bl.sig
    f[0] = bl.sig

    if bl.status_p == 1: s[0] = 1.0

    for i in range(n - 1):
        ph = -ag[i + 1] + c7 * u[i] + c8 * v[i] + c9 * a[i]
        u[i + 1] = newton_raphson(bl, kh, ph)
        v[i + 1] = c3 * (u[i + 1] - u[i]) + c4 * v[i] + dt * c5 * a[i]
        a[i + 1] = -ag[i + 1] - c * v[i + 1] - bl.sig
        f[i + 1] = bl.sig
        if bl.status_p == 1: s[i + 1] = 1.0

    return u, v, f, s


if __name__ == '__main__':
    acc0 = np.loadtxt("EQ-S-3.txt")  # 读取地震波
    dt = 0.005  # 时间间隔
    n = len(acc0)
    t0 = np.linspace(0.0, dt * (n - 1), n)

    # 对时程数据进行补零以显示地震结束后一段时间内的自由振动衰减情况
    ne = round(n * 1.2)
    t = np.linspace(0.0, dt * (ne - 1), ne)
    ag = np.zeros(ne)
    ag[:n] = acc0

    omg = 2.0 * np.pi  # 自振圆频率
    zeta = 0.05  # 固有阻尼比

    k_0 = omg * omg
    f_y = k_0 * 0.012
    bl = Bilinear(omg * omg, f_y, 0.02)

    u, v, f, s = solve_sdof_eqwave_nmk_bl(omg, zeta, ag, dt, bl)

    draw_response("Seismic Response -- newmark-β bilinear spring", t0, acc0, t, u)
    draw_hysteretic_curve(u, f)

    plt.figure("Status", (12, 5))
    plt.plot(t, s)
    plt.xlabel("时间 (s)")
    plt.ylabel("状态")
    plt.yticks((0, 1), ("弹性", "屈服"))
    plt.show()