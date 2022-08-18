import numpy as np
import matplotlib.pyplot as plt

plt.style.use("ggplot")
plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
plt.rcParams['axes.unicode_minus'] = False


class MaterialNonlinear():
    def __init__(self) -> None:
        self.sig = 0.0
        self.eps = 0.0
        self.delta_eps = 0.0
        self.sig_p = 0.0
        self.eps_p = 0.0
        self.delta_eps_p = 0.0
        self.E = 1.0
        self.E_p = 1.0
        self.status = 0
        self.status_p = 0

    def set_strain(self, eps):
        self.eps = eps
        self.delta_eps = eps - self.eps_p

    def confirm(self):
        self.E_p = self.E
        self.sig_p = self.sig
        self.eps_p = self.eps
        self.delta_eps_p = self.delta_eps
        self.status_p = self.status


class Bilinear(MaterialNonlinear):
    def __init__(self, E_0, sig_y, alpha=0.0) -> None:
        super().__init__()
        self.E_0 = E_0  # 弹性状态弹性模量
        self.sig_y = sig_y  # 屈服应力
        self.alpha = alpha  # 屈服状态弹性模量除以弹性状态弹性模量
        self.E_a = alpha * E_0  # 屈服状态弹性模量
        self.eps_y = sig_y / E_0  # 屈服应变
        self.E = E_0
        self.E_p = E_0

    def calc_stress_tangent(self):

        sig_try = self.sig_p + self.delta_eps * self.E_p  # 试应力
        sig_upper = self.sig_y + self.E_a * (self.eps - self.eps_y)  # 应力上限
        sig_lower = -self.sig_y + self.E_a * (self.eps + self.eps_y)  # 应力下限
        self.status = self.status_p

        if self.delta_eps == 0.0:  # 无应变增量
            self.sig = self.sig_p
            self.E = self.E_p
        else:  # 有应变增量
            if self.status == 0:  # 弹性状态
                if self.delta_eps > 0:  # 正应变增量
                    if sig_try > sig_upper:  # 试应力超过应力界限
                        self.sig = sig_upper
                        self.E = self.E_a
                        self.status = 1  # 变为屈服状态
                    else:
                        self.sig = sig_try
                        self.E = self.E_0
                else:  # 负应变增量
                    if sig_try < sig_lower:  # 试应力超过应力界限
                        self.sig = sig_lower
                        self.E = self.E_a
                        self.status = 1  # 变为屈服状态
                    else:
                        self.sig = sig_try
                        self.E = self.E_0
            else:  # 屈服状态
                if self.delta_eps_p * self.delta_eps > 0:  # 加载
                    self.sig = sig_try
                    self.E = self.E_a
                else:  # 卸载
                    self.E = self.E_0
                    self.sig = self.sig_p + self.E * self.delta_eps
                    self.status = 0  # 变为弹性状态
                    # 考虑卸载过大可能导致反向屈服
                    if self.sig > sig_upper:  # 试应力超过应力界限
                        self.sig = sig_upper
                        self.E = self.E_a
                        self.status = 1  # 变为屈服状态
                    elif self.sig < sig_lower:  # 试应力超过应力界限
                        self.sig = sig_lower
                        self.E = self.E_a
                        self.status = 1  # 变为屈服状态


def test_bilinear(strain_hist_type="random"):
    n = 10001
    t = np.linspace(0, 5, n)

    if strain_hist_type == "line*sine":
        eps = 2.5 * t / 5.0 * np.sin(2.0 * np.pi / 0.5 * t)
    elif strain_hist_type == "sine*sine":
        eps = 2.5 * np.sin(2.0 * np.pi / 10.0 * t) * np.sin(2.0 * np.pi / 0.5 * t)
    elif strain_hist_type == "random":
        n_e = 21
        t_e = np.linspace(0, 5, n_e)
        eps_e = 5.0 * (np.random.rand(n_e) - 0.5)
        eps = np.interp(t, t_e, eps_e) * np.sin(2.0 * np.pi / 0.5 * t)

    sig = np.zeros(n)

    bl = Bilinear(1.0, 1.0, 0.0)

    for i in range(n):
        bl.set_strain(eps[i])
        bl.calc_stress_tangent()
        bl.confirm()
        sig[i] = bl.sig

    plt.figure("Bilinear Material", (10, 10))
    ax1 = plt.subplot2grid((4, 1), (0, 0))
    ax1.plot(t, eps)
    ax1.set_title("应变历史")
    ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
    ax2.plot(eps, sig)
    ax2.set_title("应力-应变曲线")
    plt.show()


if __name__ == '__main__':
    test_bilinear("line*sine")
    test_bilinear("sine*sine")
    test_bilinear("random")
    test_bilinear("random")