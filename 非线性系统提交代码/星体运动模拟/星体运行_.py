import math
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import imageio
from mpl_toolkits.mplot3d import Axes3D


def parameter():

    n = int(input('输入星体个数（建议大于3小于10），随机产生参数并运行，或输入0自行设置参数\n'))

    step = 10000  # 步长，可变

    # 两种输入参数方法
    if n != 0:
        time = int(input('输入总运行时间（建议大于10000000）\n'))
        data = [[] for i in range(n)]
        for i in range(n):
            data[i] = [random.randint(1, 10)*1e30, random.randint(-10, 10)*1e9,
                       random.randint(-10, 10)*1e9, random.randint(-10, 10)*1e9,
                       random.randint(-5, 5)*1e3, random.randint(-5, 5)*1e3,
                       random.randint(-5, 5)*1e3]

    else:
        n = int(input('输入星体个数: '))
        data = [[0]]*n
        for i in range(n):
            print('请输入第%d个星体的质量、横坐标、纵坐标、水平方向速度、竖直方向速度，采用国际单位，空格间隔' % (i+1))
            data[i] = list(map(int, list(input().split())))

        time, step = map(int, input('输入总运行时间和时间步长，国际单位，空格间隔\n').split())

    s = []
    for i in range(n):
        s.append(data[i][0])
    for i in range(n):
        if max(s) != min(s):
            s[i] = (s[i] - min(s)) * 30 / (max(s) - min(s)) + 20
        else:
            s[i] = 40

    return data, n, time, step, s


def velocity_processing(data, n):

    velocity = [[0]*3 for i in range(n)]

    # 一阶导为速度
    for i in range(n):
        velocity[i][0] = data[i][4]
        velocity[i][1] = data[i][5]
        velocity[i][2] = data[i][6]

    return velocity


def acceleration_processing(data, n):

    g_constant = 6.67e-11

    acceleration = [[0]*3 for i in range(n)]
    # 二阶导为加速度
    for i in range(n):
        for j in range(i+1, n):

            r2 = (data[i][1] - data[j][1]) ** 2 + ((data[i][2] - data[j][2]) ** 2) + ((data[i][3] - data[j][3]) ** 2)
            f = data[i][0] * data[j][0] * g_constant / r2

            axi = (f / data[i][0]) * ((data[j][1] - data[i][1]) / math.sqrt(r2))
            axj = (f / data[j][0]) * ((data[i][1] - data[j][1]) / math.sqrt(r2))
            ayi = (f / data[i][0]) * ((data[j][2] - data[i][2]) / math.sqrt(r2))
            ayj = (f / data[j][0]) * ((data[i][2] - data[j][2]) / math.sqrt(r2))
            azi = (f / data[i][0]) * ((data[j][3] - data[i][3]) / math.sqrt(r2))
            azj = (f / data[j][0]) * ((data[i][3] - data[j][3]) / math.sqrt(r2))

            acceleration[i][0] += axi
            acceleration[i][1] += ayi
            acceleration[i][2] += azi
            acceleration[j][0] += axj
            acceleration[j][1] += ayj
            acceleration[j][2] += azj

    return acceleration


def tolist(data, k_, l_, n, step):

    li = [[]]*n
    for i in range(n):
        li[i] = [0] + list(k_[i]) + list(l_[i])

    return list(np.array(data) + step / 2 * np.array(li))


def runge_kutta(data, n, array, t_array, time, step):

    step2 = time_optimizer(data, n, step)  # 根据星体间的间隔调整步长，提高精度。在过近时终止程序。

    # 使用四阶龙格库塔法求数值解
    k_1 = velocity_processing(data, n)
    l_1 = acceleration_processing(data, n)
    m_1 = tolist(data, k_1, l_1, n, step2)
    k_2 = velocity_processing(m_1, n)
    l_2 = acceleration_processing(m_1, n)
    m_2 = tolist(data, k_2, l_2, n, step2)
    k_3 = velocity_processing(m_2, n)
    l_3 = acceleration_processing(m_2, n)
    m_3 = tolist(data, k_3, l_3, n, step2)
    k_4 = velocity_processing(m_3, n)
    l_4 = acceleration_processing(m_3, n)

    for i in range(n):
        data[i][1] += (k_1[i][0] + 2 * k_2[i][0] + 2 * k_3[i][0] + k_4[i][0]) * step2 / 6
        data[i][2] += (k_1[i][1] + 2 * k_2[i][1] + 2 * k_3[i][1] + k_4[i][1]) * step2 / 6
        data[i][3] += (k_1[i][2] + 2 * k_2[i][2] + 2 * k_3[i][2] + k_4[i][2]) * step2 / 6
        data[i][4] += (l_1[i][0] + 2 * l_2[i][0] + 2 * l_3[i][0] + l_4[i][0]) * step2 / 6
        data[i][5] += (l_1[i][1] + 2 * l_2[i][1] + 2 * l_3[i][1] + l_4[i][1]) * step2 / 6
        data[i][6] += (l_1[i][2] + 2 * l_2[i][2] + 2 * l_3[i][2] + l_4[i][2]) * step2 / 6
        array[i][0].append(data[i][1])
        array[i][1].append(data[i][2])
        array[i][2].append(data[i][3])

    t_array.append(time+step2)

    return data, step2, array, t_array


def plotting(array, n, time, step, el, az, f, s):

    sns.set(style='darkgrid')
    plt.ion()
    track_length = 2000  # 控制绘制的轨迹长度

    # time // 1000可调整，控制绘图的间隔，除数越大程序运行速度越快
    if time[-1] // 1000 > f:

        plt.clf()
        fig = plt.gcf()
        ax = fig.gca(projection='3d')
        ax.view_init(el, az)

        # 绘制星体一定运行时间内的运动轨迹
        for i in range(n):
            if len(array[i][0]) < track_length:
                ax.plot(array[i][0], array[i][1], array[i][2])
                ax.scatter(array[i][0][-1], array[i][1][-1], array[i][2][-1], s=s[i])
            else:
                index = 0
                for j in range(len(time)-1, 0):
                    if time[n] < time[-1] - track_length * step:
                        index = n
                    break
                ax.plot(array[i][0][index:-1], array[i][1][index:-1], array[i][2][index:-1])
                ax.scatter(array[i][0][-1], array[i][1][-1], array[i][2][-1], s=s[i])

        # 绘制星体的完整运动轨迹
        # for i in range(n):
            # ax.plot(array[i][0], array[i][1], array[i][2])
            # ax.scatter(array[i][0][-1], array[i][1][-1], array[i][2][-1], s=s[i])

        """ax.legend(["[" + str("%e" % array[i][0][-1]) + " , " + str("%e" % array[i][1][-1]) + " , "
                   + str("%e" % array[i][2][-1]) + "]" for i in range(n)], loc='upper right')"""
        plt.title('program has run for %.2f units of time' % time[-1])
        plt.pause(0.01)
        el, az = ax.elev, ax.azim
        plt.savefig('%d.jpg' % f)  # 保存每一次绘制的图片，可能会降低运行速度
        plt.ioff()
        f += 1

    return el, az, f


# 调整步长
def time_optimizer(data, n, step):

    g_constant = 6.67e-11
    constant = 0.00001  # 该参数可以调整，在尽可能小的情况下可以避免星体距离过近时因为步长过长导致误差较大的情况

    for i in range(n):
        for j in range(i+1, n):
            cdist = math.sqrt((data[i][1] - data[j][1]) ** 2 + (data[i][2] - data[j][2]) ** 2
                              + (data[i][3] - data[j][3]) ** 2)
            force = g_constant * data[i][0] * data[j][0] / cdist ** 2
            ai = force / data[i][0]
            aj = force / data[j][0]
            if math.sqrt(constant * cdist / ai) <= constant * step \
                    or math.sqrt(constant * cdist / aj) <= constant * step:
                print("星体之间距离过近，程序终止")
                exit()

            # 调整步长
            if math.sqrt(constant * cdist / ai) <= step:
                step = math.sqrt(constant * cdist / ai)
            if math.sqrt(constant * cdist / aj) <= step:
                step = math.sqrt(constant * cdist / aj)

    return step


# 生成gif图片
def create_gif(f):

    images = []
    for i in range(f):
        im = imageio.imread('%d.jpg' % i)
        images.append(im)
    imageio.mimsave('结果.gif', images, 'GIF', duration=0.01)


# #######################主程序###############################

# 导入星体的各项参数，程序总运行时间，时间步长
star_parameter, star_number, total_time, time_step, size = parameter()

t = 0
position_array = [[[] for j in range(3)] for i in range(star_number)]  # 记录星体的位置列表
time_array = []  # 记录时间
flag2 = 0
azim = -60  # 绘图视角
elev = 30  # 绘图视角

while t < total_time:

    # 利用龙格库塔法更新坐标
    star_parameter, t_step, position_array, time_array = runge_kutta(star_parameter, star_number,
                                                                     position_array, time_array, t, time_step)

    # 绘图
    elev, azim, flag2 = plotting(position_array, star_number, time_array, time_step, elev, azim, flag2, size)

    t += t_step

# 生成动态图，可能会降低运行速度
create_gif(flag2)

print('运行时间结束，程序终止')
