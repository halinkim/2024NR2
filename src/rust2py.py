import math
from typing import List
import matplotlib.pyplot as plt

PI = 3.14159
G = 6.67430e-8  # cm3 g-1 s-2
C = 2.99792e5  # km s-1
AU = 1.49598e13  # cm
M_SUN = 1.98854e33  # g


def to_phys_vel(v_sim: float) -> float:
    return v_sim / 100000.0 * math.sqrt(G * M_SUN / AU)


def to_sim_vel(v_phys: float) -> float:
    return v_phys * 100000.0 / math.sqrt(G * M_SUN / AU)


def to_phys_time(t_sim: float) -> float:
    return t_sim / math.sqrt(G * M_SUN / AU / AU / AU)


def to_sim_time(t_phys: float) -> float:
    return t_phys * math.sqrt(G * M_SUN / AU / AU / AU)


class Vec2:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y

    def __add__(self, other):
        return Vec2(self.x + other.x, self.y + other.y)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        return self

    def __sub__(self, other):
        return Vec2(self.x - other.x, self.y - other.y)

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        return self

    def __mul__(self, scalar: float):
        return Vec2(self.x * scalar, self.y * scalar)

    def __rmul__(self, scalar: float):
        return Vec2(self.x * scalar, self.y * scalar)

    def __imul__(self, scalar: float):
        self.x *= scalar
        self.y *= scalar
        return self

    def __truediv__(self, scalar: float):
        return Vec2(self.x / scalar, self.y / scalar)

    def __itruediv__(self, scalar: float):
        self.x /= scalar
        self.y /= scalar
        return self

    def __neg__(self):
        return Vec2(-self.x, -self.y)

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def size(self):
        return math.hypot(self.x, self.y)

    def __repr__(self):
        return f"Vec2({self.x}, {self.y})"


class Body:
    def __init__(self, mass: float, pos: Vec2, vel: Vec2):
        self.mass = mass
        self.pos = pos
        self.vel = vel


def init(m1: float, m2: float, a: float, e: float):
    total_m = m1 + m2
    mu = m1 * m2 / total_m
    k = m1 * m2
    r = a * (1.0 + e)
    l = math.sqrt(mu * k * a * (1.0 - e * e))
    theta_dot = l / mu / r / r

    x1 = -r * m2 / total_m
    x2 = r + x1
    v1 = theta_dot * x1
    v2 = theta_dot * x2

    period = math.sqrt(4.0 * PI * PI * mu / k * a * a * a)

    return x1, x2, v1, v2, period


def acc(system: List[Body], ind: int) -> Vec2:
    res = Vec2(0.0, 0.0)
    for j in range(len(system)):
        if ind != j:
            pos_rel = system[ind].pos - system[j].pos
            r = pos_rel.size()
            res += -system[j].mass * pos_rel / (r * r * r)
    return res


def acc_pn(system: List[Body], ind: int) -> Vec2:
    res = Vec2(0.0, 0.0)
    for j in range(len(system)):
        if ind != j:
            pos_rel = system[ind].pos - system[j].pos
            vel_rel = system[ind].vel - system[j].vel
            r = pos_rel.size()
            v = vel_rel.size()
            total_m = system[0].mass + system[1].mass
            eta = system[0].mass * system[1].mass / total_m ** 2
            r_dot = pos_rel.dot(vel_rel) / r
            pn_a = 8.0 / 5.0 * eta * total_m / r * r_dot * (17.0 / 3.0 * total_m / r + 3.0 * v * v)
            pn_b = -8.0 / 5.0 * eta * total_m / r * (3.0 * total_m / r + v * v)

            res += system[j].mass / (r * r) * (pn_a * pos_rel / r + pn_b * vel_rel)
    res /= to_sim_vel(C) ** 5
    return res


def acc_dot(system: List[Body], ind: int) -> Vec2:
    res = Vec2(0.0, 0.0)
    for j in range(len(system)):
        if ind != j:
            pos_rel = system[ind].pos - system[j].pos
            r = pos_rel.size()
            vel_rel = system[ind].vel - system[j].vel
            r_dot = pos_rel.dot(vel_rel) / r
            res += -system[j].mass * (vel_rel / (r ** 3) - 3.0 * r_dot * pos_rel / (r ** 4))
    return res


def acc_dot_pn(system: List[Body], ind: int) -> Vec2:
    res = Vec2(0.0, 0.0)
    for j in range(len(system)):
        if ind != j:
            pos_rel = system[ind].pos - system[j].pos
            vel_rel = system[ind].vel - system[j].vel
            acc_rel = acc(system, ind)
            r = pos_rel.size()
            v = vel_rel.size()
            total_m = system[0].mass + system[1].mass
            eta = system[0].mass * system[1].mass / total_m ** 2
            r_dot = pos_rel.dot(vel_rel) / r
            v_dot = vel_rel.dot(acc_rel) / v
            pn_a = 8.0 / 5.0 * eta * total_m / r * r_dot * (17.0 / 3.0 * total_m / r + 3.0 * v * v)
            pn_b = -8.0 / 5.0 * eta * total_m / r * (3.0 * total_m / r + v * v)
            r_ddot = (v * v + pos_rel.dot(acc_rel) - r_dot * r_dot) / r
            pn_a_dot = (
                    8.0 * total_m * eta / 5.0 / r * (17.0 * total_m / 3.0 / r + 3.0 * v * v) * r_ddot
                    - 8.0 * total_m * eta / 5.0 / r / r * (17.0 * total_m / 3.0 / r + 3.0 * v * v) * r_dot * r_dot
                    + 8.0 * total_m * eta / 5.0 / r * (-17.0 * total_m * r_dot / 3.0 / r / r + 6.0 * v * v_dot) * r_dot
            )
            pn_b_dot = (
                    8.0 / 5.0 * eta * total_m / r / r * r_dot * (3.0 * total_m / r + v * v)
                    - 8.0 / 5.0 * eta * total_m / r * (-3.0 * total_m / r / r * r_dot + 2.0 * vel_rel.dot(acc_rel))
            )
            term1 = -2.0 * system[j].mass / (r ** 3) * (pn_a * pos_rel / r + pn_b * vel_rel) * r_dot
            term2 = system[j].mass / (r * r) * (
                    -pn_a * pos_rel / (
                        r * r) * r_dot + pn_a * vel_rel / r + pn_b * acc_rel + pos_rel * pn_a_dot / r + vel_rel * pn_b_dot
            )
            res += term1 + term2
    res /= to_sim_vel(C) ** 5
    return res


def integrator(system: List[Body], dt: float, use_pn: bool) -> List[Body]:
    a_0 = []
    a_0_dot = []
    for i in range(len(system)):
        a_0.append(acc(system, i))
        a_0_dot.append(acc_dot(system, i))
    if use_pn:
        for i in range(len(system)):
            a_0[i] += acc_pn(system, i)
            a_0_dot[i] += acc_dot_pn(system, i)
    system_p = [Body(body.mass, Vec2(body.pos.x, body.pos.y), Vec2(body.vel.x, body.vel.y)) for body in system]
    for i in range(len(system)):
        system_p[i].pos += system[i].vel * dt + a_0[i] * dt * dt / 2.0 + a_0_dot[i] * dt * dt * dt / 6.0
        system_p[i].vel += a_0[i] * dt + a_0_dot[i] * dt * dt / 2.0
    a_p = []
    a_p_dot = []
    for i in range(len(system)):
        a_p.append(acc(system_p, i))
        a_p_dot.append(acc_dot(system_p, i))
    if use_pn:
        for i in range(len(system)):
            a_p[i] += acc_pn(system_p, i)
            a_p_dot[i] += acc_dot_pn(system_p, i)
    a_0_ddot = []
    a_0_dddot = []
    for i in range(len(system)):
        a_0_ddot.append(-6.0 * (a_0[i] - a_p[i]) / dt / dt - 2.0 * (2.0 * a_0_dot[i] + a_p_dot[i]) / dt)
        a_0_dddot.append(12.0 * (a_0[i] - a_p[i]) / dt ** 3 + 6.0 * (a_0_dot[i] + a_p_dot[i]) / dt ** 2)
    for i in range(len(system)):
        system_p[i].pos += a_0_ddot[i] * dt ** 4 / 24.0 + a_0_dddot[i] * dt ** 5 / 120.0
        system_p[i].vel += a_0_ddot[i] * dt ** 3 / 6.0 + a_0_dddot[i] * dt ** 4 / 24.0
    return system_p


def solver(m1: float, m2: float, a: float, e: float, n_period: float, n_step: int, use_pn: bool):
    x1, x2, v1, v2, period = init(m1, m2, a, e)
    tot_t = period * n_period
    acc_t = 0.0

    his_x1 = []
    his_y1 = []
    his_x2 = []
    his_y2 = []
    system = [Body(m1, Vec2(x1, 0.0), Vec2(0.0, v1)),
              Body(m2, Vec2(x2, 0.0), Vec2(0.0, v2))]
    sch_r = 2.0 * (m1 + m2) / to_sim_vel(C) ** 2
    print("Initial value:")
    print(
        f"\tBody1\n\t\tMass = {m1}\tX = {system[0].pos.x:.12e}\tY = {system[0].pos.y:.12e}\n\t\tVX = {system[0].vel.x:.12e}\tVY = {system[0].vel.y:.12e}")
    print(
        f"\tBody2\n\t\tMass = {m2}\tX = {system[1].pos.x:.12e}\tY = {system[1].pos.y:.12e}\n\t\tVX = {system[1].vel.x:.12e}\tVY = {system[1].vel.y:.12e}")
    print(f"\tSeparation = {system[1].pos.x - system[0].pos.x:.12e}")
    print(f"\tTarget Time = {tot_t:.12e}")
    print(f"\tR_sch of total mass = {sch_r:.12e}\n")

    dt_alpha = 30.0 * (1.0 + e) * (1.0 + int(use_pn))
    if not use_pn:
        pos = (system[0].pos - system[1].pos).size()
        vel = (system[0].vel - system[1].vel).size()
        est_dt1 = pos / vel / dt_alpha
        est_dt2 = (2.0 * a - pos) ** 2 / vel / pos / dt_alpha
        est_n1 = tot_t / (est_dt1 + est_dt2) * 2.0
        est_n2 = tot_t / math.sqrt(est_dt1 * est_dt2)
        est_n = (est_n1 + est_n2) * 0.5
        print(f"\tEst. N = {est_n}")
        if abs(est_n - n_step) / n_step * 100.0 > 30.0:
            dt_alpha *= n_step / est_n / 1.2

    print("Starting solve..")
    early_break = False

    for it in range(1, n_step + 1):
        pos = (system[0].pos - system[1].pos).size()
        if pos < sch_r:
            print(f"\tReach Sch limit at iter: {it}")
            early_break = True
            break

        vel = (system[0].vel - system[1].vel).size()

        dt = pos / vel / dt_alpha
        if acc_t + dt > tot_t:
            dt = tot_t - acc_t
        new_system = integrator(system, dt, use_pn)
        system = new_system

        if it % 1 == 0:
            his_x1.append(system[0].pos.x)
            his_y1.append(system[0].pos.y)
            his_x2.append(system[1].pos.x)
            his_y2.append(system[1].pos.y)

        acc_t += dt
        if acc_t >= tot_t:
            print(f"\tTotal t achieved at iter: {it}")
            early_break = True
            break

    if not early_break:
        print(f"\tComplete full iter: {n_step}")
    pos = (system[0].pos - system[1].pos).size()
    print("\nFinal value:")
    print(
        f"\tBody1\n\t\tMass = {m1}\tX = {system[0].pos.x:.12e}\tY = {system[0].pos.y:.12e}\n\t\tVX = {system[0].vel.x:.12e}\tVY = {system[0].vel.y:.12e}")
    print(
        f"\tBody2\n\t\tMass = {m2}\tX = {system[1].pos.x:.12e}\tY = {system[1].pos.y:.12e}\n\t\tVX = {system[1].vel.x:.12e}\tVY = {system[1].vel.y:.12e}")
    print(f"\tSeparation = {pos:.12e}")
    print(f"\tTotal Time = {acc_t:.12e} = {to_phys_time(acc_t):.12e} s")
    print(f"\tBody1 V = {to_phys_vel(system[0].vel.size()):.12e} km/s")
    print(f"\tBody2 V = {to_phys_vel(system[1].vel.size()):.12e} km/s")

    return his_x1, his_y1, his_x2, his_y2


def main():
    print("Unit Conversion")
    print(f"\tsim to phys.. 1 = {to_phys_vel(1.0):.12e} km/s,\t1 = {to_phys_time(1.0):.12e} s")
    print(f"\tphys to sim.. 1 km/s = {to_sim_vel(1.0):.12e},\t1 s = {to_sim_time(1.0):.12e}\n")

    m1 = 1.0
    m2 = 1.1163969545495692e-16
    a = 17.8
    e = 0.967
    n_period = 50.0
    n_step = 40000
    use_pn = False

    his_x1, his_y1, his_x2, his_y2 = solver(m1, m2, a, e, n_period, n_step, use_pn)

    print("Done!")

    plt.plot(his_x1, his_y1)
    plt.plot(his_x2, his_y2)
    plt.grid()
    plt.axis('equal')
    plt.show()


if __name__ == "__main__":
    main()
