const PI: f64 = 3.14159;
const G: f64 = 6.67430e-8;  // cm3 g-1 s-2
const C: f64 = 2.99792e5; // km s-1
const AU: f64 = 1.49598e13;  // cm
const M_SUN: f64 = 1.98854e33;  // g

use plotpy::{Curve, Plot, StrError};

fn main() -> Result<(), StrError>{
    println!("Unit Conversion");
    println!("\tsim to phys.. 1 = {:.12e} km/s,\t1 = {:.12e} s", to_phys_vel(1.0), to_phys_time(1.0));
    println!("\tphys to sim.. 1 km/s = {:.12e},\t 1 s = {:.12e}\n", to_sim_vel(1.0), to_sim_time(1.0));

    let args: Vec<String> = std::env::args().collect();

    if args.len() != 2 {
        println!("Use args -> cargo run -- <parameter file path>");
        std::process::exit(1);
    }

    let input_path = &args[1];
    println!("Read input file from \"{}\"", input_path);

    let mut file = std::fs::File::open(input_path).unwrap();
    let mut par = String::new();
    file.read_to_string(&mut par).unwrap();


    let pars: Vec<&str> = par.lines().collect();
    if pars.len() != 7 {
        println!("Number of Parameters should 7.\n\t(m1: float, m2: float, a: float, e: float, n_period: float, n_step: int, use_pn: bool)");
        std::process::exit(1);
    }

    println!("\t{:?}\n", pars);
    let m1: f64 = pars[0].parse().unwrap();
    let m2: f64 = pars[1].parse().unwrap();
    let a: f64 = pars[2].parse().unwrap();
    let e: f64 = pars[3].parse().unwrap();
    let n_period: f64 = pars[4].parse().unwrap();
    let n_step: i32 = pars[5].parse().unwrap();
    let use_pn: bool = parse_bool(pars[6]);

    let (his_x1, his_y1, his_x2, his_y2) = solver(m1, m2, a, e, n_period, n_step, use_pn);
    


    // Comment below to do not plot //

    println!("Save Figure..");
    let mut plot = Plot::new();
    let mut orbit1 = Curve::new();
    let mut orbit2 = Curve::new();
    orbit1.draw(&his_x1, &his_y1);
    orbit2.draw(&his_x2, &his_y2);
    plot.add(&orbit1)
        .set_equal_axes(true);
    plot.add(&orbit2)
        .set_equal_axes(true);
    let fig_path = format!("{input_path}.png");
    plot.save_and_show(fig_path.as_str())?;
    println!("Done!");
    Ok(())
}

fn to_phys_vel(v_sim: f64) -> f64 {
    v_sim / 100000.0 * (G * M_SUN / AU).sqrt()
}
fn to_sim_vel(v_phys: f64) -> f64 {
    v_phys * 100000.0 / (G * M_SUN / AU).sqrt()
}
fn to_phys_time(t_sim: f64) -> f64 {
    t_sim / (G * M_SUN / AU / AU / AU).sqrt()
}
fn to_sim_time(t_phys: f64) -> f64 {
    t_phys * (G * M_SUN / AU / AU / AU).sqrt()
}


fn init(m1: f64, m2: f64, a: f64, e: f64) -> (f64, f64, f64, f64, f64) {
    let total_m = m1 + m2;
    let mu = m1 * m2 / total_m;
    let k = m1 * m2;
    let r = a * (1.0 + e);
    let l = (mu * k * a * (1.0 - e * e)).sqrt();
    let theta_dot = l / mu / r / r;

    let x1 = -r * m2 / total_m;
    let x2 = r + x1;
    // let x2 = r * m1 / total_m;
    let v1 = theta_dot * x1;
    let v2 = theta_dot * x2;

    let period = (4.0 * PI * PI * mu / k * a * a * a).sqrt();

    (x1, x2, v1, v2, period)
}

fn acc(system: &Vec<Body>, ind: usize) -> Vec2{
    let mut res = Vec2::ZERO;
    for j in 0..system.len() {
        if ind != j {
            let pos_rel = system[ind].pos - system[j].pos;
            let r = pos_rel.size();
            res -= system[j].mass * pos_rel / r / r / r;
        }
    }
    // println!("acc: {:?}", res);
    res
}

fn acc_pn(system: &Vec<Body>, ind: usize) -> Vec2{
    let mut res = Vec2::ZERO;
    for j in 0..system.len() {
        if ind != j {
            let pos_rel = system[ind].pos - system[j].pos;
            let vel_rel = system[ind].vel - system[j].vel;
            let r = pos_rel.size();
            let v = vel_rel.size();
            let total_m = system[0].mass + system[1].mass;
            let eta = system[0].mass * system[1].mass / total_m / total_m;
            let r_dot = pos_rel.dot(vel_rel) / r;
            let pn_a = 8.0 / 5.0 * eta * total_m / r * r_dot * (17.0 / 3.0 * total_m / r + 3.0 * v * v);
            let pn_b = - 8.0 / 5.0 * eta * total_m / r * (3.0 * total_m / r + v * v);

            res += system[j].mass / r / r * (pn_a * pos_rel / r + pn_b * vel_rel);
        }
    }
    // println!("acc_pn: {:?}", res);
    res / to_sim_vel(C).powi(5)
}

fn acc_dot(system: &Vec<Body>, ind: usize) -> Vec2 {
    let mut res = Vec2::ZERO;
    for j in 0..system.len() {
        if ind != j {
            let pos_rel = system[ind].pos - system[j].pos;
            let r = pos_rel.size();
            let vel_rel = system[ind].vel - system[j].vel;
            let r_dot = pos_rel.dot(vel_rel) / r;
            res -= system[j].mass * (vel_rel / r / r / r - 3.0 * r_dot * pos_rel / r / r / r / r);
        }
    }
    // println!("acc_dot: {:?}", res);
    res
}

fn acc_dot_pn(system: &Vec<Body>, ind: usize) -> Vec2{
    let mut res = Vec2::ZERO;
    for j in 0..system.len() {
        if ind != j {
            let pos_rel = system[ind].pos - system[j].pos;
            let vel_rel = system[ind].vel - system[j].vel;
            let acc_rel = acc(&system, ind);
            let r = pos_rel.size();
            let v = vel_rel.size();
            let total_m = system[0].mass + system[1].mass;
            let eta = system[0].mass * system[1].mass / total_m / total_m;
            let r_dot = pos_rel.dot(vel_rel) / r;
            let v_dot = vel_rel.dot(acc_rel) / v;
            let pn_a = 8.0 / 5.0 * eta * total_m / r * r_dot * (17.0 / 3.0 * total_m / r + 3.0 * v * v);
            let pn_b = - 8.0 / 5.0 * eta * total_m / r * (3.0 * total_m / r + v * v);
            let r_ddot = (v * v + pos_rel.dot(acc_rel) - r_dot * r_dot) / r;
            let pn_a_dot = 8.0 * total_m * eta / 5.0 / r * (17.0 * total_m / 3.0 / r + 3.0 * v * v) * r_ddot
                - 8.0 * total_m * eta / 5.0 / r / r * (17.0 * total_m / 3.0 / r + 3.0 * v * v) * r_dot * r_dot
                + 8.0 * total_m * eta / 5.0 / r * (-17.0 * total_m * r_dot / 3.0 / r / r + 6.0 * v * v_dot) * r_dot;
            let pn_b_dot = 8.0 / 5.0 * eta * total_m / r / r * r_dot * (3.0 * total_m / r + v * v)
                            -8.0 / 5.0 * eta * total_m / r * (-3.0 * total_m / r / r * r_dot + 2.0 * vel_rel.dot(acc_rel));
            res += -2.0 * system[j].mass / r / r / r * (pn_a * pos_rel / r + pn_b * vel_rel) * r_dot
                    + system[j].mass / r / r * (-pn_a * pos_rel / r / r * r_dot + pn_a * vel_rel / r + pn_b * acc_rel + pos_rel * pn_a_dot / r + vel_rel * pn_b_dot);
        }
    }
    // println!("acc_dot_pn: {:?}", res);
    res / to_sim_vel(C).powi(5)
}

fn integrator(system: &Vec<Body>, dt: f64, use_pn: bool) -> Vec<Body> {
    let mut a_0: Vec<Vec2> = vec![];
    let mut a_0_dot: Vec<Vec2> = vec![];
    for i in 0..system.len() {
        a_0.push(acc(&system, i));
        a_0_dot.push(acc_dot(&system, i));
    }
    if use_pn {
        for i in 0..system.len() {
            a_0[i] += acc_pn(&system, i);
            a_0_dot[i] += acc_dot_pn(&system, i);
        }
    }
    let mut system_p = system.to_vec();
    for i in 0..system.len() {
        system_p[i].pos += system[i].vel * dt + a_0[i] * dt * dt / 2.0 + a_0_dot[i] * dt * dt * dt / 6.0;
        system_p[i].vel += a_0[i] * dt + a_0_dot[i] * dt * dt / 2.0;
    }
    let mut a_p: Vec<Vec2> = vec![];
    let mut a_p_dot: Vec<Vec2> = vec![];
    for i in 0..system.len() {
        a_p.push(acc(&system_p, i));
        a_p_dot.push(acc_dot(&system_p, i));
    }
    if use_pn {
        for i in 0..system.len() {
            a_p[i] += acc_pn(&system_p, i);
            a_p_dot[i] += acc_dot_pn(&system_p, i);
        }
    }
    let mut a_0_ddot: Vec<Vec2> = vec![];
    let mut a_0_dddot: Vec<Vec2> = vec![];
    for i in 0..system.len() {
        a_0_ddot.push(-6.0 * (a_0[i] - a_p[i]) / dt / dt - 2.0 * (2.0 * a_0_dot[i] + a_p_dot[i]) / dt);
        a_0_dddot.push(12.0 * (a_0[i] - a_p[i]) / dt / dt / dt + 6.0 * (a_0_dot[i] + a_p_dot[i]) / dt / dt);
    }
    for i in 0..system.len() {
        system_p[i].pos += a_0_ddot[i] * dt * dt * dt * dt / 24.0 + a_0_dddot[i] * dt * dt * dt * dt * dt / 120.0;
        system_p[i].vel += a_0_ddot[i] * dt * dt * dt / 6.0 + a_0_dddot[i] * dt * dt * dt * dt / 24.0;
    }
    system_p
}



fn solver(m1: f64, m2: f64, a: f64, e: f64, n_period: f64, n_step: i32, use_pn: bool) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (x1, x2, v1, v2, period) = init(m1, m2, a, e);
    let mut dt: f64;
    let tot_t = period * n_period;
    let mut acc_t = 0f64;

    let mut his_x1: Vec<f64> = Vec::new();
    let mut his_y1: Vec<f64> = Vec::new();
    let mut his_x2: Vec<f64> = Vec::new();
    let mut his_y2: Vec<f64> = Vec::new();
    let mut system = vec![Body::new(m1, Vec2::new(x1, 0.0), Vec2::new(0.0, v1)),
                                 Body::new(m2, Vec2::new(x2, 0.0), Vec2::new(0.0, v2))];
    let sch_r = 2.0 * (m1 + m2) / to_sim_vel(C).powi(2);
    println!("Initial value:");
    println!("\tBody1\n\t\tMass = {}\tX = {:.12e}\tY = {:.12e}\n\t\tVX = {:.12e}\tVY = {:.12e}", m1, system[0].pos.x, system[0].pos.y, system[0].vel.x, system[0].vel.y);
    println!("\tBody2\n\t\tMass = {}\tX = {:.12e}\tY = {:.12e}\n\t\tVX = {:.12e}\tVY = {:.12e}", m2, system[1].pos.x, system[1].pos.y, system[1].vel.x, system[1].vel.y);
    println!("\tSeperation = {:.12e}", system[1].pos.x - system[0].pos.x);
    println!("\tTarget Time = {:.12e}", tot_t);
    println!("\tR_sch of total mass = {:.12e}\n", sch_r);
    
    let mut dt_alpha = 30.0 * (1.0 + e) * (1.0 + use_pn as i32 as f64);
    if !use_pn {
        let pos = (system[0].pos - system[1].pos).size();
        let vel = (system[0].vel - system[1].vel).size();
        let est_dt1 = pos / vel / dt_alpha;
        let est_dt2 = (2.0 * a - pos).powi(2) / vel / pos / dt_alpha;
        let est_n1 = tot_t / (est_dt1 + est_dt2) * 2.0;
        let est_n2 = tot_t / (est_dt1 * est_dt2).powf(0.5);
        let est_n = (est_n1 + est_n2) * 0.5;
        println!("\tEst. N = {}", est_n);
        if ((est_n - n_step as f64).abs() / n_step as f64) * 100.0 > 30.0 {
            dt_alpha *= n_step as f64 / est_n / 1.2;
        }
    }

    // println!("sch_r : {}, dt: {}, tot_t: {}, a/vdt {}", sch_r, dt, tot_t, a / v1 / dt);
    // for i in 0..system.len() {
        // println!("{:?}", system[i]);
    // }

    println!("Starting solve..");
    let mut early_break = false;
    
    for it in 1..n_step+1 {
        let pos = (system[0].pos - system[1].pos).size();
        if pos < sch_r {
            println!("\tReach Sch limit at iter: {}", it);
            early_break = true;
            break
        }

        let vel = (system[0].vel - system[1].vel).size();
        
        dt = pos / vel / dt_alpha;
        if acc_t + dt > tot_t {
            dt = tot_t - acc_t;
        }
        let mut new_system = integrator(&system, dt, use_pn);
        std::mem::swap(&mut system, &mut new_system);
        
        
        if it % 1 == 0 {
            his_x1.push(system[0].pos.x);
            his_y1.push(system[0].pos.y);
            his_x2.push(system[1].pos.x);
            his_y2.push(system[1].pos.y);
        }

        acc_t += dt;
        if acc_t >= tot_t {
            println!("\tTotal t achieved at iter: {}", it);
            early_break = true;
            break;
        }
    }
    if !early_break {
        println!("\tComplete full iter: {}", n_step);
    }
    let pos = (system[0].pos - system[1].pos).size();
    println!("\nFinal value:");
    println!("\tBody1\n\t\tMass = {}\tX = {:.12e}\tY = {:.12e}\n\t\tVX = {:.12e}\tVY = {:.12e}", m1, system[0].pos.x, system[0].pos.y, system[0].vel.x, system[0].vel.y);
    println!("\tBody2\n\t\tMass = {}\tX = {:.12e}\tY = {:.12e}\n\t\tVX = {:.12e}\tVY = {:.12e}", m2, system[1].pos.x, system[1].pos.y, system[1].vel.x, system[1].vel.y);
    println!("\tSeperation = {:.12e}", pos);
    println!("\tTotal Time = {:.12e} = {:.12e} s", acc_t, to_phys_time(acc_t));
    println!("\tBody1 V = {:.12e} km/s", to_phys_vel(system[0].vel.size()));
    println!("\tBody2 V = {:.12e} km/s", to_phys_vel(system[1].vel.size()));
    // for i in 0..system.len() {
    //     println!("{:?}", system[i]);
    // }

    // let dist = ((system[0].pos.x - system[1].pos.x).powi(2) + (system[0].pos.y - system[1].pos.y).powi(2)).sqrt();
    // println!("dist: {}", dist);
    // let vel = (system[0].vel - system[1].vel).size();
    // println!("dt: {}, acc_t: {}, a/vdt {}", dt, acc_t, dist / vel / dt);

    (his_x1, his_y1, his_x2, his_y2)
}

#[derive(Debug, Clone, Copy)]
struct Body {
    mass: f64,
    pos: Vec2,
    vel: Vec2
}

impl Body {
    pub const fn new(mass: f64, pos: Vec2, vel: Vec2) -> Body {
        Body { mass, pos, vel }
    }
}

#[derive(Debug, Clone, Copy)]
struct Vec2 {
    x: f64,
    y: f64,
}

impl Vec2 {
    #[inline]
    pub const fn new(x: f64, y: f64) -> Vec2 {
        Vec2 { x, y }
    }

    pub const ZERO: Vec2 = Vec2::new(0.0, 0.0);

    pub fn dot(&self, other: Vec2) -> f64 {
        self.x * other.x + self.y * other.y
    }

    pub fn size(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}
use std::{ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign}, io::Read, fmt::format};
impl Add for Vec2 {
    type Output = Vec2;
    #[inline]
    fn add(self, other: Vec2) -> Vec2 {
        Vec2 {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl AddAssign for Vec2 {
    #[inline]
    fn add_assign(&mut self, other: Vec2) {
        *self = Vec2 {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Vec2 {
    type Output = Vec2;
    #[inline]
    fn sub(self, other: Vec2) -> Vec2 {
        Vec2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl SubAssign for Vec2 {
    #[inline]
    fn sub_assign(&mut self, other: Vec2) {
        *self = Vec2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Mul<f64> for Vec2 {
    type Output = Vec2;
    #[inline]
    fn mul(self, other: f64) -> Vec2 {
        Vec2 {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl MulAssign<f64> for Vec2 {
    #[inline]
    fn mul_assign(&mut self, other: f64) {
        *self = Vec2 {
            x: self.x * other,
            y: self.y * other,
        };
    }
}

impl Mul<Vec2> for f64 {
    type Output = Vec2;
    #[inline]
    fn mul(self, other: Vec2) -> Vec2 {
        other * self
    }
}

impl Div<f64> for Vec2 {
    type Output = Vec2;
    #[inline]
    fn div(self, other: f64) -> Vec2 {
        self * other.recip()
    }
}

impl DivAssign<f64> for Vec2 {
    #[inline]
    fn div_assign(&mut self, other: f64) {
        *self *= other.recip();
    }
}

impl Neg for Vec2 {
    type Output = Vec2;
    #[inline]
    fn neg(self) -> Vec2 {
        Vec2 {
            x: -self.x,
            y: -self.y,
        }
    }
}

fn parse_bool(s: &str) -> bool {
    match s {
        "true" | "True" | "1" => true,
        "false" | "False" | "0" => false,
        _ => panic!("Invalid boolean string"),
    }
}