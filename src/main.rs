use rand::prelude::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

type Spin = i8;
type Lattice = Vec<Spin>;

/// 不同维度的模拟参数
#[derive(Debug, Clone)]
struct DimensionParams {
    dimension: usize,
    l: usize,              // 每个维度的长度
    n: usize,              // 总自旋数
    thermalization: usize, // 热化 MCS
    measurement: usize,    // 测量 MCS
    temp_min: f64,         // 最小温度
    temp_max: f64,         // 最大温度
    dt: f64,               // 温度步长
    beta_star_min: f64,    // β*最小值
    beta_star_max: f64,    // β*最大值
    beta_star_step: f64,   // β*步长
}

impl DimensionParams {
    fn new(dimension: usize) -> Self {
        match dimension {
            1 => Self {
                dimension,
                l: 10000,  // 1D需要很长链
                n: 10000,
                thermalization: 1000,
                measurement: 5000,
                temp_min: 0.5,
                temp_max: 2.0,
                dt: 0.1,
                beta_star_min: 0.05,
                beta_star_max: 0.15,
                beta_star_step: 0.01,
            },
            2 => Self {
                dimension,
                l: 64,
                n: 64 * 64,
                thermalization: 5000,
                measurement: 20000,
                temp_min: 2.0,
                temp_max: 2.4,
                dt: 0.02,
                beta_star_min: 0.10,
                beta_star_max: 0.15,
                beta_star_step: 0.0025,
            },
            3 => Self {
                dimension,
                l: 32,
                n: 32 * 32 * 32,
                thermalization: 3000,
                measurement: 10000,
                temp_min: 4.0,
                temp_max: 5.0,
                dt: 0.05,
                beta_star_min: 0.25,
                beta_star_max: 0.40,
                beta_star_step: 0.01,
            },
            4 => Self {
                dimension,
                l: 16,     // 4D使用更小的晶格
                n: 16 * 16 * 16 * 16,
                thermalization: 2000,
                measurement: 8000,
                temp_min: 6.5,
                temp_max: 7.5,
                dt: 0.1,
                beta_star_min: 0.40,
                beta_star_max: 0.55,
                beta_star_step: 0.01,
            },
            _ => panic!("只支持1-4维"),
        }
    }
}

/// 生成随机初始位形：全部向上
fn init_lattice(n: usize) -> Lattice {
    vec![1; n]
}

/// 计算1D最近邻能量变化
fn delta_e_1d(lat: &Lattice, i: usize, l: usize) -> f64 {
    let left = if i == 0 { l - 1 } else { i - 1 };
    let right = if i == l - 1 { 0 } else { i + 1 };
    let s = lat[i] as f64;
    let nb = lat[left] + lat[right];
    2.0 * s * (nb as f64)
}

/// 计算2D最近邻能量变化
fn delta_e_2d(lat: &Lattice, i: usize, l: usize) -> f64 {
    let row = i / l;
    let col = i % l;
    let up = if row == 0 { l - 1 } else { row - 1 };
    let down = if row == l - 1 { 0 } else { row + 1 };
    let left = if col == 0 { l - 1 } else { col - 1 };
    let right = if col == l - 1 { 0 } else { col + 1 };
    let s = lat[i] as f64;
    let nb = lat[up * l + col] + lat[down * l + col]
           + lat[row * l + left] + lat[row * l + right];
    2.0 * s * (nb as f64)
}

/// 计算3D最近邻能量变化
fn delta_e_3d(lat: &Lattice, i: usize, l: usize) -> f64 {
    let z = i / (l * l);
    let yz = i % (l * l);
    let y = yz / l;
    let x = yz % l;
    
    let x_up = if x == l - 1 { 0 } else { x + 1 };
    let x_down = if x == 0 { l - 1 } else { x - 1 };
    let y_up = if y == l - 1 { 0 } else { y + 1 };
    let y_down = if y == 0 { l - 1 } else { y - 1 };
    let z_up = if z == l - 1 { 0 } else { z + 1 };
    let z_down = if z == 0 { l - 1 } else { z - 1 };
    
    let s = lat[i] as f64;
    let nb = lat[z * l * l + y * l + x_up] + lat[z * l * l + y * l + x_down]
           + lat[z * l * l + y_up * l + x] + lat[z * l * l + y_down * l + x]
           + lat[z_up * l * l + y * l + x] + lat[z_down * l * l + y * l + x];
    2.0 * s * (nb as f64)
}

/// 计算4D最近邻能量变化
fn delta_e_4d(lat: &Lattice, i: usize, l: usize) -> f64 {
    let w = i / (l * l * l);
    let xyz = i % (l * l * l);
    let z = xyz / (l * l);
    let yz = xyz % (l * l);
    let y = yz / l;
    let x = yz % l;
    
    // 4D的8个最近邻
    let x_up = if x == l - 1 { 0 } else { x + 1 };
    let x_down = if x == 0 { l - 1 } else { x - 1 };
    let y_up = if y == l - 1 { 0 } else { y + 1 };
    let y_down = if y == 0 { l - 1 } else { y - 1 };
    let z_up = if z == l - 1 { 0 } else { z + 1 };
    let z_down = if z == 0 { l - 1 } else { z - 1 };
    let w_up = if w == l - 1 { 0 } else { w + 1 };
    let w_down = if w == 0 { l - 1 } else { w - 1 };
    
    let s = lat[i] as f64;
    let nb = lat[w * l * l * l + z * l * l + y * l + x_up] + 
           lat[w * l * l * l + z * l * l + y * l + x_down] +
           lat[w * l * l * l + z * l * l + y_up * l + x] + 
           lat[w * l * l * l + z * l * l + y_down * l + x] +
           lat[w * l * l * l + z_up * l * l + y * l + x] + 
           lat[w * l * l * l + z_down * l * l + y * l + x] +
           lat[w_up * l * l * l + z * l * l + y * l + x] + 
           lat[w_down * l * l * l + z * l * l + y * l + x];
    2.0 * s * (nb as f64)
}

/// 通用能量变化计算函数
fn delta_e(lat: &Lattice, i: usize, params: &DimensionParams) -> f64 {
    match params.dimension {
        1 => delta_e_1d(lat, i, params.l),
        2 => delta_e_2d(lat, i, params.l),
        3 => delta_e_3d(lat, i, params.l),
        4 => delta_e_4d(lat, i, params.l),
        _ => panic!("不支持的维度"),
    }
}

/// 单温度模拟，返回平均 |M|/N
fn run_temp(temp: f64, params: &DimensionParams) -> f64 {
    let mut rng = thread_rng();
    let mut lat = init_lattice(params.n);
    let beta = 1.0 / temp;

    // 热化
    for _ in 0..params.thermalization {
        for _ in 0..params.n {
            let idx = rng.gen_range(0..params.n);
            let de = delta_e(&lat, idx, params);
            if de <= 0.0 || rng.gen_range(0.0..1.0) < (-beta * de).exp() {
                lat[idx] *= -1;
            }
        }
    }

    // 测量
    let mut m_sum = 0.0;
    let mut count = 0;
    for _ in 0..params.measurement {
        for _ in 0..params.n {
            let idx = rng.gen_range(0..params.n);
            let de = delta_e(&lat, idx, params);
            if de <= 0.0 || rng.gen_range(0.0..1.0) < (-beta * de).exp() {
                lat[idx] *= -1;
            }
        }
        let m: i32 = lat.iter().map(|&s| s as i32).sum();
        m_sum += (m.abs() as f64) / (params.n as f64);
        count += 1;
    }
    m_sum / (count as f64)
}

fn main() {
    let dimensions = vec![1, 2, 3, 4];
    
    // 为每个维度创建输出文件
    let mut m_vs_t_files = Vec::new();
    let mut m1_beta_star_files = Vec::new();
    let mut all_params = Vec::new();
    
    for &dim in &dimensions {
        let params = DimensionParams::new(dim);
        all_params.push(params.clone());
        
        // 创建文件名
        let m_vs_t_name = format!("m_vs_t_{}d.dat", dim);
        let m1_beta_star_name = format!("m1_beta_star_{}d.dat", dim);
        
        m_vs_t_files.push(BufWriter::new(File::create(&m_vs_t_name).unwrap()));
        m1_beta_star_files.push(BufWriter::new(File::create(&m1_beta_star_name).unwrap()));
        
        println!("开始计算{}D伊辛模型...", dim);
    }
    
    // 并行计算每个维度
    let results: Vec<_> = dimensions.par_iter().zip(all_params.par_iter())
        .map(|(&dim, params)| {
            let temps: Vec<f64> = float_range(params.temp_min, params.temp_max, params.dt);
            let m_avg: Vec<f64> = temps.par_iter().map(|&t| run_temp(t, params)).collect();
            (dim, params.clone(), temps, m_avg)
        }).collect();
    
    // 写入结果
    for (i, (dim, params, temps, m_avg)) in results.into_iter().enumerate() {
        // 写入 m_vs_t 文件
        writeln!(m_vs_t_files[i], "# T  |M|/N").unwrap();
        for (&t, &m) in temps.iter().zip(&m_avg) {
            writeln!(m_vs_t_files[i], "{:.4} {:.6}", t, m).unwrap();
        }
        
        // 写入 m1_beta_star 文件
        let beta_star_range = float_range(params.beta_star_min, params.beta_star_max, params.beta_star_step);
        write!(m1_beta_star_files[i], "# T").unwrap();
        for &b in &beta_star_range {
            write!(m1_beta_star_files[i], "  M^1/β*={:.4}", b).unwrap();
        }
        writeln!(m1_beta_star_files[i]).unwrap();
        
        for (j, &t) in temps.iter().enumerate() {
            write!(m1_beta_star_files[i], "{:.4}", t).unwrap();
            let m = m_avg[j];
            for &b in &beta_star_range {
                write!(m1_beta_star_files[i], " {:.6}", m.powf(1.0 / b)).unwrap();
            }
            writeln!(m1_beta_star_files[i]).unwrap();
        }
        
        println!("{}D计算完成！数据已写入 m_vs_t_{}d.dat 与 m1_beta_star_{}d.dat", dim, dim, dim);
    }
    
    println!("\n所有维度计算完成！");
    println!("理论β值：1D=0, 2D=0.125, 3D≈0.326, 4D=0.5（平均场论）");
    println!("用Python脚本查看结果：python3 plot_comparison.py");
}

/// 辅助：生成均匀浮点区间
fn float_range(start: f64, end: f64, step: f64) -> Vec<f64> {
    let n = ((end - start) / step).round() as usize + 1;
    (0..n).map(|i| start + (i as f64) * step).collect()
}