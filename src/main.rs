use std::fmt::Error;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::io::prelude::*;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn sub_fac(n: u64) -> u64 {
	if n < 2 {
		return 1 as u64;
	} else {
		return n * sub_fac(n - 1);
	}
}

/* 
fn sub_binom(n: u32, k: u32) -> f64  {
	let mut bin: f = 0.0;
	if n < k {
		return bin;
	} else {
		if k < n/2 {
			for i in (n-k+1)..n {
				bin *=  i as f64 / (i - n + k) as f64;
			}
		} else {
			for i in (k+1)..n {
				bin *=  i as f64 / (i - k) as f64;
			}
		}
		return bin;
	}
}
*/

fn tajd_cstes(nb: u64) -> (f64, f64) {
	let mut a1: f64 = 0.0;
	let mut a2: f64 = 0.0;
	for n in 1..nb {
		a1 += 1.0 / n as f64;
		a2 += 1.0 / (n.pow(2) as f64);
	}

	let b1: f64 = (nb + 1) as f64 / ((3 * (nb - 1)) as f64);
	let b2: f64 = ((2 * (nb.pow(2) + nb + 3)) as f64 )/ ((9 * nb * (nb - 1)) as f64);
	let c1: f64 = b1 - (1.0 / a1);
	let c2: f64 = b2 - ((nb + 2) as f64 / (a1 * nb as f64)) + (a2 / a1.powi(2));
	let e1: f64 = c1 / a1;
	let e2: f64 = c2 / (a1.powi(2) + a2);
	return (e1, e2);
}

/* 
fn nwsfs (pi: Vec<Vec<f64>>, nb: u32) {
	let mut nwsfs: Vec<f64> = Vec::new();
	for ii in 0..nb {
		let mut nwsfs_ii: f64 = 0.0;
		for jj in 1..nb - 1 {
			let multiplier = pi[jj][ii];
			nwsfs_ii += sfsbins(nb)[nb - 2][jj][ii] * multiplier;
		}
		nwsfs.push(nwsfs_ii);
	}
	return nwsfs;
}

fn sfsbins (nb: u32) {
	let mut sfsbins: Vec<Vec<Vec<f64>>> = Vec::new();
	for spl in 2..nb {
		let mut sfsbins_spl: Vec<Vec<f64>> = Vec::new();
		for u in 1..nb - 1 {
			let mut sfsbins_spl_u: Vec<f64> = Vec::new();
			for v in 0..spl {
				let binom = sub_binom(nb - u, spl - v) * sub_binom(u, v);
				let binom = binom / sub_binom(nb, spl);
				sfsbins_spl_u.push(binom);
			}
			sfsbins_spl.push(sfsbins_spl_u);
		}
		sfsbins.push(sfsbins_spl);
	}
	return sfsbins;
}
*/

fn main() {
  let all_nb: u64 = 242824;
  let min_spl: u64 = 27428;
  let mut pi_tmp: f64 = 0.0;
  let mut h_tmp: f64 = 0.0;
  let mut seg: f64 = 0.0;
  let mut tajd: f64 = 0.0;
  let det: f64 = 8276.0;

  // calculate a_1
  let mut a_1: f64 = 0.0;
  for m in 1..(min_spl+1) {
    a_1 += 1.0 / m as f64;
  }
 /*  let mut a_1: Vec<f64> = Vec::new();
  for m in 1..all_nb {
	let mut a1tmp: f64 = 0.0;
	for n in 1..m {
		a1tmp += 1.0 / n as f64;
	}
	a_1.push(a1tmp);
  }
*/
    // calculate sfsbins
	let base: f64 = (242824.0 - 27428.0) / 242824.0;
	let mut sfsbins: Vec<Vec<f64>> = vec![vec![0.0; (min_spl+1) as usize]; all_nb as usize];
	let mut binom: f64 = base;
	sfsbins[1][0] = base;
	for u in 2..(all_nb-min_spl+1) {
	  binom = binom * (242824.0 - 27428.0 + 1.0 - u as f64) / (242824.0 + 1.0 - u as f64) ;
	  sfsbins[u as usize][0] = binom;
	}
	sfsbins[2][1] = 2.0 * 27428.0 * ( 242824.0 - 27428.0) / 242824.0 / 242823.0;
	binom = sfsbins[2][1];
	println!("sfsbins[2][1]:\t{}", sfsbins[2][1]);
	for u in 3..(all_nb-min_spl) {
	  binom = binom * u as f64 * ( 242824.0 - 27428.0 + 1.0 - u as f64) / (242824.0 + 1.0 - u as f64) / (u as f64 - 1.0);
	  sfsbins[u as usize][1] = binom;
	  for v in 2..u {
		if v <= min_spl {
		  sfsbins[u as usize][v as usize] = sfsbins[u as usize][(v-1) as usize] * (u as f64 - v as f64 + 1.0) * (27428.0 - v as f64 + 1.0) / (242824.0 - 27428.0 + v as f64 -u as f64) / v as f64;
		}
	  }
	}
  
	for u in (all_nb-min_spl)..all_nb {
	  for v in 2..u {
		if v <= min_spl {
		  if (u - v) < (242824 - 27428) {
			sfsbins[u as usize][v as usize] = sfsbins[(u-1) as usize][v as usize] * u as f64 * (242824.0 - 27428.0 - u as f64 + v as f64) / (242824.0 + 1.0 - u as f64) / (u as f64 - v as f64);
		  }
		  else {
			sfsbins[u as usize][v as usize] = sfsbins[(u-1) as usize][v as usize] * u as f64 / (242824.0 + 1.0 - u as f64) / (u as f64 - v as f64);
		  }
  
		}
	  }
	}

	let mut file = File::create("foo.txt").expect("failed to open file");
	for (i, val) in sfsbins.iter().enumerate() {
		for (j, val2) in val.iter().enumerate() {
			if (i as u64) > (j as u64) {
				file.write_all(format!("{}\t{}\t{}\n", i, j, val2).as_bytes()).expect("failed to write to file");
			}
		}
	}
	//Ok(());


/*	let mut sfsbins: Vec<Vec<f64>> = vec![vec![0.0; all_nb as usize]; (min_spl+1) as usize];
	for y in 1..all_nb - 1 {
		for z in 0..min_spl {
			let binom: f64 = sub_binom(all_nb - y, min_spl - z) * sub_binom(y, z);
			let binom = binom / sub_binom(all_nb, min_spl);
			sfsbins[y as usize][z as usize] = binom;
		}
	}
*/
  // Read in pi values
  let mut nwsfs: Vec<f64> = vec![0.0; (min_spl+1) as usize ];
if let Ok(lines) = read_lines("./pykeys_new.txt") {
	for line in lines {
		if let Ok(ip) = line {
			let mut spl = ip.split_whitespace();
			let snp: u8 = spl.next().unwrap().parse::<u8>().unwrap();
			let nb = spl.next().unwrap().parse::<u32>().unwrap();
			let freq = spl.next().unwrap().parse::<u32>().unwrap();
			for i in 0..(min_spl+1){
				nwsfs[i as usize] += sfsbins[freq as usize][i as usize] * snp as f64;
			}
			println!("nwsfs{{{}}}: {} {}", freq, nwsfs[0], nwsfs[27428]);
		}
	}
}
  /* 
  // calculate nwsfs
  let mut nwsfs: Vec<f64> = vec![0.0; (min_spl+1) as usize ];
  for spl in 0..(pi_all.len()-1) {
	for jj in 1..all_nb - 1 {
	  for ii in 0..min_spl {	
		let multiplier = pi_all[spl as usize][jj as usize];
		nwsfs[ii as usize] += sfsbins[jj as usize][ii as usize] * multiplier as f64;
	  }
	}
  }	
  */


  let tipi: u64 = min_spl;
  for y in 1..tipi {
	pi_tmp += y as f64 * (tipi as f64 - y as f64) * nwsfs[y as usize];
	h_tmp += y.pow(2) as f64 * nwsfs[y as usize];
	seg += nwsfs[y as usize];
  }

  let pi = pi_tmp / ((tipi as f64 * (tipi as f64 - 1.0)) / 2.0) / det ;
  let pi_corrected = 0.75 * (1.0 - ( 4.8 * pi / 3.0)).ln();
  let theta = seg / a_1 / det;
  let thetaW_corrected = ( det / (det - seg as f64)).ln() / a_1;

  let mut d: f64 = 0.0;
  let (e1, e2) = tajd_cstes(min_spl);
  let p1p2: f64 = e1 * seg + e2 * seg * (seg - 1.0);
  if pi > 0.0 && theta > 0.0 {
	d = pi - theta;
  }

  if p1p2 > 0.0 {
	tajd = d * det / p1p2.sqrt();
  }

  let hache = h_tmp / ((tipi as f64 * (tipi as f64 - 1.0)) / 2.0);
  let final_h = pi - hache;

  println!("pi:\t{}", pi);
  println!("pi Corrected:\t{}", pi_corrected);
  println!("thetaW:\t{}", theta);
  println!("thetaW Corrected:\t{}", thetaW_corrected);
  println!("D:\t{}", d);
  println!("Tajima's D:\t{}", tajd);
  println!("H:\t{}", final_h);
}



