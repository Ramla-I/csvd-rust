#![no_std]
#![feature(alloc)]
#![feature(extern_crate_item_prelude)]

extern crate alloc;
extern crate num_complex;
extern crate libm;
// extern crate rand;

mod csvd;
mod test;

use num_complex::Complex32;
use alloc::vec::Vec;
use self::csvd::csvd;

/// Finds the pseudo-inverse of matrix using Singular Value Decomposition
/// Assumes that input_mat has dimensions mxn and inverse_mat has dimension nxm
/// Stores the return value in inverse_mat, and values of input_mat are modified
pub fn pinv(mut input_mat: &mut Vec<Complex32>, mut inverse_mat: &mut Vec<Complex32>, input_num_rows: usize, input_num_cols: usize) -> Result< (), &'static str> {
    let m = input_num_rows;
    let n = input_num_cols;

    //create S vector with dimension n
    let mut s: Vec<f32> = Vec::with_capacity(n);
    for _ in 0..n {
        s.push(0.0);
    }

    //create U matrix dimension mxm
    let mut u: Vec<Complex32> = Vec::with_capacity(m*m);
    for _ in 0..m*m {
        u.push(Complex32{re: 0.0, im: 0.0});
    }

    //create v matrix with dimension nxn
    let mut v: Vec<Complex32> = Vec::with_capacity(n*n);
    for _ in 0..n*n {
        v.push(Complex32{re: 0.0, im: 0.0});
    }

    csvd(&mut input_mat, m, n, n, m, 0, m, n, &mut s, &mut u, &mut v)?;

    find_pinv_from_svd(&mut s, &u, &v, m, n, &mut inverse_mat);


    Ok(())
}

/// Finds the pseudo-inverse of a matrix from the singular value decompositions
/// INV = V x S+ x U*
/// where S+ is found by taking the reciprocal fo all non-zero elements of S and changing the dimension from n to nxm
/// and U* is the conjugate-transpose of U
fn find_pinv_from_svd(s: &mut Vec<f32>, u: &Vec<Complex32>, v: &Vec<Complex32>, m: usize, n: usize, inv: &mut Vec<Complex32>) {
    // cut-off value for a number to be assumed to be 0
    let eps = 0.0001;
    let mut n_ = n;

    // take reciprocal of all non-zero elements in S
    for i in 0..n {
        if s[i] > eps {
            s[i] = 1.0/s[i];
        }
        else {
            s[i] = 0.0;
        }
    }

    // extend S to be of size m 
    while n_ < m {
        s.push(0.0);
        n_ += 1;
    }

    for i in 0..n {
        for j in 0..m {
            inv[i*n + j].re = 0.0;
            inv[i*n + j].im = 0.0;
            for k in 0..n {
                inv[i*n + j] = inv[i*n + j] + v[i*n + k] * s[k] * u[j*m + k].conj();
            }
        }
    }
}

pub fn matrix_mult(mat_a: &[Complex32], a_rows: usize, a_cols: usize, mat_b: &[Complex32], b_rows: usize, b_cols: usize, mat_c: &mut[Complex32]) -> Result< (), &'static str> {
    if a_cols != b_rows {
        return Err("Matrix dimension not compatible!");
    }

    for i in 0..a_rows {
        for j in 0..b_cols {
            for k in 0..a_cols{
                mat_c[i*a_rows + j] += mat_a[i*a_rows + k] * mat_b[k*b_rows + j]; 
            }
        }
    }

    Ok(())
}

