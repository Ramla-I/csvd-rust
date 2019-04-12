#![no_std]
#![feature(alloc)]
#![feature(extern_crate_item_prelude)]

#[macro_use] extern crate log;
extern crate alloc;
extern crate num_complex;
extern crate libm;
// extern crate rand;
// extern crate aligned_vec;

pub mod csvd;
pub mod test;

use num_complex::Complex32;
use alloc::vec::Vec;
use self::csvd::csvd;
// use aligned_vec::{aligned_alloc, aligned_alloc_f32_16};
use core::mem;

#[repr(align(16))]
struct Align16(u64,u64);

#[repr(align(32))]
struct Align32(u64,u64,u64,u64);

#[repr(align(64))]
struct Align64(u64,u64,u64,u64,u64,u64,u64,u64);

/// function to return a vector aligned at "alignment" no of bytes with the only valid inputs being 16, 32 and 64
/// len is the the number of Complex32 items in the vector
/// vec is where the resulting vector is stored
/// function also initializes all elements to '0.0 + 0.0i'
pub fn aligned_alloc(alignment: u8, len: usize, mut vec: &mut Vec<Complex32>) -> Result<(), &'static str> {
    
    if alignment == 16 {
        aligned_alloc_16(len, &mut vec);
    }
    else if alignment == 32 {
        aligned_alloc_32(len, &mut vec);
    }
    else if alignment == 64 {
        aligned_alloc_64(len, &mut vec);
    }
    else{
        return Err("Invalid alignment input")
    }
    
 
    Ok(())
}

fn aligned_alloc_16(len: usize, vec: &mut Vec<Complex32>) {
    let C32_in_A16 = 2;

    let buffer : Vec<Align16> = Vec::with_capacity(len/C32_in_A16);
    let buffer_ptr = buffer.as_slice().as_ptr() as usize;
    // debug!("Aligned pointer: {:#X}", buffer_ptr);
    let ptr = buffer_ptr as *mut Complex32;
    mem::forget(buffer);
    *vec = unsafe {Vec::from_raw_parts(ptr, len, len)};

    for num in vec {
        num.re = 0.0;
        num.im = 0.0;
    }
}

fn aligned_alloc_32(len: usize, vec: &mut Vec<Complex32>) {
    let C32_in_A32 = 4;

    let buffer : Vec<Align32> = Vec::with_capacity(len/C32_in_A32);
    let buffer_ptr = buffer.as_slice().as_ptr() as usize;
    // debug!("Aligned pointer: {:#X}", buffer_ptr);
    let ptr = buffer_ptr as *mut Complex32;
    mem::forget(buffer);
    *vec = unsafe {Vec::from_raw_parts(ptr, len, len)};

    for num in vec {
        num.re = 0.0;
        num.im = 0.0;
    }
}

fn aligned_alloc_64(len: usize, vec: &mut Vec<Complex32>) {
    let C32_in_A64 = 8;

    let buffer : Vec<Align64> = Vec::with_capacity(len/C32_in_A64);
    let buffer_ptr = buffer.as_slice().as_ptr() as usize;
    // debug!("Aligned pointer: {:#X}", buffer_ptr);
    let ptr = buffer_ptr as *mut Complex32;
    mem::forget(buffer);
    *vec = unsafe {Vec::from_raw_parts(ptr, len, len)};

    for num in vec {
        num.re = 0.0;
        num.im = 0.0;
    }
    
}

///TODO: need to generalize this
pub fn aligned_alloc_u8_32(len: usize, vec: &mut Vec<u8>) {
    let U8_in_A32 = 32;

    let buffer : Vec<Align16> = Vec::with_capacity(len/U8_in_A32);
    let buffer_ptr = buffer.as_slice().as_ptr() as usize;
    // debug!("Aligned pointer: {:#X}", buffer_ptr);
    let ptr = buffer_ptr as *mut u8;
    mem::forget(buffer);
    *vec = unsafe {Vec::from_raw_parts(ptr, len, len)};

    for num in vec {
        *num = 0;
    }
}

///TODO: need to generalize this
pub fn aligned_alloc_f32_32(len: usize, vec: &mut Vec<f32>) {
    let F32_in_A32 = 8;

    let buffer : Vec<Align16> = Vec::with_capacity(len/F32_in_A32);
    let buffer_ptr = buffer.as_slice().as_ptr() as usize;
    // debug!("Aligned pointer: {:#X}", buffer_ptr);
    let ptr = buffer_ptr as *mut f32;
    mem::forget(buffer);
    *vec = unsafe {Vec::from_raw_parts(ptr, len, len)};

    for num in vec {
        *num = 0.0;
    }
}

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
pub fn find_pinv_from_svd(s: &mut Vec<f32>, u: &Vec<Complex32>, v: &Vec<Complex32>, m: usize, n: usize, inv: &mut Vec<Complex32>) {

    // debug!("In find pinv from svd");
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
    let a = &mat_a[0..a_rows*a_cols];
    let b = &mat_b[0..b_rows*b_cols];
    let c = &mut mat_c[0..a_rows*b_cols];
    if a_cols != b_rows {
        return Err("Matrix dimension not compatible!");
    }

    // //transpose b
    // let mut b = Vec::new();
    // aligned_alloc_32(b_rows*b_cols, &mut b);

    // for i in 0..b_cols {
    //     for j in 0..b_rows {
    //         b[i*b_rows + j] = mat_b[j*b_cols +i];
    //     }
    // }

    // let b = &b[0..b_cols*b_rows];


    // let mut i = 0;
    // loop {
    //     if i==a_rows {break;}
    //     let mut j = 0;

    //     loop {
    //         if j == a_rows { break;}
    //         let mut k = 0;

    //         loop {
    //             if k == a_cols { break;}

    //             c[i * b_cols + j] += a[j*a_cols + k] * b[j*a_cols + k];
    //             k += 1;
    //         }

    //         j += 1;
    //     }

    //     i +=1;
    // }
    // const a_r: usize = 8;
    // const b_c: usize = 8;
    // const a_c: usize = 8;

    for i in 0..a_rows {
        for j in 0..b_cols {
            for k in 0..a_cols{
                c[i * b_cols + j] += a[i*a_cols + k] * b[k*b_cols + j]; 
            }
        }
    }

    Ok(())
}

