use num_complex::Complex32;
use libm::F32Ext;
use alloc::vec::Vec;
// use rand::Rng;

use super::csvd::csvd;
use super::pinv;

// fn print_matrix(mat: &Vec<Complex32>, rows: usize, cols: usize) {
//     for i in 0..rows {
//         for j in 0..cols{
//             print!("{} + {}i, ", mat[i*rows + j].re, mat[i*rows + j].im);
//         }
//         println!("");
//     }
//     println!("");
// }

// fn print_vector(mat: &Vec<f32>, rows: usize) {
//     for i in 0..rows {
//         println!("{}", mat[i]);
//     }
//     println!("");
// }

/// Finds the original matrix from the singular value decompositions
/// A = U x S x V*
/// stores the new matrix in a
fn find_orig_matrix_from_svd(mut a: &mut Vec<Complex32>, m: usize, n: usize) {
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

    let _ = csvd(&mut a, m, n, n, m, 0, m, n, &mut s, &mut u, &mut v);

    let min = m.min(n);
    for i in 0..m {
        for j in 0..n {
            a[i*m + j].re = 0.0;
            a[i*m + j].im = 0.0;
            for k in 0..min {
                a[i*m + j] = a[i*m + j] + u[i*m + k] * s[k] * v[j*n + k].conj();
            }
        }
    }
}
/// Verifies pinv function
/// Checks if A*Ainv*A = A where A is a mxn matrix
fn check_pinv(mut a: &mut Vec<Complex32>, m: usize, n: usize) -> bool {  
    let a_orig = a.clone();

    //create inverse matrix with dimension nxm
    let mut inv: Vec<Complex32> = Vec::with_capacity(n*m);
    for _ in 0..n*m {
        inv.push(Complex32{re: 0.0, im: 0.0});
    }

    let _ = pinv(&mut a, &mut inv, m, n);
    
    //create I matrix dimension mxm
    let mut I: Vec<Complex32> = Vec::with_capacity(m*m);
    for _ in 0..m*m {
        I.push(Complex32{re: 0.0, im: 0.0});
    }

    // I = A x Ainv
    for i in 0..m {
        for j in 0..m {
            for k in 0..n{
                I[i*m + j] = I[i*m + j] + (a_orig[i*m + k] * inv[k*n + j]); 
            }
        }
    }

    // I x A
    for i in 0..m {
        for j in 0..n {
            a[i*m + j].re = 0.0;
            a[i*m + j].im = 0.0;
            for k in 0..m{
                a[i*m + j] = a[i*m + j] + (I[i*m +k] * a_orig[k*n + j]); 
            }
        }
    }

    check_matrix_equality(&a_orig, &a, m, n)

}

/// checks that 2 complex matrices are equal by taking the square of the euclidean distance between the elements
fn check_matrix_equality(a: &Vec<Complex32>, b: &Vec<Complex32>, m: usize, n:usize) -> bool {
    let mut equal = true;

    let eps = 0.0001;

    for i in 0..m {
        for j in 0..n {
            if F32Ext::powf(a[i*m + j].re - b[i*m + j].re, 2.0) + F32Ext::powf(a[i*m + j].im - b[i*m + j].im, 2.0) > eps {
                equal = false;
            }
        }
    }

    equal
}

/// Checks the correctness of svd function in 2 ways
/// 1. multiplies decomposed matrices together to see if equal to original matrix 
/// 2. finds inverse of matrix using svd and then verifies the correctness of the inverse
/// a has dimensions m x n
fn check_svd(mut a: &mut Vec<Complex32>, m: usize, n: usize) {
    
    let a_orig  = a.clone(); 

    find_orig_matrix_from_svd(&mut a, m, n);

    // if check_matrix_equality(&a_orig, &a, m, n){
    //     println!("svd successful");
    // }
    // else {
    //     println!("svd failed");
    // }

    // *a = a_orig.clone();

    // if check_pinv(&mut a, m, n) {
    //     println!("pseudo-inverse successful");
    // }

    // else {
    //     println!("pseudo-inverse failed");
    // }
    
}

// /// A basic example to test with: https://math.stackexchange.com/questions/647321/moore-penrose-inverse-of-complex-square-matrices
// pub fn test() {

//     let m = 8;
//     let n = 8;

//     let mut a: Vec<Complex32> = Vec::with_capacity(m*n);

//     let mut rng = rand::thread_rng();

//     for _ in 0..m*n {
//         a.push(Complex32{re: rng.gen(), im: rng.gen()})
//     }
    
//     check_svd(&mut a, m, n) ;
  
// }