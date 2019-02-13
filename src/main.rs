// #![no_std]

/// CSVD computes the singular value decomposition of an M by N complex matrix.
///
/// Discussion:
///
///    This routine requires that N <= M.
///
///    The singular value decomposition of a complex M by N matrix A
///    has the form
///
///      A = U S V*
///
///    where 
///
///      U is an M by M unitary matrix,
///      S is an M by N diagonal matrix,
///      V is an N by N unitary matrix.
///
///    Moreover, the entries of S are nonnegative and occur on the diagonal
///    in descending order.
///
///    Several internal arrays are dimensioned under the assumption
///    that N <= 100.

///
///  Reference:
///
///    Peter Businger, Gene Golub,
///    Algorithm 358:
///    Singular Value Decomposition of a Complex Matrix,
///    Communications of the ACM,
///    Volume 12, Number 10, October 1969, pages 564-565.
///
///  Parameters:
///
///   Input/output, complex A(MMAX,*), the M by N matrix, which may be
///    augmented by P extra columns to which the transformation U*
///    is to be applied.  On output, A has been overwritten, and
///    if 0 < P, columns N+1 through N+P have been premultiplied by U*.
///
///    Input, integer MMAX, the leading dimension of the arrays A
///    and U.
///
///    Input, integer NMAX, the leading dimension of V, and perhaps
///    the second dimension of A and U.
///
///    Input, integer M, N, the number of rows and columns in A.
///    It must be the case that 1 <= N <= M.  Several internal arrays are
///    dimensioned under the assumption that N <= NBIG, where NBIG
///    is an internal parameter, currently set to 100.
///
///    Input, integer P, the number of vectors, stored in A(*,N+1:N+P),
///  to which the transformation U* should be applied.
///
///    Input, integer NU, the number of columns of U to compute.
///
///    Input, integer NV, the number of columns of V to compute.
///
///    Output, real S(N), the computed singular values.
///
///    Output, complex U(MMAX,NU), the first NU columns of U.
///
///    Output, complex V(NMAX,NV), the first NV columns of V.
///
///  Local Parameters:
///
///    Local, real ETA, the relative machine precision.
///    The original text uses ETA = 1.5E-8.
///
///    Local, integer NBIG, is a parameter used to dimension work arrays.
///    The size of NBIG limits the maximum possible size of N that can
///    be handled.  If you want to work with values of N that are larger,
///    simply increase the value assigned to NBIG below.
///
///    Local, real TOL, the smallest normalized positive number, divided by ETA.
///    The original test uses TOL = 1.E-31.

extern crate num_complex;
extern crate libm;

use num_complex::Complex32;
use libm::F32Ext;

const NBIG: usize = 100;

fn print_matrix(mat: &Vec<Vec<Complex32>>, rows: usize, cols: usize) {
    for i in 0..rows {
        for j in 0..cols{
            print!("{} + {}i, ", mat[i][j].re, mat[i][j].im);
        }
        println!("");
    }
    println!("");
}

fn print_vector(mat: &Vec<f32>, rows: usize) {
    for i in 0..rows {
        println!("{}", mat[i]);
    }
    println!("");
}

fn cabs(input: &Complex32) -> f32{
    F32Ext::sqrt(F32Ext::powf(input.re, 2.0) + F32Ext::powf(input.im, 2.0))
}

fn csvd(a: &mut Vec<Vec<Complex32>>, mmax: usize, nmax: usize, n: usize, m: usize, p: usize, nu: usize, nv: usize, 
        s: &mut Vec<f32>, u: &mut Vec<Vec<Complex32>>, v: &mut Vec<Vec<Complex32>>) 
        -> Result<(), &'static str> {

    //check n
    if n < 1 {
        return Err("Fatal Error = Input N < 1");
    }
    else if NBIG < n{
        return Err("Fatal Error: NBIG < N");
    }

    //check m
    if m < 1 {
        return Err("Fatal Error: Input M < 1");
    } 
    else if m < n {
        return Err("Fatal Error: M < N");
    }
    
    // Householder reduction.
    let mut c: [f32; NBIG] = [0.0; NBIG];
    c[1] = 0.0;
    let mut k = 0;
    let mut b: [f32; NBIG] = [0.0; NBIG];
    let mut k1;
    let tol = 1.5 * F32Ext::powf(10.0, -31.0);

    //10 continue for k in 0..n
    for k in 0..n {
        k1 = k + 1;

        // Elimination of A(I,K), I = K+1, ..., M.
        let mut z: f32 = 0.0;
        for i in k..m {
            z = z + F32Ext::powf(a[i][k].re, 2.0) + F32Ext::powf(a[i][k].im, 2.0);
        }

        b[k] = 0.0;

        let (mut w, mut q);
        if tol < z {

            z = F32Ext::sqrt(z);
            b[k] = z;
            w = cabs(&a[k][k]);

            if w == 0.0 {
                q = Complex32{ re: 1.0, im: 0.0};
            }
            else {
                q = a[k][k]/w;
            }

            a[k][k] = q * ( z + w );

            if k != (n - 1 + p) {
                for j in k1..(n + p){
                    q = Complex32{ re: 0.0, im: 0.0};
                    
                    for i in k..m {
                        q = q + a[i][k].conj() * a[i][j];
                    }
                    q = q / ( z * ( z + w ) );

                    for i in k..m {
                        a[i][j] = a[i][j] - q * a[i][k]
                    }
                }

                // Phase transformation.
                q = -a[k][k].conj() / cabs(&a[k][k]);

                for j in k1..(n + p) {
                    a[k][j] = q * a[k][j];
                }
            }
        }

        //Elimination of A(K,J), J = K+2, ..., N

        if k == (n - 1) {
            break;
        }

        z = 0.0;
        for j in k1..n {
            z = z + F32Ext::powf(a[k][j].re, 2.0) + F32Ext::powf(a[k][j].im, 2.0);
        }
        c[k1] = 0.0;

        if tol < z {
            z = F32Ext::sqrt(z);
            c[k1] = z;
            w = cabs(&a[k][k1]);

            if w == 0.0 {
                q = Complex32{ re: 1.0, im: 0.0};
            }
            else {
                q = a[k][k1] / w;
            }

            a[k][k1] = q * (z + w);

            for i in k1..m {
                q = Complex32{ re: 0.0, im: 0.0};

                for j in k1..n {
                    q = q + a[k][j].conj()  * a[i][j];
                }

                q = q / (z * (z + w));

                for j in k1..n {
                    a[i][j] = a[i][j] - q * a[k][j];
                }
            }
    
            // Phase transformation.
            q = -a[k][k1].conj() / cabs(&a[k][k1]);
            for i in k1..m {
                a[i][k1] = a[i][k1] * q;
            }
        }
    }

    // Tolerance for negligible elements.
    //140 continue
    let mut eps: f32 = 0.0;
    let eta: f32 = 1.1920929 * F32Ext::powf(10.0, -7.0);
    let mut t: [f32; NBIG] = [0.0; NBIG];

    for k in 0..n {
       s[k] = b[k];
       t[k] = c[k];
       eps = eps.max(s[k] + t[k]);
    }

    eps = eps * eta;

    // Initialization of U and V.
    if 0 < nu {
        for j in 0..nu {
            for i in 0..m {
                u[i][j] = Complex32{re: 0.0, im: 0.0};
            }
            u[j][j] = Complex32{re: 1.0, im: 0.0};
        }
    }

    if 0 < nv {
        for j in 0..nv {
            for i in 0..n {
                v[i][j] = Complex32{re: 0.0, im: 0.0};
            }
            v[j][j] = Complex32{re: 1.0, im: 0.0};
        }
    }

    // println!("****************");
    
    // println!("a");
    // print_matrix(a, 3, 3);

    // println!("u");
    // print_matrix(u, 3, 3);

    // println!("v");
    // print_matrix(v, 3, 3);

    // println!("****************");

    let mut l = 0;
    let mut cs;
    let mut sn;
    let mut l1;
    let mut f;
    let mut h;
    let mut w;
    let mut x;
    let mut y;
    let mut q;
    // let mut r;
    let mut g;

    // QR diagonalization.
    for kk in 0..n {
        k = n - 1 - kk;

        //Test for split.
        //220 continue
        loop {
            for ll in 0..=k {
                l = k - ll;
                if F32Ext::abs(t[l]) <= eps {
                    //go to 290
                    break;
                }

                if F32Ext::abs(s[l-1]) <= eps {
                    //go to 240
                    break;
                }

            }

            if F32Ext::abs(t[l]) <= eps {
                //go to 290
            }

            //Cancellation of E(L).
            // 240 continue
            else if F32Ext::abs(s[l-1]) <= eps {
                cs = 0.0;
                sn = 1.0;
                l1 = l - 1;

                for i in l..=k {
                    f = sn * t[i];
                    t[i] = cs * t[i];

                    if F32Ext::abs(f) <= eps {
                        //go to 290
                        break;
                    }

                    h = s[i];
                    w = F32Ext::sqrt(f * f + h * h);
                    s[i] = w;
                    cs = h / w;
                    sn = - f / w;

                    if 0 < nu {
                        for j in 0..n {
                            x = u[j][l1].re;
                            y = u[j][i].re;
                            u[j][l1] = Complex32{re: x * cs + y * sn, im: 0.0};
                            u[j][i] = Complex32{re: y * cs - x * sn, im: 0.0};
                        }
                    }

                    // if p != 0 {
                    //     for j in (n + 1)..=(n + p) {
                    //         q = a[l1][j];
                    //         r = a[i][j];
                    //         a[l1][j] = q * cs + r * sn;
                    //         a[i][j] = r * cs - q * sn;
                    //     }
                    // }
                }
            }

            // Test for convergence.
            // 290 continue
            w = s[k];

            if l == k {
                //go to 360
                break;
            }

            // Origin shift.
            x = s[l];
            y = s[k-1];
            g = t[k-1];
            h = t[k];
            f = ( ( y - w ) * ( y + w ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y );
            g = F32Ext::sqrt(f * f + 1.0);
            if f < 0.0 {
                g = -g;
            }
            f = ( ( x - w ) * ( x + w ) + ( y / ( f + g ) - h ) * h ) / x;

            // QR Step
            cs = 1.0;
            sn = 1.0;
            l1 = l + 1;

            for i in l1..=k {

                g = t[i];
                y = s[i];
                h = sn * g;
                g = cs * g;
                w = F32Ext::sqrt(h * h + f * f);
                t[i-1] = w;
                cs = f / w;
                sn = h / w;
                f = x * cs + g * sn;
                g = g * cs - x * sn;
                h = y * sn;
                y = y * cs;

                if 0 < nv {
                    for j in 0..n {
                        x = v[j][i-1].re;
                        w = v[j][i].re;
                        v[j][i-1] = Complex32{re: x * cs + w * sn, im: 0.0};
                        v[j][i] = Complex32{re: w * cs - x * sn, im: 0.0};
                    }
                }

                w = F32Ext::sqrt(h * h + f * f);
                s[i-1] = w;
                cs = f / w;
                sn = h / w;
                f = cs * g + sn * y;
                x = cs * y - sn * g;

                if 0 < nu {
                    for j in 0..n {
                        y = u[j][i-1].re;
                        w = u[j][i].re;
                        u[j][i-1] = Complex32{re: y * cs + w * sn, im: 0.0};
                        u[j][i] = Complex32{re: w * cs - y * sn, im: 0.0};
                    }
                }

                // if p != 0 {
                //     for j in (n + 1)..=(n + p) {
                //         q = a[i-1][j];
                //         r = a[i][j];
                //         a[i-1][j] = q * cs + r * sn;
                //         a[i][j] = r * cs - q * sn;
                //     }
                // }
            }

            t[l] = 0.0;
            t[k] = f;
            s[k] = x;
            //go to 220
        }

        // Convergence
        // 360 continue

        if w < 0.0 {
            s[k] = -w;

            if 0 < nv {
                for j in 0..n {
                    v[j][k] = -v[j][k];
                }
            }
        }
    }

    let mut j;
    
    // Sort the singular values.
    for k in 0..n {
        g = -1.0;
        j = k;

        for i in k..n {
            if g < s[i] { 
                g = s[i];
                j = i;
            }
        }

        if j != k {
            s[j] = s[k];
            s[k] = g;

            //Interchange V(1:N,J) and V(1:N,K).
            if 0 < nv {
               for i in 0..n {
                    q = v[i][j];
                    v[i][j] = v[i][k];
                    v[i][k] = q;
               }
            }

            // Interchange U(1:N,J) and U(1:N,K).
            if 0 < nu {
                for i in 0..n {
                    q = u[i][j];
                    u[i][j] = u[i][k];
                    u[i][k] = q;
                }
            }

            // Interchange A(J,N1:NP) and A(K,N1:NP).
            // if p != 0 {
            //     for i in (n + 1)..=(n + p) {
            //         q = a[j][i];
            //         a[j][i] = a[k][i];
            //         a[k][i] = q;
            //     }
            // }
        }
    }

    // Back transformation.
    if 0 < nu {
        for kk in 0..n {
            k = n - 1 - kk;

            if b[k] != 0.0 {
                q = -a[k][k] / cabs(&a[k][k]);

                for j in 0..nu {
                    u[k][j] = q * u[k][j];
                }

                for j in 0..nu {

                    q = Complex32{re: 0.0, im: 0.0};

                    for i in k..m {
                        q = q + a[i][k].conj() * u[i][j];
                    }

                    q = q / (cabs(&a[k][k]) * b[k]);

                    for i in k..m {
                        u[i][j] = u[i][j] - q * a[i][k];
                    }

                }

            }

        }

    }

    if 0 < nv {

        if 1 < n {

            for kk in 1..n {
                k = n - 1 - kk;
                k1 = k + 1;

                if c[k1] != 0.0 { 
                    q = -(a[k][k1].conj()) / cabs(&a[k][k1]);

                    for j in 0..nv {
                        v[k1][j] = q * v[k1][j];
                    }

                    for j in 0..nv {
                        q = Complex32{re: 0.0, im: 0.0};

                        for i in k1..n {
                            q = q + a[k][i] * v[i][j];
                        }
                        q = q / (cabs(&a[k][k1]) * c[k1]);

                        for i in k1..n {
                            v[i][j] = v[i][j] - q * a[k][i].conj();
                        }
                    }
                }
            }
        }
    }     

    Ok(())   
}

/// Finds the original matrix from the singular value decompositions
/// A = U x S x V*
fn find_orig_matrix_from_svd(s: &Vec<f32>, u: &Vec<Vec<Complex32>>, v: &Vec<Vec<Complex32>>, m: usize, n: usize, a: &mut Vec<Vec<Complex32>>) {
    let min = m.min(n);
    for i in 0..m {
        for j in 0..n {
            a[i][j].re = 0.0;
            a[i][j].im = 0.0;
            for k in 0..min {
                a[i][j] = a[i][j] + u[i][k] * s[k] * v[j][k].conj();
            }
        }
    }
}

/// Finds the pseudo-inverse of a matrix from the singular value decompositions
/// INV = V x S+ x U*
/// where S+ is found by taking the reciprocal fo all non-zero elements of S and changing the dimension from n to nxm
/// and U* is the conjugate-transpose of U
fn find_pinv_from_svd(s: &mut Vec<f32>, u: &Vec<Vec<Complex32>>, v: &Vec<Vec<Complex32>>, m: usize, n: usize, inv: &mut Vec<Vec<Complex32>>) {
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
            inv[i][j].re = 0.0;
            inv[i][j].im = 0.0;
            for k in 0..n {
                inv[i][j] = inv[i][j] + v[i][k] * s[k] * u[j][k].conj();
            }
        }
    }
}



/// Checks the svd function by comparing original matrix with inverse
/// a has dimensions m x n
fn check_svd(mut a: &mut Vec<Vec<Complex32>>, mut inv: &mut Vec<Vec<Complex32>>, m: usize, n: usize) {
    
    let a_orig  = a.clone(); 

    let mut u: Vec<Vec<Complex32>> = vec![
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}]];

    let mut v: Vec<Vec<Complex32>> = vec![
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}]];
    
    let mut I: Vec<Vec<Complex32>> = vec![
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}]];

    let mut s: Vec<f32> = vec![0.0, 0.0, 0.0];

    let _ = csvd(&mut a, m, n, n, m, 0, m, n, &mut s, &mut u, &mut v);
    
    find_orig_matrix_from_svd(&s, &u, &v, m, n, &mut a);

    let mut equal = true;

    let eps = 0.0001;

    for i in 0..m {
        for j in 0..n {
            if F32Ext::powf(a_orig[i][j].re - a[i][j].re, 2.0) + F32Ext::powf(a_orig[i][j].im - a[i][j].im, 2.0) > eps {
                equal = false;
            }
        }
    }

    if equal == false {
        println!("SVD failed!");
    }
    else {
        println!("SVD successful!");
    }

    print_matrix(&a_orig, m, n);
    print_matrix(&a, m, n);
    
    find_pinv_from_svd(&mut s, &u, &v, m, n, &mut inv);

    print_matrix(&inv, n, m);

    // Find A * Ainv * A
    // should be equal to A
    for i in 0..m {
        for j in 0..m {
            for k in 0..n{
                I[i][j] = I[i][j] + (a_orig[i][k] * inv[k][j]); 
            }
        }
    }
    print_matrix(&I, m, m);

    for i in 0..m {
        for j in 0..n {
            a[i][j].re = 0.0;
            a[i][j].im = 0.0;
            for k in 0..m{
                a[i][j] = a[i][j] + (I[i][k] * a_orig[k][j]); 
            }
        }
    }
    for i in 0..m {
        for j in 0..n {
            if F32Ext::powf(a_orig[i][j].re - a[i][j].re, 2.0) + F32Ext::powf(a_orig[i][j].im - a[i][j].im, 2.0) > eps {
                equal = false;
            }
        }
    }

    if equal == false {
        println!("Pinv failed!");
    }
    else {
        println!("Pinv successful!");
    }

    print_matrix(&a_orig, m, m);
    print_matrix(&a, m, n);
}

fn main() {

    // let mut a: Vec<Vec<Complex32>> = vec![
    //                                     vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
    //                                     vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.4032, im:0.0876}, Complex32{re: 0.1678, im:0.0390}, Complex32{re: 0.5425, im:0.5118}], 
    //                                     vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.3174, im:0.3352}, Complex32{re: 0.9784, im:0.4514}, Complex32{re: -0.4416, im:-1.3188}],
    //                                     vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.4008, im:-0.0504}, Complex32{re: 0.0979, im:-0.2558}, Complex32{re: 0.2983, im:0.7800}]];

    let mut a: Vec<Vec<Complex32>> = vec![
                                        vec![Complex32{re: 0.4032, im:0.0876}, Complex32{re: 0.1678, im:0.0390}, Complex32{re: 0.5425, im:0.5118}], 
                                        vec![Complex32{re: 0.3174, im:0.3352}, Complex32{re: 0.9784, im:0.4514}, Complex32{re: -0.4416, im:-1.3188}],
                                        vec![Complex32{re: 0.4008, im:-0.0504}, Complex32{re: 0.0979, im:-0.2558}, Complex32{re: 0.2983, im:0.7800}],
                                        vec![Complex32{re: 0.1395, im:-0.6213}, Complex32{re: 0.012, im:-0.3587}, Complex32{re: 0.7536, im:0.4729}]];
    let mut inv: Vec<Vec<Complex32>> = vec![
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}]];
    
    check_svd(&mut a, &mut inv, 3, 3) ;
  

}
