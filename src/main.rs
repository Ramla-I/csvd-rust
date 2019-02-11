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
extern crate num_traits;
extern crate backtrace;

use num_complex::Complex32;
use num_traits::float::FloatCore;

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


fn cabs(input: &Complex32) -> f32{
    (input.re.powi(2) + input.im.powi(2)).sqrt()
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
    let tol = 1.5 * (10.0 as f32).powi(-31);

    //10 continue for k in 0..n
    for k in 1..=n {
        k1 = k + 1;

        // Elimination of A(I,K), I = K+1, ..., M.
        let mut z: f32 = 0.0;
        for i in k..=m {
            z = z + (a[i][k].re).powi(2) + (a[i][k].im).powi(2);
        }

        b[k] = 0.0;

        let (mut w, mut q);
        if tol < z {

            z = z.sqrt();
            b[k] = z;
            w = cabs(&a[k][k]);

            if w == 0.0 {
                q = Complex32{ re: 1.0, im: 0.0};
            }
            else {
                q = a[k][k]/w;
            }

            a[k][k] = q * ( z + w );

            if k != (n + p) {
                for j in k1..=(n + p){
                    q = Complex32{ re: 0.0, im: 0.0};
                    
                    for i in k..=m {
                        q = q + a[i][k].conj() * a[i][j];
                    }
                    q = q / ( z * ( z + w ) );

                    for i in k..=m {
                        a[i][j] = a[i][j] - q * a[i][k]
                    }
                }

                // Phase transformation.
                q = -a[k][k].conj() / cabs(&a[k][k]);

                for j in k1..=(n + p) {
                    a[k][j] = q * a[k][j];
                }
            }
        }

        //Elimination of A(K,J), J = K+2, ..., N

        if k == n {
            break;
        }

        z = 0.0;
        for j in k1..=n {
            z = z + a[k][j].re.powi(2) + a[k][j].im.powi(2)
        }
        c[k1] = 0.0;

        if tol < z {
            z = z.sqrt();
            c[k1] = z;
            w = cabs(&a[k][k1]);

            if w == 0.0 {
                q = Complex32{ re: 1.0, im: 0.0};
            }
            else {
                q = a[k][k1] / w;
            }

            a[k][k1] = q * (z + w);

            for i in k1..=m {
                q = Complex32{ re: 0.0, im: 0.0};

                for j in k1..=n {
                    q = q + a[k][j].conj()  * a[i][j];
                }

                q = q / (z * (z + w));

                for j in k1..=n {
                    a[i][j] = a[i][j] - q * a[k][j];
                }
            }
    
            // Phase transformation.
            q = -a[k][k1].conj() / cabs(&a[k][k1]);
            for i in k1..=m {
                a[i][k1] = a[i][k1] * q;
            }
        }
    }

    // Tolerance for negligible elements.
    //140 continue
    let mut eps: f32 = 0.0;
    let eta: f32 = 1.1920929 * (10.0 as f32).powi(-7);
    let mut t: [f32; NBIG] = [0.0; NBIG];

    for k in 1..=n {
       s[k] = b[k];
       t[k] = c[k];
       eps = eps.max(s[k] + t[k]);
    }

    eps = eps * eta;

    // Initialization of U and V.
    if 0 < nu {
        for j in 1..=nu {
            for i in 1..=m {
                u[i][j] = Complex32{re: 0.0, im: 0.0};
            }
            u[j][j] = Complex32{re: 1.0, im: 0.0};
        }
    }

    if 0 < nv {
        for j in 1..=nv {
            for i in 1..=n {
                v[i][j] = Complex32{re: 0.0, im: 0.0};
            }
            v[j][j] = Complex32{re: 1.0, im: 0.0};
        }
    }

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
    let mut r;
    let mut g;

    // QR diagonalization.
    for kk in 1..=n {
        k = n + 1 - kk;

        //Test for split.
        //220 continue
        loop {
            for ll in 1..=k {
                l = k + 1 - ll;
                if t[l].abs() <= eps {
                    //go to 290
                    break;
                }

                if s[l-1].abs() <= eps {
                    //go to 240
                    break;
                }

            }

            //Cancellation of E(L).
            // 240 continue
            if s[l-1].abs() <= eps {
                cs = 0.0;
                sn = 1.0;
                l1 = l - 1;

                for i in l..=k {
                    f = sn * t[i];
                    t[i] = cs * t[i];

                    if f.abs() <= eps {
                        //go to 290
                        break;
                    }

                    h = s[i];
                    w = (f * f + h * h).sqrt();
                    s[i] = w;
                    cs = h / w;
                    sn = - f / w;

                    if 0 < nu {
                        for j in 1..=n {
                            x = u[j][l1].re;
                            y = u[j][i].re;
                            u[j][l1] = Complex32{re: x * cs + y * sn, im: 0.0};
                            u[j][i] = Complex32{re: y * cs - x * sn, im: 0.0};
                        }
                    }

                    if p != 0 {
                        for j in (n + 1)..=(n + p) {
                            q = a[l1][j];
                            r = a[i][j];
                            a[l1][j] = q * cs + r * sn;
                            a[i][j] = r * cs - q * sn;
                        }
                    }
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
            g = (f * f + 1.0).sqrt();
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
                w = (h * h + f * f).sqrt();
                t[i-1] = w;
                cs = f / w;
                sn = h / w;
                f = x * cs + g * sn;
                g = g * cs - x * sn;
                h = y * sn;
                y = y * cs;

                if 0 < nv {
                    for j in 1..=n {
                        x = v[j][i-1].re;
                        w = v[j][i].re;
                        v[j][i-1] = Complex32{re: x * cs + w * sn, im: 0.0};
                        v[j][i] = Complex32{re: w * cs - x * sn, im: 0.0};
                    }
                }

                w = (h * h + f * f).sqrt();
                s[i-1] = w;
                cs = f / w;
                sn = h / w;
                f = cs * g + sn * y;
                x = cs * y - sn * g;

                if 0 < nu {
                    for j in 1..=n {
                        y = u[j][i-1].re;
                        w = u[j][i].re;
                        u[j][i-1] = Complex32{re: y * cs + w * sn, im: 0.0};
                        u[j][i] = Complex32{re: w * cs - y * sn, im: 0.0};
                    }
                }

                if p != 0 {
                    for j in (n + 1)..=(n + p) {
                        q = a[i-1][j];
                        r = a[i][j];
                        a[i-1][j] = q * cs + r * sn;
                        a[i][j] = r * cs - q * sn;
                    }
                }
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
                for j in 1..=n {
                    v[j][k] = -v[j][k];
                }
            }
        }
    }

    let mut j;
    
    // Sort the singular values.
    for k in 1..=n {
        g = -1.0;
        j = k;

        for i in k..=n {
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
               for i in 1..=n {
                    q = v[i][j];
                    v[i][j] = v[i][k];
                    v[i][k] = q;
               }
            }

            // Interchange U(1:N,J) and U(1:N,K).
            if 0 < nu {
                for i in 1..=n {
                    q = u[i][j];
                    u[i][j] = u[i][k];
                    u[i][k] = q;
                }
            }

            // Interchange A(J,N1:NP) and A(K,N1:NP).
            if p != 0 {
                for i in (n + 1)..=(n + p) {
                    q = a[j][i];
                    a[j][i] = a[k][i];
                    a[k][i] = q;
                }
            }
        }
    }

    // Back transformation.
    if 0 < nu {
        for kk in 1..=n {
            k = n + 1 - kk;

            if b[k] != 0.0 {
                q = -a[k][k] / cabs(&a[k][k]);

                for j in 1..=nu {
                    u[k][j] = q * u[k][j];
                }

                for j in 1..=nu {

                    q = Complex32{re: 0.0, im: 0.0};

                    for i in k..=m {
                        q = q + a[i][k].conj() * u[i][j];
                    }

                    q = q / (cabs(&a[k][k]) * b[k]);

                    for i in k..=m {
                        u[i][j] = u[i][j] - q * a[i][k];
                    }

                }

            }

        }

    }

    if 0 < nv {

        if 1 < n {

            for kk in 2..=n {
                k = n + 1 - kk;
                k1 = k + 1;

                if c[k1] != 0.0 { 
                    q = -(a[k][k1].conj()) / cabs(&a[k][k1]);

                    for j in 1..=nv {
                        v[k1][j] = q * v[k1][j];
                    }

                    for j in 1..=nv {
                        q = Complex32{re: 0.0, im: 0.0};

                        for i in k1..=n {
                            q = q + a[k][i] * v[i][j];
                        }
                        q = q / (cabs(&a[k][k1]) * c[k1]);

                        for i in k1..=n {
                            v[i][j] = v[i][j] - q * a[k][i].conj();
                        }
                    }
                }
            }
        }
    }     

    Ok(())   
}

/// Finds the pseudo-inverse of a matrix from the singular value decompositions
/// pinv = U x S x V*
/// Assumes that inv is initialized to 0.0 + 0.0j
fn find_matrix_inv_from_svd(s: &Vec<f32>, u: &Vec<Vec<Complex32>>, v: &Vec<Vec<Complex32>>, m: usize, n: usize, inv: &mut Vec<Vec<Complex32>>) {
    let min = m.min(n);
    for i in 1..=m {
        for j in 1..=n {
            for k in 1..=min {
                inv[i][j] = inv[i][j] + u[i][k] * s[k] * v[j][j].conj();
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
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}]];



    let mut s: Vec<f32> = vec![0.0, 0.0, 0.0, 0.0];

    let _ = csvd(&mut a, 3, 3, 3, 3, 0, 3, 3, &mut s, &mut u, &mut v);
    
    find_matrix_inv_from_svd(&s, &u, &v, 3, 3, &mut inv);

    let mut equal = true;

    for i in 1..=3 {
        for j in 1..=3 {
            if a_orig[i][j] != inv[i][j] {
                equal = false;
            }
        }
    }

    if equal == false {
        println!("SVD pseudo-inverse failed!");
    }

    print_matrix(&a_orig, m, n);
    print_matrix(&inv, m, n);
}

fn main() {

    let mut a: Vec<Vec<Complex32>> = vec![
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.4032, im:0.0876}, Complex32{re: 0.1678, im:0.0390}, Complex32{re: 0.5425, im:0.5118}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.3174, im:0.3352}, Complex32{re: 0.9784, im:0.4514}, Complex32{re: -0.4416, im:-1.3188}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.4008, im:-0.0504}, Complex32{re: 0.0979, im:-0.2558}, Complex32{re: 0.2983, im:0.7800}]];

    let mut inv: Vec<Vec<Complex32>> = vec![
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}], 
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}],
                                        vec![Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}, Complex32{re: 0.0, im:0.0}]];
    
    check_svd(&mut a, &mut inv, 4, 4) ;

    

}
