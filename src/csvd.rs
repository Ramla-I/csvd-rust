use alloc::vec::Vec;

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


use num_complex::Complex32;
use libm::F32Ext;

const NBIG: usize = 150;

fn sqrt(input: f32) -> f32 {
    F32Ext::sqrt(input)
}

fn powf(input: f32, power: f32) -> f32 {
    F32Ext::powf(input, power)
}

fn abs(input: f32) -> f32 {
    F32Ext::abs(input)
}

fn cabs(input: &Complex32) -> f32{
    let a = powf(input.re, 2.0);
    let b = powf(input.im, 2.0); 
    sqrt(a + b)
}

pub fn csvd(a: &mut Vec<Complex32>, mmax: usize, nmax: usize, n: usize, m: usize, p: usize, nu: usize, nv: usize, 
        s: &mut Vec<f32>, u: &mut Vec<Complex32>, v: &mut Vec<Complex32>) 
        -> Result<(), &'static str> {
    
    // debug!("In csvd");

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
    let tol = 1.5 * powf(10.0, -31.0);

    //10 continue for k in 0..n
    for k in 0..n {
        k1 = k + 1;

        // Elimination of A(I,K), I = K+1, ..., M.
        let mut z: f32 = 0.0;
        for i in k..m {
            z = z + powf(a[i*m + k].re, 2.0) + powf(a[i*m + k].im, 2.0);
        }

        b[k] = 0.0;

        let (mut w, mut q);
        if tol < z {

            z = sqrt(z);
            b[k] = z;
            w = cabs(&a[k*m + k]);

            if w == 0.0 {
                q = Complex32{ re: 1.0, im: 0.0};
            }
            else {
                q = a[k*m + k]/w;
            }

            a[k*m + k] = q * ( z + w );

            if k != (n - 1 + p) {
                for j in k1..(n + p){
                    q = Complex32{ re: 0.0, im: 0.0};
                    
                    for i in k..m {
                        q = q + a[i*m + k].conj() * a[i*m + j];
                    }
                    q = q / ( z * ( z + w ) );

                    for i in k..m {
                        a[i*m + j] = a[i*m + j] - q * a[i*m + k];
                    }
                }

                // Phase transformation.
                q = -a[k*m + k].conj() / cabs(&a[k*m + k]);

                for j in k1..(n + p) {
                    a[k*m + j] = q * a[k*m + j];
                }
            }
        }

        //Elimination of A(K,J), J = K+2, ..., N

        if k == (n - 1) {
            break;
        }

        z = 0.0;
        for j in k1..n {
            z = z + powf(a[k*m + j].re, 2.0) + powf(a[k*m + j].im, 2.0);
        }
        c[k1] = 0.0;

        if tol < z {
            z = sqrt(z);
            c[k1] = z;
            w = cabs(&a[k*m + k1]);

            if w == 0.0 {
                q = Complex32{ re: 1.0, im: 0.0};
            }
            else {
                q = a[k*m + k1] / w;
            }

            a[k*m + k1] = q * (z + w);

            for i in k1..m {
                q = Complex32{ re: 0.0, im: 0.0};

                for j in k1..n {
                    q = q + a[k*m + j].conj()  * a[i*m + j];
                }

                q = q / (z * (z + w));

                for j in k1..n {
                    a[i*m + j] = a[i*m + j] - q * a[k*m + j];
                }
            }
    
            // Phase transformation.
            q = -a[k*m + k1].conj() / cabs(&a[k*m + k1]);
            for i in k1..m {
                a[i*m + k1] = a[i* m + k1] * q;
            }
        }
    }

    // Tolerance for negligible elements.
    //140 continue
    let mut eps: f32 = 0.0;
    let eta: f32 = 1.1920929 * powf(10.0, -7.0);
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
                u[i*m + j] = Complex32{re: 0.0, im: 0.0};
            }
            u[j*m + j] = Complex32{re: 1.0, im: 0.0};
        }
    }

    if 0 < nv {
        for j in 0..nv {
            for i in 0..n {
                v[i*n + j] = Complex32{re: 0.0, im: 0.0};
            }
            v[j*n + j] = Complex32{re: 1.0, im: 0.0};
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
                if abs(t[l]) <= eps {
                    //go to 290
                    break;
                }

                if abs(s[l-1]) <= eps {
                    //go to 240
                    break;
                }

            }

            if abs(t[l]) <= eps {
                //go to 290
            }

            //Cancellation of E(L).
            // 240 continue
            else if abs(s[l-1]) <= eps {
                cs = 0.0;
                sn = 1.0;
                l1 = l - 1;

                for i in l..=k {
                    f = sn * t[i];
                    t[i] = cs * t[i];

                    if abs(f) <= eps {
                        //go to 290
                        break;
                    }

                    h = s[i];
                    w = sqrt(f * f + h * h);
                    s[i] = w;
                    cs = h / w;
                    sn = - f / w;

                    if 0 < nu {
                        for j in 0..n {
                            x = u[j*m + l1].re;
                            y = u[j*m + i].re;
                            u[j*m + l1] = Complex32{re: x * cs + y * sn, im: 0.0};
                            u[j*m + i] = Complex32{re: y * cs - x * sn, im: 0.0};
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
            g = sqrt(f * f + 1.0);
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
                w = sqrt(h * h + f * f);
                t[i-1] = w;
                cs = f / w;
                sn = h / w;
                f = x * cs + g * sn;
                g = g * cs - x * sn;
                h = y * sn;
                y = y * cs;

                if 0 < nv {
                    for j in 0..n {
                        x = v[j*n + i-1].re;
                        w = v[j*n + i].re;
                        v[j*n + i-1] = Complex32{re: x * cs + w * sn, im: 0.0};
                        v[j*n + i] = Complex32{re: w * cs - x * sn, im: 0.0};
                    }
                }

                w = sqrt(h * h + f * f);
                s[i-1] = w;
                cs = f / w;
                sn = h / w;
                f = cs * g + sn * y;
                x = cs * y - sn * g;

                if 0 < nu {
                    for j in 0..n {
                        y = u[j*m + i-1].re;
                        w = u[j*m + i].re;
                        u[j*m + i-1] = Complex32{re: y * cs + w * sn, im: 0.0};
                        u[j*m + i] = Complex32{re: w * cs - y * sn, im: 0.0};
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
                    v[j*n + k] = -v[j*n + k];
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
                    q = v[i*n + j];
                    v[i*n + j] = v[i*n + k];
                    v[i*n + k] = q;
               }
            }

            // Interchange U(1:N,J) and U(1:N,K).
            if 0 < nu {
                for i in 0..n {
                    q = u[i*m + j];
                    u[i*m + j] = u[i*m + k];
                    u[i*m + k] = q;
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
                q = -a[k*m + k] / cabs(&a[k*m + k]);

                for j in 0..nu {
                    u[k*m + j] = q * u[k*m + j];
                }

                for j in 0..nu {

                    q = Complex32{re: 0.0, im: 0.0};

                    for i in k..m {
                        q = q + a[i*m + k].conj() * u[i*m + j];
                    }

                    q = q / (cabs(&a[k*m + k]) * b[k]);

                    for i in k..m {
                        u[i*m + j] = u[i*m + j] - q * a[i*m + k];
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
                    q = -(a[k*m + k1].conj()) / cabs(&a[k*m + k1]);

                    for j in 0..nv {
                        v[k1*n + j] = q * v[k1*n + j];
                    }

                    for j in 0..nv {
                        q = Complex32{re: 0.0, im: 0.0};

                        for i in k1..n {
                            q = q + a[k*m + i] * v[i*n + j];
                        }
                        q = q / (cabs(&a[k*m + k1]) * c[k1]);

                        for i in k1..n {
                            v[i*n + j] = v[i*n + j] - q * a[k*m + i].conj();
                        }
                    }
                }
            }
        }
    }     

    Ok(())   
}
