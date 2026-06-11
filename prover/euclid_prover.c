#define _POSIX_C_SOURCE 200809L

/* euclid_prover.c -- rigorous interval-arithmetic prover for unit-edge
 * Euclidean triangulations of S^2.
 *
 * Build:  cc -O3 -std=c11 -Wall -Wextra -pedantic euclid_prover.c -lm -o euclid_prover
 * Usage:  euclid_prover input.obj
 *
 * Do not compile with -ffast-math. The verification layer assumes
 * ordinary IEEE-754 binary64 operations and nextafter padding.
 *
 * Three certificates:
 *
 *   1. EXISTENCE  -- rho_upper < sigma_lower^2 / (16 sqrt(E)).
 *      rho_upper:  interval bound on sum of (||edge||^2 - 1)^2.
 *      sigma_lower:  certified lower bound on smallest singular value
 *      of the Jacobian J of the squared-length map, via Cholesky
 *      witness on JJ^T - s^2 I.
 *
 *   2. EMBEDDING (Ellison sqrt(V))  --  sqrt(V) * motion_upper < collision_lower.
 *      motion_upper:  Ellison radius derived from (sigma_lower, rho_upper, E).
 *      collision_lower:  certified lower bound on min distance between
 *      disjoint simplices (vertex/edge/face pairs), via separating
 *      hyperplane on the GJK-style closest point pair.
 *
 *   3. UNDENTED  --  certified sin(T/2) > 0 at every vertex link, where
 *      T is the total geodesic turning around the spherical link.
 *      Computed without atan2 via half-angle products of (cos t_i,
 *      sin t_i) along the link cycle.
 *
 * ACCEPT iff all three certificates pass; otherwise REJECT.
 *
 * Pulled math from /Users/doyle/Dropbox/neo/euclid_prover/src/euclid_prover.c
 * (the original cruft-laden version). Diagnostic timers, legacy Jacobi
 * sigma path, dead counters (phase1_had_failure, skipped_pairs), and
 * vestigial CLI flags removed.
 */

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

/* unistd.h's int dup(int) collides with our static double dup(double),
 * so forward-declare just the one POSIX function we need. */
extern pid_t getpid(void);

typedef struct { double x, y, z; } Vec3;
typedef struct { int a, b, c; } Face;
typedef struct { int a, b; } Edge;
typedef struct { double lo, hi; } Interval;
typedef struct { int n; int v[3]; } Simplex;
typedef struct { Interval x, y, z; } IVec3;
typedef struct { int degree; int *nbrs; } Cycle;

#define MAT(A, n, i, j) ((A)[(size_t)(i)*(size_t)(n) + (size_t)(j)])

static double ddown(double x) { return nextafter(x, -INFINITY); }
static double dup  (double x) { return nextafter(x,  INFINITY); }

static void die(const char *msg) {
    fprintf(stderr, "error: %s\n", msg);
    exit(1);
}
static void *xcalloc(size_t n, size_t size) {
    void *p = calloc(n, size);
    if (!p) die("out of memory");
    return p;
}
static void *xrealloc(void *p, size_t size) {
    void *q = realloc(p, size);
    if (!q) die("out of memory");
    return q;
}

/* ============================================================
 *   Interval arithmetic with directed rounding via nextafter.
 * ============================================================ */

static Interval ipoint(double x) {
    if (isnan(x)) die("NaN interval point");
    return (Interval){x, x};
}
static Interval izero(void) { return ipoint(0.0); }

static Interval iadd(Interval a, Interval b) {
    Interval r = { ddown(a.lo + b.lo), dup(a.hi + b.hi) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval add");
    return r;
}
static Interval isub(Interval a, Interval b) {
    Interval r = { ddown(a.lo - b.hi), dup(a.hi - b.lo) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval sub");
    return r;
}
static Interval imul(Interval a, Interval b) {
    double v0 = a.lo * b.lo, v1 = a.lo * b.hi;
    double v2 = a.hi * b.lo, v3 = a.hi * b.hi;
    double lo = fmin(fmin(v0, v1), fmin(v2, v3));
    double hi = fmax(fmax(v0, v1), fmax(v2, v3));
    Interval r = { ddown(lo), dup(hi) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval mul");
    return r;
}
static Interval idiv(Interval a, Interval b) {
    if (b.lo <= 0.0 && 0.0 <= b.hi)
        die("interval division by interval containing zero");
    double v0 = a.lo / b.lo, v1 = a.lo / b.hi;
    double v2 = a.hi / b.lo, v3 = a.hi / b.hi;
    double lo = fmin(fmin(v0, v1), fmin(v2, v3));
    double hi = fmax(fmax(v0, v1), fmax(v2, v3));
    Interval r = { ddown(lo), dup(hi) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval div");
    return r;
}
static Interval isquare(Interval a) {
    if (a.lo <= 0.0 && 0.0 <= a.hi) {
        double hi = fmax(a.lo * a.lo, a.hi * a.hi);
        return (Interval){ 0.0, dup(hi) };
    }
    double v0 = a.lo * a.lo, v1 = a.hi * a.hi;
    double lo = ddown(fmin(v0, v1));
    if (lo < 0.0) lo = 0.0;   /* squares are non-negative */
    return (Interval){ lo, dup(fmax(v0, v1)) };
}
static Interval clamp_nonneg(Interval a) {
    if (a.lo < 0.0) a.lo = 0.0;
    if (a.hi < 0.0) a.hi = 0.0;
    return a;
}
static Interval isqrt_interval(Interval a) {
    if (a.lo < 0.0) die("sqrt interval with negative lower endpoint");
    return (Interval){ ddown(sqrt(a.lo)), dup(sqrt(a.hi)) };
}
static bool isqrt_interval_checked(Interval a, Interval *out, char *msg, size_t msglen) {
    if (a.hi < 0.0) {
        snprintf(msg, msglen, "sqrt of strictly-negative interval: [%.17g, %.17g]", a.lo, a.hi);
        return false;
    }
    if (a.lo < 0.0) {
        snprintf(msg, msglen, "sqrt interval has negative lower endpoint: [%.17g, %.17g]", a.lo, a.hi);
        return false;
    }
    out->lo = ddown(sqrt(a.lo));
    out->hi = dup  (sqrt(a.hi));
    return true;
}
static double iabs_upper(Interval a) {
    return dup(fmax(fabs(a.lo), fabs(a.hi)));
}
static bool frob_upper_intervals_checked(const Interval *x, int n, double *out, char *msg, size_t msglen) {
    Interval total = izero();
    for (int i = 0; i < n; i++) {
        Interval a = ipoint(iabs_upper(x[i]));
        total = iadd(total, isquare(a));
    }
    Interval root;
    if (!isqrt_interval_checked(clamp_nonneg(total), &root, msg, msglen)) return false;
    *out = root.hi;
    return true;
}
static bool frob_upper_points_checked(const double *M, int n, double *out, char *msg, size_t msglen) {
    Interval total = izero();
    for (int i = 0; i < n*n; i++) {
        Interval a = ipoint(dup(fabs(M[i])));
        total = iadd(total, isquare(a));
    }
    Interval root;
    if (!isqrt_interval_checked(clamp_nonneg(total), &root, msg, msglen)) return false;
    *out = root.hi;
    return true;
}

/* ============================================================
 *   Vec3 utilities (plain double).
 * ============================================================ */

static double vecdot   (Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static Vec3   vecsub   (Vec3 a, Vec3 b) { return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static Vec3   vecadd   (Vec3 a, Vec3 b) { return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z}; }
static Vec3   vecscale (Vec3 a, double s) { return (Vec3){a.x*s, a.y*s, a.z*s}; }
static double vecnorm  (Vec3 a) { return sqrt(vecdot(a, a)); }

/* ============================================================
 *   OBJ parser + edge enumeration.
 * ============================================================ */

static bool parse_obj(const char *path, Vec3 **vout, int *nvout, Face **fout, int *nfout, char *msg, size_t msglen) {
    *vout = NULL;
    *nvout = 0;
    *fout = NULL;
    *nfout = 0;

    FILE *fp = fopen(path, "r");
    if (!fp) {
        snprintf(msg, msglen, "%s: %s", path, strerror(errno));
        return false;
    }

    int vc = 0, vcap = 64;
    int fc = 0, fcap = 64;
    Vec3 *verts = (Vec3 *)xcalloc((size_t)vcap, sizeof(Vec3));
    Face *faces = (Face *)xcalloc((size_t)fcap, sizeof(Face));
    char line[4096];

    while (fgets(line, sizeof(line), fp)) {
        char *s = line;
        while (isspace((unsigned char)*s)) s++;
        if (*s == '\0' || *s == '#') continue;

        if (s[0] == 'v' && isspace((unsigned char)s[1])) {
            char *q = s + 1;
            double x = strtod(q, &q);
            double y = strtod(q, &q);
            double z = strtod(q, &q);
            if (vc == vcap) {
                vcap *= 2;
                verts = (Vec3 *)xrealloc(verts, (size_t)vcap * sizeof(Vec3));
            }
            verts[vc++] = (Vec3){x, y, z};
        } else if (s[0] == 'f' && isspace((unsigned char)s[1])) {
            int ids[3];
            char *q = s + 1;
            for (int k = 0; k < 3; k++) {
                while (isspace((unsigned char)*q)) q++;
                long idx = strtol(q, &q, 10);
                if (idx < 0) idx = vc + idx + 1;
                ids[k] = (int)idx - 1;
                while (*q && !isspace((unsigned char)*q)) q++;
            }
            if (fc == fcap) {
                fcap *= 2;
                faces = (Face *)xrealloc(faces, (size_t)fcap * sizeof(Face));
            }
            faces[fc++] = (Face){ids[0], ids[1], ids[2]};
        }
    }

    if (ferror(fp)) {
        snprintf(msg, msglen, "%s: read error", path);
        fclose(fp);
        free(verts);
        free(faces);
        return false;
    }
    if (fclose(fp) != 0) {
        snprintf(msg, msglen, "%s: close error: %s", path, strerror(errno));
        free(verts);
        free(faces);
        return false;
    }

    *vout = verts;
    *nvout = vc;
    *fout = faces;
    *nfout = fc;
    return true;
}

static int edge_cmp(const void *pa, const void *pb) {
    const Edge *a = (const Edge *)pa, *b = (const Edge *)pb;
    if (a->a != b->a) return a->a - b->a;
    return a->b - b->b;
}
static Edge make_edge(int a, int b) {
    if (a > b) { int t = a; a = b; b = t; }
    return (Edge){a, b};
}
static Edge *edges_from_faces(const Face *faces, int nf, int *neout) {
    int nraw = 3 * nf;
    Edge *raw = (Edge *)xcalloc((size_t)nraw, sizeof(Edge));
    for (int i = 0; i < nf; i++) {
        raw[3*i+0] = make_edge(faces[i].a, faces[i].b);
        raw[3*i+1] = make_edge(faces[i].b, faces[i].c);
        raw[3*i+2] = make_edge(faces[i].c, faces[i].a);
    }
    qsort(raw, (size_t)nraw, sizeof(Edge), edge_cmp);
    int ne = 0;
    for (int i = 0; i < nraw; i++) {
        if (ne == 0 || raw[i].a != raw[ne-1].a || raw[i].b != raw[ne-1].b)
            raw[ne++] = raw[i];
    }
    *neout = ne;
    return raw;
}

/* ============================================================
 *   Existence certificate: rho_upper and sigma_lower.
 * ============================================================ */

static double *build_jacobian(const Vec3 *v, int nv, const Edge *edges, int ne, int *colsout) {
    int cols = 3 * nv;
    double *J = (double *)xcalloc((size_t)ne * (size_t)cols, sizeof(double));
    for (int r = 0; r < ne; r++) {
        int a = edges[r].a, b = edges[r].b;
        Vec3 d = vecscale(vecsub(v[a], v[b]), 2.0);
        J[(size_t)r*cols + 3*a + 0] =  d.x;
        J[(size_t)r*cols + 3*a + 1] =  d.y;
        J[(size_t)r*cols + 3*a + 2] =  d.z;
        J[(size_t)r*cols + 3*b + 0] = -d.x;
        J[(size_t)r*cols + 3*b + 1] = -d.y;
        J[(size_t)r*cols + 3*b + 2] = -d.z;
    }
    *colsout = cols;
    return J;
}

static Interval interval_length_squared(const Vec3 *v, int a, int b) {
    Interval total = izero();
    Interval dx = isub(ipoint(v[a].x), ipoint(v[b].x)); total = iadd(total, isquare(dx));
    Interval dy = isub(ipoint(v[a].y), ipoint(v[b].y)); total = iadd(total, isquare(dy));
    Interval dz = isub(ipoint(v[a].z), ipoint(v[b].z)); total = iadd(total, isquare(dz));
    return total;
}

static double certify_rho_upper(const Vec3 *v, const Edge *edges, int ne) {
    Interval total = izero();
    Interval one = ipoint(1.0);
    for (int i = 0; i < ne; i++) {
        Interval l2  = interval_length_squared(v, edges[i].a, edges[i].b);
        Interval err = isub(l2, one);
        total = iadd(total, isquare(err));
    }
    return isqrt_interval(clamp_nonneg(total)).hi;
}

static double *compute_jjt(const double *J, int m, int cols) {
    double *G = (double *)xcalloc((size_t)m * (size_t)m, sizeof(double));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j <= i; j++) {
            double s = 0.0;
            const double *ri = J + (size_t)i * cols;
            const double *rj = J + (size_t)j * cols;
            for (int k = 0; k < cols; k++) s += ri[k] * rj[k];
            MAT(G, m, i, j) = MAT(G, m, j, i) = s;
        }
    }
    return G;
}

static bool cholesky_lower(const double *A, int n, double *C) {
    memset(C, 0, (size_t)n * (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = MAT(A, n, i, j);
            for (int k = 0; k < j; k++) sum -= MAT(C, n, i, k) * MAT(C, n, j, k);
            if (i == j) {
                if (!(sum > 0.0) || !isfinite(sum)) return false;
                MAT(C, n, i, j) = sqrt(sum);
            } else {
                MAT(C, n, i, j) = sum / MAT(C, n, j, j);
            }
        }
    }
    return true;
}

static bool lower_triangular_inverse(const double *L, int n, double *B) {
    memset(B, 0, (size_t)n * (size_t)n * sizeof(double));
    for (int col = 0; col < n; col++) {
        for (int i = 0; i < n; i++) {
            double rhs = (i == col) ? 1.0 : 0.0;
            for (int k = 0; k < i; k++) rhs -= MAT(L, n, i, k) * MAT(B, n, k, col);
            double diag = MAT(L, n, i, i);
            if (!(fabs(diag) > 0.0) || !isfinite(diag)) return false;
            MAT(B, n, i, col) = rhs / diag;
        }
    }
    return true;
}

static Interval interval_bc_entry(const double *B, const double *C, int n, int i, int j) {
    Interval total = izero();
    for (int k = 0; k < n; k++)
        total = iadd(total, imul(ipoint(MAT(B, n, i, k)), ipoint(MAT(C, n, k, j))));
    return total;
}
static Interval interval_cct_entry(const double *C, int n, int i, int j) {
    Interval total = izero();
    for (int k = 0; k < n; k++)
        total = iadd(total, imul(ipoint(MAT(C, n, i, k)), ipoint(MAT(C, n, j, k))));
    return total;
}
static Interval interval_a_entry(const double *J, int cols, double s, int i, int j) {
    Interval total = izero();
    for (int k = 0; k < cols; k++) {
        double x = J[(size_t)i*cols + k], y = J[(size_t)j*cols + k];
        if (x != 0.0 && y != 0.0) total = iadd(total, imul(ipoint(x), ipoint(y)));
    }
    if (i == j) total = isub(total, isquare(ipoint(s)));
    return total;
}

/* Cholesky witness on JJ^T - s^2 I:
 *   - compute C = chol(A), B = C^{-1} (lower triangular)
 *   - interval-verify ||I - BC||_F < 1  ⇒  sigma_min(C) ≥ (1 - delta) / ||B||_F > 0
 *   - interval-verify ||A - CC^T||_F < sigma_min(C)^2  ⇒ A is positive-definite
 *     hence A = JJ^T - s^2 I > 0, hence s < sigma_min(J).
 */
static bool certify_factor_witness(const double *J, int m, int cols, double s, const double *C, char *msg, size_t msglen) {
    double *B = (double *)xcalloc((size_t)m * (size_t)m, sizeof(double));
    if (!lower_triangular_inverse(C, m, B)) {
        snprintf(msg, msglen, "could not invert Cholesky witness");
        free(B); return false;
    }
    Interval *errs = (Interval *)xcalloc((size_t)m * (size_t)m, sizeof(Interval));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
            Interval bc = interval_bc_entry(B, C, m, i, j);
            Interval target = ipoint((i == j) ? 1.0 : 0.0);
            errs[(size_t)i*m + j] = isub(target, bc);
        }
    double delta;
    if (!frob_upper_intervals_checked(errs, m*m, &delta, msg, msglen)) {
        free(errs); free(B); return false;
    }
    free(errs);
    if (!(delta < 1.0)) {
        snprintf(msg, msglen, "inverse witness failed: ||I-BC||_F <= %.17g", delta);
        free(B); return false;
    }
    double bnorm;
    if (!frob_upper_points_checked(B, m, &bnorm, msg, msglen)) {
        free(B); return false;
    }
    double factor_sigma_lower = ddown((1.0 - delta) / bnorm);
    if (!(factor_sigma_lower > 0.0)) {
        snprintf(msg, msglen, "factor singular-value bound is not positive");
        free(B); return false;
    }
    Interval *res = (Interval *)xcalloc((size_t)m * (size_t)m, sizeof(Interval));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
            Interval aij = interval_a_entry(J, cols, s, i, j);
            Interval cij = interval_cct_entry(C, m, i, j);
            res[(size_t)i*m + j] = isub(aij, cij);
        }
    double residual;
    if (!frob_upper_intervals_checked(res, m*m, &residual, msg, msglen)) {
        free(res); free(B); return false;
    }
    free(res);
    double margin = ddown(factor_sigma_lower * factor_sigma_lower - residual);
    if (margin > 0.0) {
        snprintf(msg, msglen,
                 "witness verified: ||I-BC||_F <= %.17g, sigma_min(C) >= %.17g, "
                 "||A-CCT||_F <= %.17g, margin >= %.17g",
                 delta, factor_sigma_lower, residual, margin);
        free(B); return true;
    }
    snprintf(msg, msglen,
             "factor witness margin failed: sigma_min(C)^2 lower %.17g, residual upper %.17g",
             ddown(factor_sigma_lower * factor_sigma_lower), residual);
    free(B); return false;
}

/* Two-phase certified sigma_lower search.
 *   Phase 1: exponential descent from s=1.0, halving on failure, until
 *            the first s for which Cholesky AND the witness both pass.
 *   Phase 2: geometric-mean bisection in [s_first, 2*s_first] for up to
 *            max_refine steps; tightens s_lower without changing the
 *            soundness of phase 1's certificate. */
static bool certify_sigma_lower_crude(const double *J, int m, int cols, int max_descend, int max_refine, double *sout, char *msg, size_t msglen) {
    double *G = compute_jjt(J, m, cols);
    double *A = (double *)xcalloc((size_t)m * (size_t)m, sizeof(double));
    double *C = (double *)xcalloc((size_t)m * (size_t)m, sizeof(double));
    char last[2048] = {0};
    char success_msg[2048] = {0};

    double s = 1.0;
    double s_first = 0.0;
    bool found = false;
    for (int i = 0; i < max_descend; i++) {
        memcpy(A, G, (size_t)m * (size_t)m * sizeof(double));
        for (int j = 0; j < m; j++) MAT(A, m, j, j) -= s * s;
        if (cholesky_lower(A, m, C)) {
            if (certify_factor_witness(J, m, cols, s, C, msg, msglen)) {
                s_first = s;
                snprintf(success_msg, sizeof(success_msg), "%s", msg);
                found = true;
                break;
            }
            snprintf(last, sizeof(last), "%s", msg);
        } else {
            snprintf(last, sizeof(last), "ordinary Cholesky failed at s=%.17g", s);
        }
        s = ddown(s * 0.5);
    }
    if (!found) {
        snprintf(msg, msglen, "could not certify sigma lower (phase 1); last failure: %s", last);
        free(G); free(A); free(C); return false;
    }

    double s_lo = s_first;
    double s_hi = ddown(2.0 * s_first);
    for (int i = 0; i < max_refine; i++) {
        double s_mid = ddown(sqrt(s_lo * s_hi));
        if (!(s_mid > s_lo) || !(s_mid < s_hi)) break;
        memcpy(A, G, (size_t)m * (size_t)m * sizeof(double));
        for (int j = 0; j < m; j++) MAT(A, m, j, j) -= s_mid * s_mid;
        bool ok = false;
        if (cholesky_lower(A, m, C))
            ok = certify_factor_witness(J, m, cols, s_mid, C, msg, msglen);
        if (ok) {
            s_lo = s_mid;
            snprintf(success_msg, sizeof(success_msg), "%s", msg);
        } else {
            s_hi = s_mid;
        }
    }
    snprintf(msg, msglen, "%s", success_msg);
    *sout = s_lo;
    free(G); free(A); free(C);
    return true;
}

/* ============================================================
 *   Embedding certificate: collision_lower.
 * ============================================================ */

static bool solve_linear(double *A, double *b, int n, double *x) {
    int N = n;
    double *M = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    double *rhs = (double *)xcalloc((size_t)N, sizeof(double));
    memcpy(M,   A, (size_t)N * (size_t)N * sizeof(double));
    memcpy(rhs, b, (size_t)N * sizeof(double));
    for (int k = 0; k < N; k++) {
        int piv = k;
        double best = fabs(MAT(M, N, k, k));
        for (int i = k+1; i < N; i++) {
            double val = fabs(MAT(M, N, i, k));
            if (val > best) { best = val; piv = i; }
        }
        if (!(best > 1e-14)) { free(M); free(rhs); return false; }
        if (piv != k) {
            for (int j = k; j < N; j++) {
                double tmp = MAT(M, N, k, j);
                MAT(M, N, k, j) = MAT(M, N, piv, j);
                MAT(M, N, piv, j) = tmp;
            }
            double tr = rhs[k]; rhs[k] = rhs[piv]; rhs[piv] = tr;
        }
        for (int i = k+1; i < N; i++) {
            double f = MAT(M, N, i, k) / MAT(M, N, k, k);
            MAT(M, N, i, k) = 0.0;
            for (int j = k+1; j < N; j++) MAT(M, N, i, j) -= f * MAT(M, N, k, j);
            rhs[i] -= f * rhs[k];
        }
    }
    for (int i = N-1; i >= 0; i--) {
        double s = rhs[i];
        for (int j = i+1; j < N; j++) s -= MAT(M, N, i, j) * x[j];
        x[i] = s / MAT(M, N, i, i);
    }
    free(M); free(rhs); return true;
}

static Vec3 affine_combo(const Vec3 *p, const int *idx, const double *w, int n) {
    Vec3 out = {0, 0, 0};
    for (int i = 0; i < n; i++) out = vecadd(out, vecscale(p[idx[i]], w[i]));
    return out;
}

static void closest_points_active_sets(const Vec3 *verts, const Simplex *A, const Simplex *B, Vec3 *bestp, Vec3 *bestq) {
    double best = INFINITY;
    *bestp = verts[A->v[0]];
    *bestq = verts[B->v[0]];
    int maxmaskA = 1 << A->n, maxmaskB = 1 << B->n;
    for (int ma = 1; ma < maxmaskA; ma++) {
        int ia[3], ca = 0;
        for (int i = 0; i < A->n; i++) if ((ma >> i) & 1) ia[ca++] = A->v[i];
        for (int mb = 1; mb < maxmaskB; mb++) {
            int ib[3], cb = 0;
            for (int i = 0; i < B->n; i++) if ((mb >> i) & 1) ib[cb++] = B->v[i];
            int vars = ca + cb;
            int N = vars + 2;
            double M[64]; double rhs[8]; double sol[8];
            memset(M, 0, sizeof(M));
            memset(rhs, 0, sizeof(rhs));
            memset(sol, 0, sizeof(sol));
            for (int i = 0; i < ca; i++)
                for (int j = 0; j < ca; j++)
                    MAT(M, N, i, j) =  2.0 * vecdot(verts[ia[i]], verts[ia[j]]);
            for (int i = 0; i < ca; i++)
                for (int j = 0; j < cb; j++)
                    MAT(M, N, i, ca+j) = -2.0 * vecdot(verts[ia[i]], verts[ib[j]]);
            for (int i = 0; i < cb; i++)
                for (int j = 0; j < ca; j++)
                    MAT(M, N, ca+i, j) = -2.0 * vecdot(verts[ib[i]], verts[ia[j]]);
            for (int i = 0; i < cb; i++)
                for (int j = 0; j < cb; j++)
                    MAT(M, N, ca+i, ca+j) = 2.0 * vecdot(verts[ib[i]], verts[ib[j]]);
            for (int i = 0; i < ca; i++) {
                MAT(M, N, i, vars)   = 1.0;
                MAT(M, N, vars, i)   = 1.0;
            }
            for (int i = 0; i < cb; i++) {
                MAT(M, N, ca+i, vars+1) = 1.0;
                MAT(M, N, vars+1, ca+i) = 1.0;
            }
            rhs[vars] = 1.0; rhs[vars+1] = 1.0;
            if (!solve_linear(M, rhs, N, sol)) continue;
            bool ok = true;
            for (int i = 0; i < vars; i++) if (sol[i] < -1e-9) ok = false;
            if (!ok) continue;
            double la[3], lb[3], suma = 0.0, sumb = 0.0;
            for (int i = 0; i < ca; i++) { la[i] = fmax(0.0, sol[i]);      suma += la[i]; }
            for (int i = 0; i < cb; i++) { lb[i] = fmax(0.0, sol[ca + i]); sumb += lb[i]; }
            if (suma == 0.0 || sumb == 0.0) continue;
            for (int i = 0; i < ca; i++) la[i] /= suma;
            for (int i = 0; i < cb; i++) lb[i] /= sumb;
            Vec3 p = affine_combo(verts, ia, la, ca);
            Vec3 q = affine_combo(verts, ib, lb, cb);
            double d = vecnorm(vecsub(p, q));
            if (d < best) { best = d; *bestp = p; *bestq = q; }
        }
    }
}

static Interval interval_dot_vec(Vec3 n, Vec3 v) {
    Interval total = izero();
    total = iadd(total, imul(ipoint(n.x), ipoint(v.x)));
    total = iadd(total, imul(ipoint(n.y), ipoint(v.y)));
    total = iadd(total, imul(ipoint(n.z), ipoint(v.z)));
    return total;
}
static double norm_upper_vec(Vec3 n) {
    Interval total = izero();
    total = iadd(total, isquare(ipoint(n.x)));
    total = iadd(total, isquare(ipoint(n.y)));
    total = iadd(total, isquare(ipoint(n.z)));
    return isqrt_interval(clamp_nonneg(total)).hi;
}

static double separation_lower_for_normal(const Vec3 *verts, const Simplex *A, const Simplex *B, Vec3 n) {
    double nup = norm_upper_vec(n);
    if (!(nup > 0.0) || !isfinite(nup)) return -INFINITY;
    double minA = INFINITY, maxA = -INFINITY, minB = INFINITY, maxB = -INFINITY;
    for (int i = 0; i < A->n; i++) {
        Interval d = interval_dot_vec(n, verts[A->v[i]]);
        if (d.lo < minA) minA = d.lo;
        if (d.hi > maxA) maxA = d.hi;
    }
    for (int i = 0; i < B->n; i++) {
        Interval d = interval_dot_vec(n, verts[B->v[i]]);
        if (d.lo < minB) minB = d.lo;
        if (d.hi > maxB) maxB = d.hi;
    }
    double gap1 = ddown(minA - maxB);
    double gap2 = ddown(minB - maxA);
    double gap = fmax(gap1, gap2);
    if (gap <= 0.0) return -INFINITY;
    return ddown(gap / nup);
}

static Simplex *build_simplices(int nv, const Edge *edges, int ne, const Face *faces, int nf, int *nsout) {
    int ns = nv + ne + nf;
    Simplex *s = (Simplex *)xcalloc((size_t)ns, sizeof(Simplex));
    int k = 0;
    for (int i = 0; i < nv; i++) { s[k].n = 1; s[k].v[0] = i; k++; }
    for (int i = 0; i < ne; i++) {
        s[k].n = 2; s[k].v[0] = edges[i].a; s[k].v[1] = edges[i].b; k++;
    }
    for (int i = 0; i < nf; i++) {
        s[k].n = 3;
        s[k].v[0] = faces[i].a; s[k].v[1] = faces[i].b; s[k].v[2] = faces[i].c;
        k++;
    }
    *nsout = ns;
    return s;
}

static bool simplices_disjoint(const Simplex *a, const Simplex *b) {
    for (int i = 0; i < a->n; i++)
        for (int j = 0; j < b->n; j++)
            if (a->v[i] == b->v[j]) return false;
    return true;
}

static Vec3 simplex_centroid(const Vec3 *verts, const Simplex *s) {
    Vec3 c = {0, 0, 0};
    for (int i = 0; i < s->n; i++) c = vecadd(c, verts[s->v[i]]);
    return vecscale(c, 1.0 / (double)s->n);
}

static bool certify_collision_lower(const Vec3 *verts, int nv, const Edge *edges, int ne, const Face *faces, int nf, double *Dout, Simplex *pa, Simplex *pb, int *checkedout) {
    int ns = 0;
    Simplex *simp = build_simplices(nv, edges, ne, faces, nf, &ns);
    double best_lower = INFINITY;
    int checked = 0;
    for (int i = 0; i < ns; i++) {
        for (int j = i+1; j < ns; j++) {
            if (!simplices_disjoint(&simp[i], &simp[j])) continue;
            checked++;
            Vec3 p, q;
            closest_points_active_sets(verts, &simp[i], &simp[j], &p, &q);
            Vec3 candidates[16];
            int nc = 0;
            Vec3 diff = vecsub(p, q);
            if (vecnorm(diff) > 0.0) candidates[nc++] = diff;
            Vec3 cdiff = vecsub(simplex_centroid(verts, &simp[i]),
                                simplex_centroid(verts, &simp[j]));
            if (vecnorm(cdiff) > 0.0 && nc < 16) candidates[nc++] = cdiff;
            for (int a = 0; a < simp[i].n && nc < 16; a++) {
                for (int b = 0; b < simp[j].n && nc < 16; b++) {
                    Vec3 vd = vecsub(verts[simp[i].v[a]], verts[simp[j].v[b]]);
                    if (vecnorm(vd) > 0.0) candidates[nc++] = vd;
                }
            }
            double pair_lower = -INFINITY;
            for (int c = 0; c < nc; c++) {
                double lb = separation_lower_for_normal(verts, &simp[i], &simp[j], candidates[c]);
                if (lb > pair_lower) pair_lower = lb;
            }
            if (!(pair_lower > 0.0)) {
                fprintf(stderr, "failed to certify separation for simplex pair %d,%d\n", i, j);
                free(simp); return false;
            }
            if (pair_lower < best_lower) {
                best_lower = pair_lower;
                *pa = simp[i]; *pb = simp[j];
            }
        }
    }
    *Dout = best_lower;
    *checkedout = checked;
    free(simp);
    return true;
}

/* ============================================================
 *   Undented certificate: sin(T/2) > 0 at every vertex link.
 * ============================================================ */

static IVec3 ivec_point_box(Vec3 p, double r) {
    return (IVec3){
        { ddown(p.x - r), dup(p.x + r) },
        { ddown(p.y - r), dup(p.y + r) },
        { ddown(p.z - r), dup(p.z + r) },
    };
}
static IVec3 ivec_sub(IVec3 a, IVec3 b) {
    return (IVec3){ isub(a.x, b.x), isub(a.y, b.y), isub(a.z, b.z) };
}
static Interval idot_vec(IVec3 a, IVec3 b) {
    Interval total = izero();
    total = iadd(total, imul(a.x, b.x));
    total = iadd(total, imul(a.y, b.y));
    total = iadd(total, imul(a.z, b.z));
    return total;
}
static IVec3 icross_vec(IVec3 a, IVec3 b) {
    return (IVec3){
        isub(imul(a.y, b.z), imul(a.z, b.y)),
        isub(imul(a.z, b.x), imul(a.x, b.z)),
        isub(imul(a.x, b.y), imul(a.y, b.x)),
    };
}
static Interval inorm_vec(IVec3 a) {
    Interval total = izero();
    total = iadd(total, isquare(a.x));
    total = iadd(total, isquare(a.y));
    total = iadd(total, isquare(a.z));
    return clamp_nonneg(isqrt_interval(clamp_nonneg(total)));
}
static bool inormalize_vec(IVec3 a, IVec3 *out, const char *what, char *msg, size_t msglen) {
    Interval n = inorm_vec(a);
    if (!(n.lo > 0.0)) {
        snprintf(msg, msglen,
                 "cannot normalize interval vector for %s: norm interval [%.17g, %.17g]",
                 what, n.lo, n.hi);
        return false;
    }
    out->x = idiv(a.x, n);
    out->y = idiv(a.y, n);
    out->z = idiv(a.z, n);
    return true;
}

static bool turning_sincos_interval(IVec3 uprev, IVec3 ucur, IVec3 unext, Interval *sout, Interval *cout, char *msg, size_t msglen) {
    IVec3 nprev, nnext, tin, tout;
    if (!inormalize_vec(icross_vec(uprev, ucur), &nprev, "previous great-circle normal", msg, msglen)) return false;
    if (!inormalize_vec(icross_vec(ucur,  unext), &nnext, "next great-circle normal",     msg, msglen)) return false;
    if (!inormalize_vec(icross_vec(nprev, ucur),  &tin,   "incoming tangent",              msg, msglen)) return false;
    if (!inormalize_vec(icross_vec(nnext, ucur),  &tout,  "outgoing tangent",              msg, msglen)) return false;
    *sout = idot_vec(ucur, icross_vec(tin, tout));
    *cout = idot_vec(tin, tout);
    return true;
}

static void complex_mul_i(Interval ar, Interval ai, Interval br, Interval bi, Interval *rr, Interval *ri) {
    *rr = isub(imul(ar, br), imul(ai, bi));
    *ri = iadd(imul(ar, bi), imul(ai, br));
}

static void free_cycles(Cycle *cycles, int nv);

/* Build the cyclically-ordered link at each vertex.  The OBJ is assumed
 * to be a properly oriented sphere triangulation; euclid_check validates
 * that separately when needed. */
static Cycle *build_vertex_cycles_ccw_outside(int nv, const Face *faces, int nf) {
    int *nxt = (int *)xcalloc((size_t)nv * (size_t)nv, sizeof(int));
    for (int i = 0; i < nv*nv; i++) nxt[i] = -1;

    for (int f = 0; f < nf; f++) {
        int tri[3] = {faces[f].a, faces[f].b, faces[f].c};
        for (int k = 0; k < 3; k++) {
            int v = tri[k];
            int x = tri[(k+2) % 3];
            int y = tri[(k+1) % 3];
            nxt[(size_t)v*(size_t)nv + (size_t)x] = y;
        }
    }

    Cycle *cycles = (Cycle *)xcalloc((size_t)nv, sizeof(Cycle));
    for (int v = 0; v < nv; v++) {
        int degree = 0;
        int start = -1;
        for (int j = 0; j < nv; j++) {
            if (nxt[(size_t)v*(size_t)nv + (size_t)j] >= 0) {
                if (start < 0) start = j;
                degree++;
            }
        }

        cycles[v].degree = degree;
        cycles[v].nbrs = (int *)xcalloc((size_t)degree, sizeof(int));
        int cur = start;
        for (int k = 0; k < degree; k++) {
            cycles[v].nbrs[k] = cur;
            cur = nxt[(size_t)v*(size_t)nv + (size_t)cur];
        }
    }

    free(nxt);
    return cycles;
}

static void free_cycles(Cycle *cycles, int nv) {
    if (!cycles) return;
    for (int i = 0; i < nv; i++) free(cycles[i].nbrs);
    free(cycles);
}

/* For one vertex link: each local geodesic turn t_i gives (cos t_i, sin t_i)
 * from dot/cross products on the unit sphere. cos(t_i/2) = sqrt((1+cos)/2)
 * > 0 by the principal-branch choice; sin(t_i/2) = sin / (2 cos(t_i/2)).
 * Composing the half-angle factors gives (cos(T/2), sin(T/2)) where T is
 * the total link turning. If the spherical link is itself embedded then
 * −2π < T < 2π and T > 0 ⟺ sin(T/2) > 0. */
static bool total_turning_sin_half_interval_for_cycle(const IVec3 *boxes, int v, const Cycle *cycle, Interval *out, char *msg, size_t msglen) {
    int d = cycle->degree;
    IVec3 *dirs = (IVec3 *)xcalloc((size_t)d, sizeof(IVec3));
    for (int i = 0; i < d; i++) {
        int a = cycle->nbrs[i];
        if (!inormalize_vec(ivec_sub(boxes[a], boxes[v]), &dirs[i], "vertex direction", msg, msglen)) {
            free(dirs); return false;
        }
    }
    Interval re  = ipoint(1.0);
    Interval im  = ipoint(0.0);
    Interval one = ipoint(1.0);
    Interval two = ipoint(2.0);
    for (int i = 0; i < d; i++) {
        Interval s, c;
        if (!turning_sincos_interval(dirs[(i+d-1)%d], dirs[i], dirs[(i+1)%d], &s, &c, msg, msglen)) {
            free(dirs); return false;
        }
        Interval half_cos_sq = idiv(iadd(one, c), two);
        if (!(half_cos_sq.lo > 0.0)) {
            snprintf(msg, msglen,
                     "cannot certify positive half-angle cosine at link corner %d: (1+c)/2 interval=[%.17g, %.17g]",
                     i+1, half_cos_sq.lo, half_cos_sq.hi);
            free(dirs); return false;
        }
        Interval ch = isqrt_interval(half_cos_sq);
        Interval sh = idiv(s, imul(two, ch));
        Interval nr, ni;
        complex_mul_i(re, im, ch, sh, &nr, &ni);
        re = nr; im = ni;
    }
    *out = im;
    free(dirs);
    return true;
}

static bool certify_undented(const Vec3 *verts, int nv, const Face *faces, int nf, double motion_radius, double *min_lower_out, int *min_vertex_out, char *msg, size_t msglen) {
    Cycle *cycles = build_vertex_cycles_ccw_outside(nv, faces, nf);
    IVec3 *boxes  = (IVec3 *)xcalloc((size_t)nv, sizeof(IVec3));
    for (int i = 0; i < nv; i++) boxes[i] = ivec_point_box(verts[i], motion_radius);
    double min_lower = INFINITY;
    int    min_vertex = -1;
    bool   ok = true;
    char   first_failure[1024] = {0};
    for (int v = 0; v < nv; v++) {
        int d = cycles[v].degree;
        Interval sin_half = { -INFINITY, INFINITY };
        char local_msg[512] = {0};
        bool vertex_ok = total_turning_sin_half_interval_for_cycle(boxes, v, &cycles[v], &sin_half, local_msg, sizeof(local_msg));
        double lower = vertex_ok ? sin_half.lo : -INFINITY;
        double upper = vertex_ok ? sin_half.hi : INFINITY;
        bool pass = vertex_ok && lower > 0.0;
        if (lower < min_lower) { min_lower = lower; min_vertex = v + 1; }
        if (!pass && ok) {
            ok = false;
            if (vertex_ok)
                snprintf(first_failure, sizeof(first_failure),
                         "vertex %d: sin(T/2) lower bound not positive; degree=%d, sin_half_interval=[%.17g, %.17g]",
                         v+1, d, lower, upper);
            else
                snprintf(first_failure, sizeof(first_failure),
                         "vertex %d: sin(T/2) interval evaluation failed: %s",
                         v+1, local_msg);
        }
    }
    *min_lower_out  = min_lower;
    *min_vertex_out = min_vertex;
    if (ok) snprintf(msg, msglen, "all %d vertex sin(T/2) lower bounds are positive under ccw-outside face convention", nv);
    else    snprintf(msg, msglen, "%s", first_failure);
    free(boxes);
    free_cycles(cycles, nv);
    return ok;
}

/* ============================================================
 *   Motion bound: Ellison-style upper bound on vertex motion given
 *   sigma_lower, rho_upper, edge count E.
 * ============================================================ */

static double upper_motion_bound(double sigma_lower, double rho_upper, int edge_count) {
    double sqrtE = dup(sqrt((double)edge_count));
    double disc  = ddown(sigma_lower * sigma_lower - dup(16.0 * rho_upper * sqrtE));
    if (!(disc > 0.0)) return INFINITY;
    double root        = ddown(sqrt(disc));
    double numerator   = dup(sigma_lower - root);
    double denominator = ddown(8.0 * sqrtE);
    return dup(numerator / denominator);
}

static void print_simplex(FILE *out, Simplex s) {
    fprintf(out, "(");
    for (int i = 0; i < s.n; i++) {
        if (i) fprintf(out, ",");
        fprintf(out, "%d", s.v[i] + 1);
    }
    fprintf(out, ")");
}

/* ============================================================
 *   One-item driver and batch helpers.
 * ============================================================ */

static void set_errmsg(char *errmsg, size_t errmsgsize, const char *msg)
{
    if (errmsgsize > 0) snprintf(errmsg, errmsgsize, "%s", msg);
}

static int prove_one_obj(const char *path, FILE *reportout,
                         char *errmsg, size_t errmsgsize)
{
    Vec3 *verts = NULL;
    Face *faces = NULL;
    Edge *edges = NULL;
    double *J = NULL;
    int nv = 0, nf = 0, ne = 0, cols = 0;
    int ret = 1;

    if (errmsgsize > 0) errmsg[0] = '\0';

    if (!parse_obj(path, &verts, &nv, &faces, &nf, errmsg, errmsgsize)) {
        return 1;
    }
    edges = edges_from_faces(faces, nf, &ne);
    J = build_jacobian(verts, nv, edges, ne, &cols);

    /* --- existence certificate --- */
    double rho_upper = certify_rho_upper(verts, edges, ne);
    double sigma_lower = 0.0;
    char   sigma_msg[1024];
    bool   sigma_ok = certify_sigma_lower_crude(J, ne, cols, /*max_descend*/80, /*max_refine*/4,
                                                &sigma_lower, sigma_msg, sizeof(sigma_msg));

    if (!sigma_ok) {
        fprintf(reportout, "vertices: %d\nfaces: %d\nedges: %d\n", nv, nf, ne);
        fprintf(reportout, "rho_upper: %.17g\n", rho_upper);
        fprintf(reportout, "sigma_certificate: FAILED: %s\n", sigma_msg);
        fprintf(reportout, "final: REJECT\n");
        set_errmsg(errmsg, errmsgsize, "sigma certificate failed");
        ret = 2;
        goto done;
    }

    /* --- embedding certificate (Ellison sqrt(V) motion test) --- */
    Simplex ca = {0, {0,0,0}}, cb = {0, {0,0,0}};
    double  collision_lower = 0.0;
    int     checked = 0;
    bool    collision_ok = certify_collision_lower(verts, nv, edges, ne, faces, nf,
                                                   &collision_lower, &ca, &cb, &checked);
    if (!collision_ok) {
        fprintf(reportout, "vertices: %d\nfaces: %d\nedges: %d\n", nv, nf, ne);
        fprintf(reportout, "rho_upper: %.17g\n", rho_upper);
        fprintf(reportout, "sigma_lower_certified: %.17g\n", sigma_lower);
        fprintf(reportout, "collision_certificate: FAILED\n");
        fprintf(reportout, "final: REJECT\n");
        set_errmsg(errmsg, errmsgsize, "collision certificate failed");
        ret = 2;
        goto done;
    }

    double sqrtE         = dup(sqrt((double)ne));
    double existence_rhs = ddown((sigma_lower * sigma_lower) / dup(16.0 * sqrtE));
    bool   existence_ok  = rho_upper < existence_rhs;
    double motion_upper  = upper_motion_bound(sigma_lower, rho_upper, ne);
    double ellison_value = dup(dup(sqrt((double)nv)) * motion_upper);
    bool   embedding_ok  = ellison_value < collision_lower;

    /* --- undented certificate --- */
    double turning_min_lower = -INFINITY;
    int    turning_min_vertex = -1;
    char   turning_msg[1024];
    bool   undented_ok = certify_undented(verts, nv, faces, nf, motion_upper,
                                          &turning_min_lower, &turning_min_vertex,
                                          turning_msg, sizeof(turning_msg));

    bool accepted = existence_ok && embedding_ok && undented_ok;

    fprintf(reportout, "vertices: %d\nfaces: %d\nedges: %d\n", nv, nf, ne);
    fprintf(reportout, "rho_upper: %.17g\n", rho_upper);
    fprintf(reportout, "sigma_lower_certified: %.17g\n", sigma_lower);
    fprintf(reportout, "sigma_certificate: %s\n", sigma_msg);
    fprintf(reportout, "existence_test: %s  rhs=sigma^2/(16 sqrt E)=%.17g\n",
            existence_ok ? "PASS" : "FAIL", existence_rhs);
    fprintf(reportout, "motion_upper: %.17g\n", motion_upper);
    fprintf(reportout, "collision_lower: %.17g\n", collision_lower);
    fprintf(reportout, "closest_certified_pair: "); print_simplex(reportout, ca); fprintf(reportout, " "); print_simplex(reportout, cb); fprintf(reportout, "\n");
    fprintf(reportout, "collision_certificate: checked %d disjoint simplex pairs\n", checked);
    fprintf(reportout, "embedding_test: %s  sqrt(V)*motion_upper=%.17g < collision_lower=%.17g\n",
            embedding_ok ? "PASS" : "FAIL", ellison_value, collision_lower);
    fprintf(reportout, "undented_test: %s  min_sin_half_at_v%d_lower=%.17g  %s\n",
            undented_ok ? "PASS" : "FAIL", turning_min_vertex, turning_min_lower, turning_msg);
    fprintf(reportout, "final: %s\n", accepted ? "ACCEPT" : "REJECT");

    if (accepted) {
        set_errmsg(errmsg, errmsgsize, "ACCEPT");
        ret = 0;
    } else {
        set_errmsg(errmsg, errmsgsize, "final REJECT");
        ret = 2;
    }

done:
    if (fflush(reportout) != 0) {
        set_errmsg(errmsg, errmsgsize, "report write failed");
        ret = 1;
    }
    free(verts); free(faces); free(edges); free(J);
    return ret;
}

static void chomp(char *s)
{
    size_t n = strlen(s);
    while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r')) s[--n] = '\0';
}

static uint64_t fnv1a64(const char *s)
{
    uint64_t h = UINT64_C(1469598103934665603);
    while (*s) {
        h ^= (unsigned char)*s++;
        h *= UINT64_C(1099511628211);
    }
    return h;
}

static int ensure_directory_exists(const char *path)
{
    struct stat st;
    if (stat(path, &st) == 0) return S_ISDIR(st.st_mode) ? 0 : -1;
    if (errno == ENOENT && mkdir(path, 0777) == 0) return 0;
    return -1;
}

static int make_output_paths(const char *outdir, long index, uint64_t hash,
                             char *finalpath, size_t finalsize,
                             char *tmppath, size_t tmpsize)
{
    int n1 = snprintf(finalpath, finalsize, "%s/n%08ld_%016llx.report",
                      outdir, index, (unsigned long long)hash);
    int n2 = snprintf(tmppath, tmpsize, "%s/.tmp_n%08ld_%016llx_%ld.report",
                      outdir, index, (unsigned long long)hash, (long)getpid());
    if (n1 < 0 || n2 < 0) return -1;
    if ((size_t)n1 >= finalsize || (size_t)n2 >= tmpsize) return -1;
    return 0;
}

static int run_batch(const char *outdir)
{
    if (ensure_directory_exists(outdir) < 0) {
        fprintf(stderr, "cannot use output directory %s: %s\n", outdir, strerror(errno));
        return 1;
    }

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    long index = 0;
    int saw_reject = 0;
    int saw_fail = 0;

    while ((len = getline(&line, &cap, stdin)) != -1) {
        (void)len;
        chomp(line);
        if (line[0] == '\0') continue;

        uint64_t hash = fnv1a64(line);
        char finalpath[4096], tmppath[4096], errmsg[1024];
        errmsg[0] = '\0';

        if (make_output_paths(outdir, index, hash, finalpath, sizeof(finalpath),
                              tmppath, sizeof(tmppath)) < 0) {
            snprintf(errmsg, sizeof(errmsg), "output path too long");
            fprintf(stdout, "%ld\tfail\t%s\t\t%s\n", index, line, errmsg);
            saw_fail = 1;
            index++;
            continue;
        }

        FILE *fp = fopen(tmppath, "w");
        if (!fp) {
            snprintf(errmsg, sizeof(errmsg), "could not open temp report: %s", strerror(errno));
            fprintf(stdout, "%ld\tfail\t%s\t\t%s\n", index, line, errmsg);
            saw_fail = 1;
            index++;
            continue;
        }

        int rc = prove_one_obj(line, fp, errmsg, sizeof(errmsg));
        int close_failed = 0;
        if (fflush(fp) != 0) close_failed = 1;
        if (fclose(fp) != 0) close_failed = 1;

        if (rc != 1 && close_failed) {
            snprintf(errmsg, sizeof(errmsg), "could not close temp report cleanly: %s", strerror(errno));
            rc = 1;
        }

        if (rc == 0 || rc == 2) {
            if (rename(tmppath, finalpath) != 0) {
                snprintf(errmsg, sizeof(errmsg), "rename failed: %s", strerror(errno));
                remove(tmppath);
                fprintf(stdout, "%ld\tfail\t%s\t\t%s\n", index, line, errmsg);
                saw_fail = 1;
            } else if (rc == 0) {
                fprintf(stdout, "%ld\taccept\t%s\t%s\t\n", index, line, finalpath);
            } else {
                fprintf(stdout, "%ld\treject\t%s\t%s\t%s\n", index, line, finalpath,
                        errmsg[0] ? errmsg : "final REJECT");
                saw_reject = 1;
            }
        } else {
            remove(tmppath);
            if (errmsg[0] == '\0') snprintf(errmsg, sizeof(errmsg), "prover failed");
            fprintf(stderr, "[%ld] %s\n", index, errmsg);
            fprintf(stdout, "%ld\tfail\t%s\t\t%s\n", index, line, errmsg);
            saw_fail = 1;
        }

        index++;
    }

    free(line);
    if (ferror(stdin)) {
        fprintf(stderr, "error while reading stdin\n");
        return 1;
    }
    if (saw_fail) return 1;
    if (saw_reject) return 2;
    return 0;
}

static void usage(const char *prog)
{
    fprintf(stderr, "usage:\n");
    fprintf(stderr, "  %s OBJFILE\n", prog);
    fprintf(stderr, "  %s --batch --outdir DIR < objfiles.txt\n", prog);
}

/* ============================================================
 *   main
 * ============================================================ */

int main(int argc, char **argv) {
    if (argc == 2 && strcmp(argv[1], "--batch") != 0) {
        char errmsg[1024];
        int rc = prove_one_obj(argv[1], stdout, errmsg, sizeof(errmsg));
        if (rc == 1) fprintf(stderr, "error: %s\n", errmsg[0] ? errmsg : "prover failed");
        return rc;
    }

    if (argc == 4 && strcmp(argv[1], "--batch") == 0 && strcmp(argv[2], "--outdir") == 0) {
        return run_batch(argv[3]);
    }

    usage(argv[0]);
    return 1;
}
