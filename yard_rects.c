/* yard_rects.c
 *
 * Updates from previous patch:
 *  - rect_completely_within uses the half-space test (no ray-cast)
 *  - quad vertex winding enforced CCW (reversed if clockwise)
 *  - quad_axes/offsets computed robustly as min dot over vertices (order-insensitive)
 *  - removed redundant read/write of thread_seeds in pure fitness evaluation loop
 *
 * Compile:
 *   gcc -O3 -march=native -fopenmp yard_rects.c -o yard_rects -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <stdint.h>

#define POPULATION_SIZE 1000
#define MAX_GENERATIONS 20000
#define CONVERGENCE_THRESHOLD 400
#define TOP_K_DIVISOR 10
#define ELITE_COUNT 8
#define RANDOM_INJECTION_RATE 0.08f
#define INITIAL_SIGMA_SCALE 1.0f
#define MIN_SIGMA 0.12f
#define ANNEALING_END_SCALE 0.6f
#define REPAIR_SHRINK_FACTOR 0.96f
#define REPAIR_MAX_STEPS 18
#define HILL_CLIMB_STEPS 10
#define HILL_PERTURB_SCALE 0.06f

typedef struct { float x,y; } Point;
typedef struct { Point vertices[4]; } IrregularQuadrilateral;
typedef struct { Point center; float width,height,angle; } RotatedRectangle;

typedef struct {
    RotatedRectangle rect;
    Point corners[4];
    float minx, maxx, miny, maxy;
} SmallRectPre;

/* ----- fast thread-local xorshift RNG ----- */
static inline uint32_t xorshift32(uint32_t *state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return *state = x ? x : 0x9e3779b1u;
}
static inline float rndf(uint32_t *state, float a, float b) {
    uint32_t r = xorshift32(state);
    return a + ((r / (float)UINT32_MAX) * (b - a));
}
static inline float normal_approx(uint32_t *state) {
    float s = 0.f;
    for (int i=0;i<6;i++) s += rndf(state, -1.f, 1.f);
    return s * 0.40824829f;
}

/* geometry helpers */
static inline float rect_area(const RotatedRectangle *r) { return r->width * r->height; }

/* compute rectangle corners */
static inline void get_rectangle_corners(const RotatedRectangle *rect, Point corners[4]) {
    float ca = cosf(rect->angle), sa = sinf(rect->angle);
    float hw = rect->width * 0.5f, hh = rect->height * 0.5f;
    corners[0].x = rect->center.x + (-hw * ca - -hh * sa); corners[0].y = rect->center.y + (-hw * sa + -hh * ca);
    corners[1].x = rect->center.x + ( hw * ca - -hh * sa); corners[1].y = rect->center.y + ( hw * sa + -hh * ca);
    corners[2].x = rect->center.x + ( hw * ca -  hh * sa); corners[2].y = rect->center.y + ( hw * sa +  hh * ca);
    corners[3].x = rect->center.x + (-hw * ca -  hh * sa); corners[3].y = rect->center.y + (-hw * sa +  hh * ca);
}

/* compute AABB for a rotated rect given its corners */
static inline void compute_aabb_from_corners(const Point c[4], float *minx, float *maxx, float *miny, float *maxy) {
    *minx = *maxx = c[0].x; *miny = *maxy = c[0].y;
    for (int i=1;i<4;i++){
        if (c[i].x < *minx) *minx = c[i].x;
        if (c[i].x > *maxx) *maxx = c[i].x;
        if (c[i].y < *miny) *miny = c[i].y;
        if (c[i].y > *maxy) *maxy = c[i].y;
    }
}

/* AABB overlap test */
static inline bool aabb_overlap(float a_minx,float a_maxx,float a_miny,float a_maxy,
                                float b_minx,float b_maxx,float b_miny,float b_maxy) {
    if (a_maxx < b_minx || b_maxx < a_minx) return false;
    if (a_maxy < b_miny || b_maxy < a_miny) return false;
    return true;
}

/* SAT overlap test (axes are unnormalized; projections still work) */
static inline void project_onto_axis(const Point verts[], int n, Point axis, float *minp, float *maxp) {
    float p = verts[0].x*axis.x + verts[0].y*axis.y;
    *minp = *maxp = p;
    for (int i=1;i<n;i++){
        p = verts[i].x*axis.x + verts[i].y*axis.y;
        if (p < *minp) *minp = p;
        if (p > *maxp) *maxp = p;
    }
}
bool rectangles_overlap_sat(const RotatedRectangle *r1, const RotatedRectangle *r2, const Point r1c[4], const Point r2c[4]) {
    Point axes[4];
    axes[0].x = -(r1c[1].y - r1c[0].y); axes[0].y = (r1c[1].x - r1c[0].x);
    axes[1].x = -(r1c[2].y - r1c[1].y); axes[1].y = (r1c[2].x - r1c[1].x);
    axes[2].x = -(r2c[1].y - r2c[0].y); axes[2].y = (r2c[1].x - r2c[0].x);
    axes[3].x = -(r2c[2].y - r2c[1].y); axes[3].y = (r2c[2].x - r2c[1].x);
    for (int i=0;i<4;i++){
        float min1,max1,min2,max2;
        project_onto_axis(r1c,4,axes[i],&min1,&max1);
        project_onto_axis(r2c,4,axes[i],&min2,&max2);
        if (max1 < min2 || max2 < min1) return false;
    }
    return true;
}

/* half-space point-in-quad test using precomputed quad axes and offsets */
static inline bool point_in_quad_fast(const Point *p, const Point quad_axes[4], const float quad_offsets[4]) {
    for (int i=0;i<4;i++){
        if (quad_axes[i].x * p->x + quad_axes[i].y * p->y < quad_offsets[i]) return false;
    }
    return true;
}

/* check fully inside large quad with fast AABB precheck using point_in_quad_fast */
bool rect_completely_within(const RotatedRectangle *r, const IrregularQuadrilateral *quad, float quad_minx,float quad_maxx,float quad_miny,float quad_maxy, const Point quad_axes[4], const float quad_offsets[4]) {
    Point corners[4];
    get_rectangle_corners(r, corners);
    float rminx,rmaxx,rminy,rmaxy;
    compute_aabb_from_corners(corners,&rminx,&rmaxx,&rminy,&rmaxy);
    if (rminx < quad_minx || rmaxx > quad_maxx || rminy < quad_miny || rmaxy > quad_maxy) return false;
    for (int i=0;i<4;i++) if (!point_in_quad_fast(&corners[i], quad_axes, quad_offsets)) return false;
    return true;
}

/* Evaluate fitness: area if valid; uses AABB prechecks vs small rects
   New signature accepts quad_axes and quad_offsets for fast point-in-quad. */
float evaluate_fitness_fast(const RotatedRectangle *r, const IrregularQuadrilateral *quad,
                            const SmallRectPre smalls[], int nsm,
                            float quad_minx,float quad_maxx,float quad_miny,float quad_maxy,
                            const Point quad_axes[4], const float quad_offsets[4]) {
    if (r->width <= 0.f || r->height <= 0.f) return 0.f;
    Point corners[4];
    get_rectangle_corners(r, corners);
    float rminx,rmaxx,rminy,rmaxy;
    compute_aabb_from_corners(corners,&rminx,&rmaxx,&rminy,&rmaxy);
    if (rminx < quad_minx || rmaxx > quad_maxx || rminy < quad_miny || rmaxy > quad_maxy) return 0.f;
    for (int i=0;i<4;i++) if (!point_in_quad_fast(&corners[i], quad_axes, quad_offsets)) return 0.f;
    for (int i=0;i<nsm;i++){
        if (!aabb_overlap(rminx,rmaxx,rminy,rmaxy, smalls[i].minx,smalls[i].maxx,smalls[i].miny,smalls[i].maxy)) continue;
        if (rectangles_overlap_sat(r, &smalls[i].rect, corners, (Point*)smalls[i].corners)) return 0.f;
    }
    return rect_area(r);
}

/* repair (shrink + small nudges), lighter than before
   Now forwards quad axes/offsets into evaluate_fitness_fast */
RotatedRectangle repair_light(RotatedRectangle r, const IrregularQuadrilateral *quad,
                             const SmallRectPre smalls[], int nsm,
                             float quad_minx,float quad_maxx,float quad_miny,float quad_maxy,
                             const Point quad_axes[4], const float quad_offsets[4],
                             uint32_t *rng_state) {
    float nudge = fmaxf(0.3f, fminf(r.width, r.height)*0.18f);
    for (int dx=-1; dx<=1; dx++){
        for (int dy=-1; dy<=1; dy++){
            RotatedRectangle c = r;
            c.center.x += dx * nudge * rndf(rng_state, 0.0f, 1.0f);
            c.center.y += dy * nudge * rndf(rng_state, 0.0f, 1.0f);
            float f = evaluate_fitness_fast(&c, quad, smalls, nsm, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets);
            if (f > 0.f) return c;
        }
    }
    RotatedRectangle s = r;
    for (int step=0; step < REPAIR_MAX_STEPS; step++){
        s.width *= REPAIR_SHRINK_FACTOR;
        s.height *= REPAIR_SHRINK_FACTOR;
        s.center.x += rndf(rng_state, -0.3f, 0.3f);
        s.center.y += rndf(rng_state, -0.3f, 0.3f);
        float f = evaluate_fitness_fast(&s, quad, smalls, nsm, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets);
        if (f > 0.f) return s;
    }
    s.width = s.height = 0.f;
    return s;
}

/* light hill-climb (very few steps) */
RotatedRectangle hill_climb_light(RotatedRectangle r, const IrregularQuadrilateral *quad,
                                 const SmallRectPre smalls[], int nsm,
                                 float quad_minx,float quad_maxx,float quad_miny,float quad_maxy,
                                 const Point quad_axes[4], const float quad_offsets[4],
                                 uint32_t *rng_state) {
    float bestf = evaluate_fitness_fast(&r, quad, smalls, nsm, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets);
    RotatedRectangle best = r;
    for (int i=0;i<HILL_CLIMB_STEPS;i++){
        RotatedRectangle c = best;
        c.center.x += normal_approx(rng_state) * (HILL_PERTURB_SCALE * best.width);
        c.center.y += normal_approx(rng_state) * (HILL_PERTURB_SCALE * best.height);
        c.width += normal_approx(rng_state) * (HILL_PERTURB_SCALE * best.width);
        c.height += normal_approx(rng_state) * (HILL_PERTURB_SCALE * best.height);
        c.angle += normal_approx(rng_state) * 0.02f;
        if (c.width <= 0.1f || c.height <= 0.1f) continue;
        c = repair_light(c, quad, smalls, nsm, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets, rng_state);
        float f = evaluate_fitness_fast(&c, quad, smalls, nsm, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets);
        if (f > bestf) { best = c; bestf = f; }
    }
    return best;
}

/* SVG write */
void generate_svg(const char *fn, const IrregularQuadrilateral *q, const SmallRectPre smalls[], int nsm, const RotatedRectangle *best) {
    FILE *f = fopen(fn, "w");
    if (!f) return;
    float minx=q->vertices[0].x, maxx=q->vertices[0].x, miny=q->vertices[0].y, maxy=q->vertices[0].y;
    for (int i=1;i<4;i++){ if(q->vertices[i].x<minx)minx=q->vertices[i].x; if(q->vertices[i].x>maxx)maxx=q->vertices[i].x;
                          if(q->vertices[i].y<miny)miny=q->vertices[i].y; if(q->vertices[i].y>maxy)maxy=q->vertices[i].y; }
    fprintf(f,"<svg xmlns='http://www.w3.org/2000/svg' width='800' height='800' viewBox='%.2f %.2f %.2f %.2f'>\n",
            minx-2,miny-2,(maxx-minx)+4,(maxy-miny)+4);
    fprintf(f,"<style> .q{fill:none;stroke:black;stroke-width:0.6;} .sm{fill:rgba(255,0,0,0.6);} .b{fill:rgba(0,255,0,0.4);} </style>\n");
    fprintf(f,"<polygon class='q' points='");
    for (int i=0;i<4;i++) fprintf(f,"%.2f,%.2f ", q->vertices[i].x, q->vertices[i].y);
    fprintf(f,"' />\n");
    for (int i=0;i<nsm;i++){ fprintf(f,"<polygon class='sm' points='"); for (int j=0;j<4;j++) fprintf(f,"%.2f,%.2f ", smalls[i].corners[j].x, smalls[i].corners[j].y); fprintf(f,"' />\n"); }
    Point bc[4]; get_rectangle_corners(best, bc); fprintf(f,"<polygon class='b' points='"); for(int j=0;j<4;j++) fprintf(f,"%.2f,%.2f ",bc[j].x,bc[j].y); fprintf(f,"' />\n");
    fprintf(f,"</svg>\n"); fclose(f);
    printf("SVG -> %s\n", fn);
}

/* helper: compute signed area to detect winding (positive = CCW) */
static inline float polygon_signed_area_ccw(const Point v[4]) {
    float a = 0.f;
    for (int i=0;i<4;i++){
        int j = (i+1)&3;
        a += v[i].x * v[j].y - v[j].x * v[i].y;
    }
    return a * 0.5f;
}

/* ---------------- main program ---------------- */
int main() {
    srand((unsigned)time(NULL));
    IrregularQuadrilateral quad = { .vertices = {{0,0},{79.224f,70.317f},{55.0f,129.33f},{0.55f,123.627f}} };
    RotatedRectangle smalls_in[] = {
// This solution assumes that bent dead tree is gone.
        {{64.311f,  77.657f},  2.700f,  2.700f, 0.000f},    // tree1, nearest house.
        {{40.365f,  81.172f},  3.354f,  3.354f, 0.000f},    // tree2, middle tree.
//      {{21.601f, 101.087f},  3.354f,  3.354f, 0.000f},    // tree3, bent tree.
        {{35.382f, 115.751f}, 20.000f, 12.000f, 0.101f},    // chicken coop
        {{19.550f,  71.470f}, 37.970f, 74.670f, 3.136f},    // biggest rectangle solution
//      {{52.470f,  83.190f}, 54.510f, 20.820f, 1.579f},    // 2nd biggest rectangle solution
//      {{12.690f, 116.280f}, 14.670f, 24.240f, 1.57f}      // 3rd biggest rectangle solution
//      {{ 9.570f,  25.530f}, 18.830f, 17.310f, 6.277f}     // 4th biggest rectangle solution
    };
    const int num_small = sizeof(smalls_in)/sizeof(smalls_in[0]);

    /* ensure quad is CCW; reverse if necessary */
    if (polygon_signed_area_ccw(quad.vertices) < 0.f) {
        Point tmp[4];
        for (int i=0;i<4;i++) tmp[i] = quad.vertices[3-i];
        for (int i=0;i<4;i++) quad.vertices[i] = tmp[i];
    }

    /* precompute small rect corners and aabbs */
    SmallRectPre smalls[16];
    for (int i=0;i<num_small;i++){
        smalls[i].rect = smalls_in[i];
        get_rectangle_corners(&smalls[i].rect, smalls[i].corners);
        compute_aabb_from_corners(smalls[i].corners, &smalls[i].minx, &smalls[i].maxx, &smalls[i].miny, &smalls[i].maxy);
    }

    /* quad bounding box */
    float quad_minx = quad.vertices[0].x, quad_maxx = quad.vertices[0].x, quad_miny = quad.vertices[0].y, quad_maxy = quad.vertices[0].y;
    for (int i=1;i<4;i++){ if (quad.vertices[i].x < quad_minx) quad_minx = quad.vertices[i].x; if (quad.vertices[i].x > quad_maxx) quad_maxx = quad.vertices[i].x;
                          if (quad.vertices[i].y < quad_miny) quad_miny = quad.vertices[i].y; if (quad.vertices[i].y > quad_maxy) quad_maxy = quad.vertices[i].y; }
    float max_dim = fmaxf(quad_maxx-quad_minx, quad_maxy-quad_miny);

    /* precompute quad half-space axes and offsets robustly:
       compute axis = (-edge.y, edge.x) then offsets = min(dot(axis, vertex_j)) so ordering doesn't matter */
    Point quad_axes[4];
    float quad_offsets[4];
    for (int i=0;i<4;i++){
        Point a = quad.vertices[i];
        Point b = quad.vertices[(i+1)&3];
        Point edge = { b.x - a.x, b.y - a.y };
        quad_axes[i].x = -edge.y;
        quad_axes[i].y =  edge.x;
        /* compute min dot over all vertices so offset is robust to ordering */
        float mo = quad_axes[i].x * quad.vertices[0].x + quad_axes[i].y * quad.vertices[0].y;
        for (int j=1;j<4;j++){
            float d = quad_axes[i].x * quad.vertices[j].x + quad_axes[i].y * quad.vertices[j].y;
            if (d < mo) mo = d;
        }
        quad_offsets[i] = mo;
    }

    /* allocate populations */
    RotatedRectangle *pop = malloc(sizeof(RotatedRectangle)*POPULATION_SIZE);
    RotatedRectangle *next = malloc(sizeof(RotatedRectangle)*POPULATION_SIZE);
    float *fitness = malloc(sizeof(float)*POPULATION_SIZE);
    static int idx[POPULATION_SIZE];

    /* prepare per-thread RNG seeds once */
    int omp_max_threads = omp_get_max_threads();
    uint32_t *thread_seeds = malloc(sizeof(uint32_t) * omp_max_threads);
    uint32_t seed_base = (uint32_t)time(NULL) ^ 0xA5A5A5A5u;
    for (int t=0; t<omp_max_threads; ++t) thread_seeds[t] = seed_base ^ (uint32_t)(t * 2654435761u);

    /* init pop using per-thread seeds */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        uint32_t st = thread_seeds[tid];
        #pragma omp for schedule(static)
        for (int i=0;i<POPULATION_SIZE;i++){
            pop[i].center.x = rndf(&st, quad_minx, quad_maxx);
            pop[i].center.y = rndf(&st, quad_miny, quad_maxy);
            pop[i].width    = rndf(&st, 1.0f, max_dim*0.5f);
            pop[i].height   = rndf(&st, 1.0f, max_dim*0.5f);
            pop[i].angle    = rndf(&st, 0.f, 2.f*M_PI);
        }
        thread_seeds[tid] = st;
    }

    RotatedRectangle best = pop[0]; float bestf = 0.f; int gens_no_improve = 0;
    int TOP_K = POPULATION_SIZE / TOP_K_DIVISOR; if (TOP_K < 2) TOP_K = 2;

    for (int gen=0; gen < MAX_GENERATIONS; gen++) {
        /* evaluate fitness in parallel
           NOTE: evaluate_fitness_fast does not use RNG, so we do not touch per-thread seed here */
        #pragma omp parallel for schedule(static)
        for (int i=0;i<POPULATION_SIZE;i++){
            fitness[i] = evaluate_fitness_fast(&pop[i], &quad, smalls, num_small, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets);
        }

        /* compute top-K indices (partial selection sort) */
        for (int i=0;i<POPULATION_SIZE;i++) idx[i] = i;
        for (int a=0;a<TOP_K;a++){
            int bi = a;
            for (int b=a+1;b<POPULATION_SIZE;b++) if (fitness[idx[b]] > fitness[idx[bi]]) bi = b;
            int tmp = idx[a]; idx[a] = idx[bi]; idx[bi] = tmp;
        }

        /* update best */
        for (int t=0;t<TOP_K;t++){
            int i = idx[t];
            if (fitness[i] > bestf) { bestf = fitness[i]; best = pop[i]; gens_no_improve = 0; }
        }
        gens_no_improve++;
        if (gens_no_improve > CONVERGENCE_THRESHOLD) { printf("Converged at gen %d best area %.3f\n", gen, rect_area(&best)); break; }

        /* compute allele means & std (angle via circular) */
        float mean_cx=0, mean_cy=0, mean_w=0, mean_h=0, sum_sin=0, sum_cos=0;
        for (int t=0;t<TOP_K;t++){
            RotatedRectangle *r = &pop[idx[t]];
            mean_cx += r->center.x; mean_cy += r->center.y; mean_w += r->width; mean_h += r->height;
            sum_sin += sinf(r->angle); sum_cos += cosf(r->angle);
        }
        mean_cx/=TOP_K; mean_cy/=TOP_K; mean_w/=TOP_K; mean_h/=TOP_K;
        float mean_angle = atan2f(sum_sin, sum_cos);
        float var_cx=0,var_cy=0,var_w=0,var_h=0,var_a=0;
        for (int t=0;t<TOP_K;t++){
            RotatedRectangle *r = &pop[idx[t]];
            float dx = r->center.x - mean_cx; var_cx += dx*dx;
            float dy = r->center.y - mean_cy; var_cy += dy*dy;
            float dw = r->width - mean_w; var_w += dw*dw;
            float dh = r->height - mean_h; var_h += dh*dh;
            float da = r->angle - mean_angle; while (da > M_PI) da -= 2*M_PI; while (da < -M_PI) da += 2*M_PI;
            var_a += da*da;
        }
        float std_cx = sqrtf(var_cx / TOP_K) * INITIAL_SIGMA_SCALE + MIN_SIGMA;
        float std_cy = sqrtf(var_cy / TOP_K) * INITIAL_SIGMA_SCALE + MIN_SIGMA;
        float std_w  = sqrtf(var_w  / TOP_K) * INITIAL_SIGMA_SCALE + MIN_SIGMA;
        float std_h  = sqrtf(var_h  / TOP_K) * INITIAL_SIGMA_SCALE + MIN_SIGMA;
        float std_a  = sqrtf(var_a  / TOP_K) * INITIAL_SIGMA_SCALE + MIN_SIGMA;

        /* anneal sigma */
        float tfrac = (float)gen / (float)MAX_GENERATIONS;
        float anneal = 1.f - tfrac * (1.f - ANNEALING_END_SCALE);
        std_cx *= anneal; std_cy *= anneal; std_w *= anneal; std_h *= anneal; std_a *= anneal;

        /* elitism: copy top ELITE_COUNT */
        int elite = ELITE_COUNT; if (elite > POPULATION_SIZE) elite = POPULATION_SIZE;
        for (int e=0;e<elite;e++) next[e] = pop[idx[e]];

        /* sample rest in parallel using per-thread seeds */
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            uint32_t st = thread_seeds[tid];
            #pragma omp for schedule(static)
            for (int i=elite;i<POPULATION_SIZE;i++){
                if (rndf(&st, 0.f, 1.f) < RANDOM_INJECTION_RATE) {
                    next[i].center.x = rndf(&st, quad_minx, quad_maxx);
                    next[i].center.y = rndf(&st, quad_miny, quad_maxy);
                    next[i].width    = rndf(&st, 1.0f, max_dim*0.5f);
                    next[i].height   = rndf(&st, 1.0f, max_dim*0.5f);
                    next[i].angle    = rndf(&st, 0.f, 2.f*M_PI);
                    continue;
                }
                RotatedRectangle child;
                child.center.x = mean_cx + normal_approx(&st) * std_cx;
                child.center.y = mean_cy + normal_approx(&st) * std_cy;
                child.width    = mean_w  + normal_approx(&st) * std_w;
                child.height   = mean_h  + normal_approx(&st) * std_h;
                child.angle    = mean_angle + normal_approx(&st) * std_a;
                if (child.center.x < quad_minx) child.center.x = quad_minx + rndf(&st, 0.f, std_cx);
                if (child.center.x > quad_maxx) child.center.x = quad_maxx - rndf(&st, 0.f, std_cx);
                if (child.center.y < quad_miny) child.center.y = quad_miny + rndf(&st, 0.f, std_cy);
                if (child.center.y > quad_maxy) child.center.y = quad_maxy - rndf(&st, 0.f, std_cy);
                if (child.width < 0.4f) child.width = 0.4f + rndf(&st, 0.f, 0.6f);
                if (child.height < 0.4f) child.height = 0.4f + rndf(&st, 0.f, 0.6f);
                child.angle = fmodf(child.angle, 2.f*(float)M_PI);
                if (child.angle < 0.f) child.angle += 2.f*(float)M_PI;

                RotatedRectangle repaired = repair_light(child, &quad, smalls, num_small, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets, &st);
                if (repaired.width <= 0.f || repaired.height <= 0.f) {
                    RotatedRectangle fallback = pop[idx[(int)(rndf(&st, 0.f, (float)TOP_K-1.f))]];
                    repaired = repair_light(fallback, &quad, smalls, num_small, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets, &st);
                    if (repaired.width <= 0.f || repaired.height <= 0.f) { repaired.width = repaired.height = 0.f; }
                }
                if (repaired.width > 0.f && repaired.height > 0.f) {
                    repaired = hill_climb_light(repaired, &quad, smalls, num_small, quad_minx, quad_maxx, quad_miny, quad_maxy, quad_axes, quad_offsets, &st);
                }
                next[i] = repaired;
            }
            thread_seeds[tid] = st;
        } // end parallel

        /* swap */
        RotatedRectangle *tmp = pop; pop = next; next = tmp;

        if ((gen & 255) == 0) {
            printf("gen %d best_area %.3f topk %d stds cx%.3f w%.3f a%.3f\n", gen, rect_area(&best), TOP_K, std_cx, std_w, std_a);
        }
    } // end gen loop

    printf("BEST area %.3f center %.2f %.2f w %.2f h %.2f angle %.3f\n", rect_area(&best), best.center.x, best.center.y, best.width, best.height, best.angle);
    generate_svg("output_yard_rects.svg", &quad, smalls, num_small, &best);

    free(pop); free(next); free(fitness); free(thread_seeds);
    return 0;
}

