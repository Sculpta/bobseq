/* Molecule definition used by the benchmark (definitions C, D, E).
 *
 * A molecule is (well, contig, strand, UMI cluster at edit distance <= MAXED, reads within TOL bp
 * of the molecule's START). The window is anchored at the molecule start, not sliding, so a run
 * of reads cannot chain transitively across a whole gene. MAXED=0 is exact-UMI; MAXED=1 merges a
 * UMI one substitution away from a molecule already open in the same (well, contig, strand) block.
 *
 * UMIs up to 32 nt (2 bits per base in a 64-bit key) are held in an open-addressing hash table
 * that is reset per block by a generation counter, so a 32N UMI (the poly-dT libraries) fits
 * without a 4^L index array.
 *
 * Input  (stdin) : well \t umi \t contig \t pos \t gene [\t strand]
 *                  sorted by well, contig, strand, pos; pos is the read's 5'-end anchor
 * Output (stdout): well \t gene                           one row per molecule
 * Usage: dedup MAXED TOL [perwell.tsv] [hist.tsv]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXL 32
typedef struct { uint64_t key; int gen; int anchor; long long molidx; } slot_t;
static slot_t *tab; static uint64_t tabsz = 1ULL << 20, tabmask; static uint64_t used = 0;
static int curgen = 1;

static inline uint64_t hsh(uint64_t k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return k; }

static void tab_init(uint64_t n) {
    tabsz = n; tabmask = n - 1; used = 0;
    tab = calloc(tabsz, sizeof(slot_t));
    if (!tab) { fprintf(stderr, "alloc failed\n"); exit(1); }
}
static slot_t *tab_find(uint64_t key) {            /* live slot for key in this generation, or NULL */
    uint64_t i = hsh(key) & tabmask;
    for (;;) {
        slot_t *s = &tab[i];
        if (s->gen == 0) return NULL;                /* never used */
        if (s->key == key) return s->gen == curgen ? s : NULL;   /* stale = not in this block */
        i = (i + 1) & tabmask;
    }
}
static void tab_grow(void);
static slot_t *tab_put(uint64_t key) {             /* slot to (over)write for key */
    if ((used + 1) * 2 > tabsz) tab_grow();
    uint64_t i = hsh(key) & tabmask;
    for (;;) {
        slot_t *s = &tab[i];
        if (s->gen == 0) { used++; s->key = key; return s; }
        if (s->key == key) return s;                 /* reuse (stale or live) */
        i = (i + 1) & tabmask;
    }
}
static void tab_grow(void) {
    slot_t *old = tab; uint64_t oldsz = tabsz;
    tab_init(tabsz * 2);
    for (uint64_t i = 0; i < oldsz; i++) if (old[i].gen == curgen) { slot_t *s = tab_put(old[i].key); *s = old[i]; }
    free(old);
}

static inline int enc(const char *s, int L, uint64_t *out) {
    uint64_t c = 0;
    for (int i = 0; i < L; i++) {
        int b;
        switch (s[i]) { case 'A': b=0; break; case 'C': b=1; break;
                        case 'G': b=2; break; case 'T': b=3; break; default: return -1; }
        c = (c << 2) | (uint64_t)b;
    }
    *out = c; return 0;
}

int main(int argc, char **argv) {
    int MAXED = argc > 1 ? atoi(argv[1]) : 1;
    int TOL   = argc > 2 ? atoi(argv[2]) : 1000;
    const char *PERWELL = argc > 3 ? argv[3] : NULL;
    const char *HIST    = argc > 4 ? argv[4] : NULL;   /* reads-per-molecule histogram */

    char *line = NULL; size_t cap = 0; ssize_t len;
    int L = -1;
    char pwell[64] = {0}, pcontig[64] = {0}, pstrand[8] = {0};
    long long nread = 0, nmol = 0, nmerge_exact = 0, nmerge_nbr = 0, nskip = 0;

    char (*wname)[64] = malloc(600000 * 64); long long *wcount = calloc(600000, 8);
    long long *hist = calloc(4096, 8);
    int *rcount = NULL; long long ncur = 0, capcur = 0;
    int nw = 0; char lastw[64] = {0};
    tab_init(tabsz);

    while ((len = getline(&line, &cap, stdin)) > 0) {
        char *w = line, *t1 = strchr(w, '\t');  if (!t1) continue;  *t1 = 0;
        char *u = t1 + 1, *t2 = strchr(u, '\t'); if (!t2) continue;  *t2 = 0;
        char *c = t2 + 1, *t3 = strchr(c, '\t'); if (!t3) continue;  *t3 = 0;
        char *p = t3 + 1, *t4 = strchr(p, '\t'); if (!t4) continue;  *t4 = 0;
        char *g = t4 + 1, *t5 = strchr(g, '\t');
        char *sd = "";
        if (t5) { *t5 = 0; sd = t5 + 1; }
        char *nl = strchr(t5 ? sd : g, '\n'); if (nl) *nl = 0;

        /* UMI length is a property of the WELL (constant within one well, different between the BOBseq code sets:
         * 8 / 11 / 14 nt with the 6N+spacer UMI); input is sorted by well, so re-derive it whenever the well
         * changes. Taking one length from the first row would silently skip every well of another length. */
        if (L < 0 || strcmp(w, lastw)) {
            L = (int)strlen(u);
            if (L < 4 || L > MAXL) { fprintf(stderr, "bad UMI length %d (4..%d) in well %s\n", L, MAXL, w); return 1; }
        }
        if ((int)strlen(u) != L) { nskip++; continue; }
        uint64_t code;
        if (enc(u, L, &code) < 0) { nskip++; continue; }   /* N in UMI */
        int pos = atoi(p);
        nread++;

        if (strcmp(w, pwell) || strcmp(c, pcontig) || strcmp(sd, pstrand)) {
            curgen++;                                   /* new (well, contig, strand) block: all slots stale */
            snprintf(pwell, 64, "%s", w); snprintf(pcontig, 64, "%s", c);
            snprintf(pstrand, 8, "%s", sd);
        }
        if (strcmp(w, lastw)) {
            snprintf(wname[nw], 64, "%s", w); nw++; snprintf(lastw, 64, "%s", w);
        }

        int hit = -1; long long hitidx = -1;
        slot_t *s = tab_find(code);
        if (s && pos - s->anchor <= TOL) { hit = s->anchor; hitidx = s->molidx; nmerge_exact++; }
        else if (MAXED > 0) {
            for (int i = 0; i < L && hit < 0; i++) {
                int shift = 2 * (L - 1 - i); int cur = (int)((code >> shift) & 3ULL);
                for (int b = 0; b < 4; b++) {
                    if (b == cur) continue;
                    uint64_t nb = (code & ~(3ULL << shift)) | ((uint64_t)b << shift);
                    slot_t *t = tab_find(nb);
                    if (t && pos - t->anchor <= TOL) { hit = t->anchor; hitidx = t->molidx; nmerge_nbr++; break; }
                }
            }
        }
        if (hit < 0) {                                  /* a new molecule */
            hit = pos; nmol++; wcount[nw - 1]++;
            printf("%s\t%s\n", w, g);
            if (ncur == capcur) { capcur = capcur ? capcur * 2 : 1 << 20; rcount = realloc(rcount, capcur * sizeof(int)); }
            rcount[ncur++] = 0;
            hitidx = ncur - 1;
        }
        slot_t *e = tab_put(code); e->gen = curgen; e->anchor = hit; e->molidx = hitidx;
        if (hitidx >= 0 && hitidx < ncur) rcount[hitidx]++;
    }
    for (long long i = 0; i < ncur; i++) { int k = rcount[i]; if (k < 1) k = 1; if (k > 4095) k = 4095; hist[k]++; }
    if (HIST) { FILE *f = fopen(HIST, "w"); fprintf(f, "reads_per_molecule\tmolecules\n");
        for (int k = 1; k < 4096; k++) if (hist[k]) fprintf(f, "%d\t%lld\n", k, hist[k]); fclose(f); }
    if (PERWELL) { FILE *f = fopen(PERWELL, "w"); fprintf(f, "well\tmolecules\n");
        for (int i = 0; i < nw; i++) if (wcount[i]) fprintf(f, "%s\t%lld\n", wname[i], wcount[i]); fclose(f); }
    fprintf(stderr, "  MAXED=%d TOL=%d\n  reads %lld | molecules %lld | wells %d\n"
                    "  merged: exact %lld, ed<=1 neighbour %lld | skipped %lld\n",
            MAXED, TOL, nread, nmol, nw, nmerge_exact, nmerge_nbr, nskip);
    return 0;
}
