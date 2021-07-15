/*
 * provide means for multi-dimensional distance class based soft constraints
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/constraints/soft.h>

#include "dist_class_sc.h"


static size_t *
compute_distancies(int              i,
                   int              j,
                   int              k,
                   int              l,
                   unsigned char    decomp,
                   sc_dist_class_t  *d)


// rechne mir induzierte Distanzen, 
// konvertiere das ganze in Energien


// über alle Referenzstrukturen (A und B)
// berechnet die induzierte Distanz, wenn man hairpin / multibranch zerlegung,
// mit koordinaten i und j, l (Positionen des einschließenden bps)
// k,l: Koordinaten von Subkordinaten, kann auch Ende eines Bereiches sein 
// (Intervall i bis j in 2 Teile zerlegen)

// decomp: Vienna RNA kümmert sich darum, übergibt identifiers
// z.B. multiloop Zerlegung. 

// wie entferne ich mich von Struktur A oder B 

// splitten der Teilintervalle, wir entfernen uns min. 1 bp von Struktur A

// aus 2D fold Paper... 

// erlaubt beliebig viele Strukturen A oder B, wir brauchen nur 2
// performanter: Cleanup, auf 2 Strukturen runterbrechen...


{
  short* pt1 = {120, 0, 0, 0, 0, 0, 0, 0, 0, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 0, 0, 0, 0, 0, 0, 0, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 0, 0, 0, 0, 0, 0, 80, 79, 78, 77, 76, 0, 75, 74, 0, 72, 70, 69, 68, 67, 66, 0, 65, 64, 0, 0, 0, 0, 59, 58, 56, 55, 54, 53, 52, 0, 51, 0, 49, 48, 46, 45, 44, 43, 42, 0, 0, 95, 94, 93, 92, 0, 0, 0, 0, 0, 86, 85, 84, 83, 0, 0, 0, 0, 0, 0, 0, 119, 118, 117, 115, 114, 113, 0, 0, 0, 0, 108, 107, 106, 0, 105, 104, 103, 0};
  short* pt2 = {120, 0, 0, 0, 0, 0, 0, 0, 0, 35, 34, 33, 32, 31, 30, 28, 27, 26, 25, 0, 0, 0, 0, 0, 0, 18, 17, 16, 15, 0, 14, 13, 12, 11, 10, 9, 86, 85, 84, 0, 0, 0, 80, 79, 0, 78, 77, 76, 75, 74, 0, 72, 71, 69, 68, 67, 66, 0, 65, 64, 0, 0, 0, 0, 59, 58, 56, 55, 54, 53, 0, 52, 51, 0, 49, 48, 47, 46, 45, 43, 42, 0, 0, 0, 38, 37, 36, 0, 120, 119, 118, 116, 115, 0, 114, 113, 112, 111, 110, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 98, 97, 96, 95, 94, 92, 91, 0, 90, 89, 88};

  size_t        *dist;
  unsigned int  r, numberOfRefs, **referenceBPs;
  int           ij, kl, *base_dx, *idx;
  short         **references_pt;

  references_pt = d->ref_pts;
  numberOfRefs  = d->ref_num;
  referenceBPs  = d->ref_bps;
  idx           = d->idx;

  ij  = idx[i] - j;
  kl  = idx[k] - l;


  // base_dx?
  // referenceBPs?

  // idx?

  base_dx = (int *)vrna_alloc(sizeof(int) * (numberOfRefs + 1));
  dist    = (size_t *)vrna_alloc(sizeof(size_t) * (numberOfRefs + 1));

  for (r = 0; r < numberOfRefs; r++)
    base_dx[r] = ((int)references_pt[r][i] != j) ? 1 : -1;

  printf("ij: %d kl: %d / %d %d\n", ij, kl, idx[i], j);
  printf("base_dx: %d %d\n", base_dx[0], base_dx[0]);
  printf("ref0: %d %d %d %d %d\n", references_pt[0][0], references_pt[0][1], references_pt[0][2], references_pt[0][3], references_pt[0][4]);
  
  switch (decomp) {
    /* cases where we actually introduce a base pair */
    case VRNA_DECOMP_PAIR_HP:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = base_dx[r] +
                  referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_PAIR_IL:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = base_dx[r] +
                  referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_PAIR_ML:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = base_dx[r] +
                  referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    /* cases where we split a segment into one or two subsegments */

    case VRNA_DECOMP_ML_STEM:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_ML_ML:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_ML_ML_ML:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_ML_UP:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_ML_ML_STEM:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_ML_COAXIAL: /* (i,j) stacks onto (k,l), lets ignore this case for now */
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = 0;
      break;

    case VRNA_DECOMP_EXT_EXT:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_EXT_UP:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_EXT_STEM:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_EXT_EXT_STEM: /* fall through */
    case VRNA_DECOMP_EXT_STEM_EXT:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - (j - 1)];
      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      for (r = 0; r < numberOfRefs; r++) {
        dist[r] = referenceBPs[r][ij] -
                   referenceBPs[r][kl];
        if (k > i)
          dist[r] -= referenceBPs[r][idx[i] - (k - 1)];

        if (l < j)
          dist[r] -= referenceBPs[r][idx[l + 1] - j];
      }
      break;

    default:
      fprintf(stderr, "default sc\n");
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = 0;
      break;
  }

  free(base_dx);

  printf("compute_distancies: i:%d k:%d l:%d j:%d / %c / dist0: %d / dist1: %d\n", i, k, l, j, decomp, dist[0], dist[1]);

  return dist;
}

// Energie evaluation, wir an RNAfold weiterreichen
// das sind die eigentlichen soft constraint FUnktinoen, die 
// RNAfold aufruft, wenn die Berechnungen im RNAfold stattfinden

// Distanzen bestimmen -> anhand den Distanzen eine Energieverzerrung
// Guiding potential / penalties
// Distanzvektor, dann 

// Datenstruktur d, hat funktion F

// f: Function pointer, irgendwelche energien zu generieren.
// das ist der Punkt, wo man unterschiedliche Varianten einsetzt
// Penalties: zu weit weg oder zu nah dran,
// deswegen abgekapselt als separate Funktion.

// 2 Funktionen: MFE oder andere Variante e^-kT
// ansonsten Funktionen identisch

FLT_OR_DBL
sc_exp_f_dist_class(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char decomp,
                    void          *data)
{
  size_t          *distancies;
  double          kT;
  FLT_OR_DBL      result;
  sc_dist_class_t *d;

  kT  = ((sc_dist_class_t *)data)->kT;
  d   = ((sc_dist_class_t *)data);

  distancies  = compute_distancies(i, j, k, l, decomp, d);
  result      = (FLT_OR_DBL)exp(-10. * (d->f(i, j, k, l, decomp, distancies, d)) / kT);
  // das ist der Unterschied zur anderen Funktion


  free(distancies);
  return result;
}

// unused?
int
sc_f_dist_class(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char decomp,
                void          *data)
{
  size_t          *distancies;
  int             result;
  sc_dist_class_t *d;

  d       = ((sc_dist_class_t *)data);
  result  = 0;

  distancies  = compute_distancies(i, j, k, l, decomp, d);

  // function pointer invocation
  result      = (int)d->f(i, j, k, l, decomp, distancies, d);

  free(distancies);

  return result;
}


// initializer der Datenstrukturen

sc_dist_class_t *
sc_dist_class_init(vrna_fold_compound_t *fc)
{
  sc_dist_class_t *d;

  d = (sc_dist_class_t *)vrna_alloc(sizeof(sc_dist_class_t));

  d->idx  = fc->iindx;
  d->kT   = fc->exp_params->kT;

  d->ref_num    = 0;
  d->references = NULL;
  d->ref_bps    = NULL;
  d->ref_pts    = NULL;

  d->f      = NULL;
  d->f_data = NULL;
  d->f_free = NULL;

  return d;
}

// destroy, free 

void
sc_dist_class_destroy(void *data)
{
  int             i;
  sc_dist_class_t *d = (sc_dist_class_t *)data;

  for (i = 0; i < d->ref_num; i++) {
    free(d->references[i]);
    free(d->ref_pts[i]);
    free(d->ref_bps[i]);
  }

  free(d->references);
  free(d->ref_bps);
  free(d->ref_pts);

  if (d->f_free)
    d->f_free(d->f_data);

  free(d);
}


// Funktionswrapper, für jede Struktur

void
sc_dist_class_add_ref(sc_dist_class_t *d,
                      const char      *ref_struct)
{
  int n = d->ref_num;

  d->references     = (char **)vrna_realloc(d->references, sizeof(char *) * (n + 1));
  d->references[n]  = strdup(ref_struct);
  d->ref_pts        = (short **)vrna_realloc(d->ref_pts, sizeof(short *) * (n + 1));
  d->ref_pts[n]     = vrna_ptable(ref_struct);
  d->ref_bps        = (unsigned int **)vrna_realloc(d->ref_bps, sizeof(unsigned int *) * (n + 1));
  d->ref_bps[n]     = vrna_refBPcnt_matrix(d->ref_pts[n], TURN);

  d->ref_num++;
}
