#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/boltzmann_sampling.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/mm.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>

#include "dist_class_sc.h"

typedef struct {
    size_t*       maxDistances;
    const double* distortions;
} kl_sc_dat_t;

// unnötig viel Code, dinge aus RNAexplorer zusammenkopiert
//

static size_t max_bp_dist(const char* sequence, const char* structure)
{
    short*        pt;
    size_t        max_dist, i, n, bps;
    unsigned int* mm1;
    int *         iindx, idx_1n;

    n      = strlen(sequence);
    pt     = vrna_ptable(structure);
    mm1    = maximumMatchingConstraint(sequence, pt);
    iindx  = vrna_idx_row_wise(n);
    idx_1n = iindx[1] - n;

    /* get number of bp in structure */
    for (bps = 0, i = 1; i < n; i++)
        if (pt[i] > i) bps++;

    /* compute maximum distance to this structure */
    max_dist = bps + (size_t)mm1[idx_1n];

    free(pt);
    free(mm1);
    free(iindx);

    return max_dist;
}

// das ist der springende Punkt,
// die wird an Functino pointer angebunden

// d -> f

// über alle Referenzstrukturen übergehen
// Distortion default, wie verändern wir die Energien, wie weit wir von A und B entfernt sind
// Problem: Wie will man das überhaupt machen?
// Hauptsächlich entlant des direkten Pfaden samplen,
// sobald wir uns entferen, haben wir eine Art Penalty

// was wichtig wäre: Bp distance A und B, das
// wir hier noch nirgends gemacht
// und wenn ich den decomposition step A und B,
// daraus etwas herauskristallisieren

// wir haben aber noch keine Struktur C, wir haben uns nur so und so viel von A oder B entfernt

// solange decomposition passen, und noch keine Struktur da ist,

// induzierte Distanzen A und zu B aufaddieren, dann kommen wir irgendwann zu dem
// Punkt, das wir auf der Mitte liegen, idealierweise die in der Mitte
// während wir soft constraints generieren
// erst sobald wir größere Distanzen induzieren, dann können wir dagegensteuern,
// alles was noch weiter wegliegt wir penalized.
// Überlegen: Wie können wir am besten eine Funktion aufschreiben.

// d(s, A) + d(s, B) - d(A,B) <= theta (theta als upper bound, detour distance upper bound)
// je weiter man von der idealen Gerade wegsind, umso mehr penalties

// einfache Energiefunktion
// [ d(s, A) + d(s, B) - d(A,b) + max_detour ] * penalty
// 3 Funktionen, max detour wie weit sind wir maximal weg, mit einem pena

// Teilkomponenten während der Berechnung (Grafik)

// Funktion müsste auch bei Komponentenberechnung stimmen k->l
// wir betrachten immer nur ein Intervall, betrachten dann müsste die Formel noch funktionieren.

// d wird schon berechnet

static double distortion_default(int i, int j, int k, int l, unsigned char decomp, size_t* distance,
                                 sc_dist_class_t* d)

// Komponentenzerlegung k->l wie auf Grafik
// Funktion bekommt Distanzen zu Struktur A und B, die
// induziert werden in der Zerlegung

// d([i,j], A[i,j])
// bzw.
// d([i,j], B[i,j])

// Was ist die Distanz im Intervall i->j die ich erzeuge zur Struktur A und B,
// wenn ich die Zerlegung mache
// 2 Komponenten in der Formel
// das einzige was ich noch nicht weiß?

// was ist d(A[i,j], B[i,j])

// wir brauchen die 3 Informationen, Dreiecksungleichung
// wenn wir das haben: [ d(s, A) + d(s, B) - d(A,b) + max_detour ] * penalty
// Problem: siehe Grafik mit Kreisen: Wenn wir nicht auf Strukturebene sondern Komponentenebene
// wir uns was anschauen, dadurch schauen wir uns alles an in einem Blob um A und B rundherum.

// Variante im RNA explorer:
// B -> Kreise, je weiter wir uns an A annähern - viele Iterationen
// ziemlich aufwändig, viele Partition Functions ausrechnen
// Schrittweises annähern von A nach B
// in jedem Schritt

// extended options -e NBRSEV
// attractions
//  -e, --extended_opt=STRING     Some extended options:
//                                   N    normal distortion (no shift)
//                                   B    alter both potentials at once
//                                   R    relax potential instead of increasing it
//                                   S    shift potential to other structure
//                                   F    shift to first structure
//                                   V    verbose

// String ist S für shift, shift von A -> B Schrittweise, nicht gut beschrieben
// weils eben nicht gut funktioniert hat
// iterativer Shift ist viel zu aufwendig - einmal max detour definieren, maximale  Entfernung
// straight line, penalty Faktor - herumspielen. muss wahrscheinlich ein großer Energiebetrag sein,
// das wir Dinge außerhalb der Ellipse stark abschwächen

{
    size_t        num_refs;
    double        result;
    const double* distortions;

    distortions = ((kl_sc_dat_t*)d->f_data)->distortions;
    num_refs    = d->ref_num; // =2
    

    int pt1[] = {120, 0, 0, 0, 0, 0, 0, 0, 0, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 0, 0, 0, 0, 0, 0, 0, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 0, 0, 0, 0, 0, 0, 80, 79, 78, 77, 76, 0, 75, 74, 0, 72, 70, 69, 68, 67, 66, 0, 65, 64, 0, 0, 0, 0, 59, 58, 56, 55, 54, 53, 52, 0, 51, 0, 49, 48, 46, 45, 44, 43, 42, 0, 0, 95, 94, 93, 92, 0, 0, 0, 0, 0, 86, 85, 84, 83, 0, 0, 0, 0, 0, 0, 0, 119, 118, 117, 115, 114, 113, 0, 0, 0, 0, 108, 107, 106, 0, 105, 104, 103, 0};
    int pt2[] = {120, 0, 0, 0, 0, 0, 0, 0, 0, 35, 34, 33, 32, 31, 30, 28, 27, 26, 25, 0, 0, 0, 0, 0, 0, 18, 17, 16, 15, 0, 14, 13, 12, 11, 10, 9, 86, 85, 84, 0, 0, 0, 80, 79, 0, 78, 77, 76, 75, 74, 0, 72, 71, 69, 68, 67, 66, 0, 65, 64, 0, 0, 0, 0, 59, 58, 56, 55, 54, 53, 0, 52, 51, 0, 49, 48, 47, 46, 45, 43, 42, 0, 0, 0, 38, 37, 36, 0, 120, 119, 118, 116, 115, 0, 114, 113, 112, 111, 110, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 98, 97, 96, 95, 94, 92, 91, 0, 90, 89, 88};


    // printf("i: %d j: %d / %d %d\n", i,j, k,l);



    result = 0.;

    // vorfaktoren etc. adjust
    for (size_t r = 0; r < num_refs; r++) {
        result += distortions[r] * (double)distance[r];
    }

    // printf("%d %d %d %d / %d / d0: %d / d1: %d\n", i, j, k, l, num_refs, distance[0], distance[1]);

    result = result * 100.;
    
    
    if ((int)pt1[i] == j || (int)pt2[i] == j) {   
        // printf("i: %d j: %d / %d %d\n", i,j, k,l);
        // printf("%d %d %d %d / %d / d0: %d / d1: %d / result: %f\n", i, j, k, l, num_refs, distance[0], distance[1], result);

        result = result/50;
    }
    
    return result;
}

static double distortion_repel(int i, int j, int k, int l, unsigned char decomp, size_t* distance,
                               sc_dist_class_t* d)
{
    double        result;
    const double* distortions;
    size_t        num_refs, *maxDist;

    distortions = ((kl_sc_dat_t*)d->f_data)->distortions;
    maxDist     = ((kl_sc_dat_t*)d->f_data)->maxDistances;
    num_refs    = d->ref_num;

    result = 0.;
    for (int r = 0; r < num_refs; r++)
        result += distortions[r] * ((double)maxDist[r] - distance[r]);

    result = result * 100.;
    return result;
}

static void kl_sc_dat_destroy(void* data)
{
    kl_sc_dat_t* d = (kl_sc_dat_t*)data;

    free(d->maxDistances);
    free(d);
}

static sc_dist_class_t* init_sc(vrna_fold_compound_t* fc, const char** referenceStructures,
                                size_t num_ref, const double* distortions)
{
    char*            s;
    sc_dist_class_t* d;
    kl_sc_dat_t*     data;

    s = fc->sequence;
    d = sc_dist_class_init(fc);

    for (size_t i = 0; i < num_ref; i++) sc_dist_class_add_ref(d, referenceStructures[i]);

    d->f = distortion_default;  // set the function pointer
    // d->f = distortion_repel; // set the function pointer

    data               = (kl_sc_dat_t*)vrna_alloc(sizeof(kl_sc_dat_t));
    data->distortions  = distortions;
    data->maxDistances = vrna_alloc(sizeof(size_t*) * num_ref);
    for (size_t i = 0; i < num_ref; i++)
        data->maxDistances[i] = max_bp_dist(s, referenceStructures[i]);

    d->f_data = (void*)data;
    d->f_free = &kl_sc_dat_destroy;

    return d;
}

// Stochastic backtracking

// main function: Random String, MFE für random string
// ausrechnen der Partition function für stochastic backtracking
// dann 2 Sample Strukturen

// Zeile 189: Strukturen A und B

// Penalty vom 0.1 kcal/mol

// dann init_sc aus dist_call_sc
// RNAfold wrapper
// erzeugt die Distortion Datenstruktur
// dann kommt das Binden der Soft constraints nach RNAfold

int main()
{
    /* initialize random number generator */
    vrna_init_rand();

    /* Generate a random sequence of 50 nucleotides */
    char* seq = vrna_random_string(120, "ACGU");

    seq =
        "UAAGGAAACUGAUGAGGGCAAAGUCUCUUCAUUGGCGCAAAACGGGGUAGAUGUCGGCCUGCAGGUGAUGAUGUAUCCGUCCGGCGAUCA"
        "ACGCCUAAUUCGCAUCUCUUUAGGGUAUGA";

    /* Create a fold compound for the sequence */
    vrna_fold_compound_t* fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

    /* allocate memory for MFE structure (length + 1) */
    char* structure = (char*)vrna_alloc(sizeof(char) * 121);

    /* predict Minmum Free Energy and corresponding secondary structure */
    double mfe = (double)vrna_mfe(fc, structure);

    /* print seqeunce, structure and MFE */
    printf("%s\n%s [ %6.2f ] <-- MFE\n", seq, structure, mfe);

    fc->params->model_details.compute_bpp = 0; /* deactivate base pair probability computation */
    fc->params->model_details.uniq_ML     = 1; /* required for Boltzmann sampling */

    /* Rescale Boltzmann factors to avoid overflows */
    vrna_exp_params_rescale(fc, &mfe);

    /* compute partition function  */
    (void)vrna_pf(fc, NULL);

    /* generate two structure samples */
    char** samples = vrna_pbacktrack_num(fc, 2, VRNA_PBACKTRACK_NON_REDUNDANT);

    samples[0] =
        "........((((((((((.......))))))))))......(((((.((.((((((.((....))))))).).)))))))..((((...."
        ".)))).......((((((....))).))).";
    samples[1] =
        "........((((((((((......)))).))))))(((...((.(((((.((((((.((....)))))).)).)))))))...))).((("
        "((.(((((...........))))))).)))";

    printf("--------------------------------\n");

    /* print both structures */
    printf("%s <-- S1\n%s <-- S2\n", samples[0], samples[1]);

    printf("--------------------------------\n");


    /* now, perform distortion of energies according to both references */
    double distortions[2];

    distortions[0] = 0.5; // original value: 0.1
    distortions[1] = 0.5;

    sc_dist_class_t* data = init_sc(fc, (const char**)samples, 2, (const double*)distortions);

    //

    vrna_sc_init(fc); /*  to remove old soft constraints */
    vrna_sc_add_data(fc, (void*)data, &sc_dist_class_destroy);
    // Datenstrukturen, wird vom dist classes sc gebraucht

    // soft constraint function point
    //

    // add the custom function pointer / soft constraint
    vrna_sc_add_f(fc, &sc_exp_f_dist_class);
    vrna_sc_add_exp_f(fc, &sc_exp_f_dist_class);

    // jetzt kümmert sich RNAfold um das, callback Mechanismus

    // geänderte Energielandschaft, jetzt berechnen wir uns 50 Strukturen
    // ohne Redundanz
    // Stochastic backtracing ohne Strukturen zurücklegen,
    // Mist ausgeben

    /* recompute MFE and partition function under soft constraints */
    mfe = (double)vrna_mfe(fc, NULL);
    vrna_exp_params_rescale(fc, &mfe);
    (void)vrna_pf(fc, NULL);

    /* sample 50 structures from distorted landscape */
    // char** d_samples = vrna_pbacktrack_num(fc, 50, VRNA_PBACKTRACK_NON_REDUNDANT);
    char** d_samples = vrna_pbacktrack_num(fc, 10, VRNA_PBACKTRACK_NON_REDUNDANT);

    for (size_t i = 0; d_samples[i]; i++) {
        printf("%s - d_sample %u\n", d_samples[i], i + 1);
        free(d_samples[i]);
    }
    free(d_samples);

    return 0;

    /* cleanup */
    for (size_t i = 0; i < 2; i++) free(samples[i]);
    free(samples);

    /* cleanup memory */
    free(seq);
    free(structure);
    vrna_fold_compound_free(fc);

    return 0;
}

// non redundant Boltzmann sampling
// Energien gewichtet
// Stochastic backtracking, n^3.
// partition function berechnung braucht läner als MFE
