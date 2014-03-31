/* map.c: calculate distances of objects to unis
   This version: Sept 26, 2006
   Author: Ron Wehrens

	Modified Aug 16, 2013
	Alan P Boyle
*/

#include <R.h>

void mapKohonen(double *data, double *codes,
		Sint *pncodes, Sint *pnd, Sint *pp, double *dists)
{
  int ncodes = *pncodes, nd = *pnd, p = *pp;
  int i, j, cd, distcounter, NAcounter;
  double tmp;

	//Added variables - APB
	double dists_tmp;
  double dists_min;
	double min_unit;
	double MAX_VAL = 100000000.0;

  /* i is a counter over objects in data, cd  is a counter over SOM
     units, and j is a counter over variables. */
  distcounter = -1;
  for (i = 0; i < nd; i++) {
		dists_min = MAX_VAL;  // APB
		min_unit = 0.0;	//APB

    for (cd = 0; cd < ncodes; cd++) {

      NAcounter = 0;

			// APB
			// Make this now a pre-processed entry per data item instead of all possible combinations
			//dists[++distcounter] = 0.0;
			dists_tmp = 0.0;

      for (j = 0; j < p; j++) {
				if (!ISNAN(data[i + j*nd])) {
	  			tmp = data[i + j*nd] - codes[cd + j*ncodes];
	  			dists_tmp += tmp * tmp;
				} else {
	  			NAcounter++;
				}
      }

      if (NAcounter == p) {
				dists_tmp = NA_REAL;
      } else { /* correct distance for lower number of non-NA variables */
				if (NAcounter > 0) dists_tmp *= p/(p - NAcounter);
      }

			if(!ISNAN(dists_tmp)) {
				if(dists_tmp < dists_min) {
					dists_min = dists_tmp;
					min_unit = (double)cd;
				}
			}

    }

		if(dists_min == MAX_VAL) {
			dists[++distcounter] = NA_REAL;
			dists[++distcounter] = NA_REAL;
		} else {
			dists[++distcounter] = dists_min;
			dists[++distcounter] = min_unit;
		}
  }
}
