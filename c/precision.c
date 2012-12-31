/**************************************************************************
#
#  This software/database is "United States Government Work" under the terms of
#  the United States Copyright Act.  It was written as part of the authors'
#  official duties for the United States Government and thus cannot be
#  copyrighted.  This software/database is freely available to the public for
#  use without a copyright notice.  Restrictions cannot be placed on its present
#  or future use.
# 
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and data, the National Human Genome Research
#  Institute (NHGRI) and the U.S. Government does not and cannot warrant the
#  performance or results that may be obtained by using this software or data.
#  NHGRI and the U.S.  Government disclaims all warranties as to performance,
#  merchantability or fitness for any particular purpose.
# 
#  In any work or product derived from this material, proper attribution of the
#  authors as the source of the software or data should be made, using "NHGRI
#  Genome Technology Branch" as the citation.
#
**************************************************************************/

#include "bard.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

extern Params *parameters;

void small_double_check(double *minexparg)
{
    int checkno;
    double expval, logval;

    for (checkno = 0; checkno > -1000; checkno--) {
        expval = exp(1.0*checkno);
        logval = log(expval);

        /* fprintf(stderr, "%d %.10g %.10g\n", checkno, expval, logval); */
        if (logval != (1.0*checkno)) {
            *minexparg = 1.0*(checkno + 8); /* Adding 8 just to be safe */
            return;
        }
    }
}
