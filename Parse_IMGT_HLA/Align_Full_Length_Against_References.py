# This file is part of HLA-allele-analysis.
#
# HLA-allele-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HLA-allele-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with HLA-allele-analysis. If not, see <http://www.gnu.org/licenses/>.

import sys
import os

if __name__=='__main__':
    try:
        # Parse Args. Get Input files Get output Directory.


        #if not os.path.isdir(outputDirectory):
        #    os.mkdir(outputDirectory)


        print('Aligning Full-length sequences against reference')

        # Parse Reference Sequences

        # Parse all full-length sequences

        # For each reference sequence
            # For each full-length Sequence
                # If the allele locus and group match add to list

            # Print Single Reference Sequence

            # Print matching full-length alleles

            # Align!

            # Do something else?

        print('Done.  Ben did a great job.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print('Unexpected problem during execution:')
        print(sys.exc_info()[1])
        raise
