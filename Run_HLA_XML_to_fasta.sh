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

source activate minionvironment

cd Parse_IMGT_HLA

python HLA_XML_to_fasta.py /home/minion/MinIONData/HLA_XML/hla.xml/hla.xml /home/ben/MUMCScripts/HlaAlleles/OutputDirectory/

source deactivate
