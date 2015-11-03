#!/bin/sh

# stats_to_tsv.sh
#  
# Copyright 2015 Genome Research Limited
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Created by Martin Pollard on 21/10/2015.
#

(echo -e "library\tfiltered_reads"; grep -H 'SN'$'\t''sequences:' *.stats | sed -e 's/:/'$'\t''/g' -e 's/\.stats//g' | cut -f 1,5) > total.list
(echo -e "library\tdup_reads"; grep -H 'reads duplicated:' *.stats | sed -e 's/:/'$'\t''/g' -e 's/\.stats//g' | cut -f 1,5) > dup.list
(echo -e "library\tmapped_reads"; grep -H 'reads mapped:' *.stats | sed -e 's/:/'$'\t''/g' -e 's/\.stats//g' | cut -f 1,5) > mapped.list
(echo -e "library\tmappedpaired_reads"; grep -H 'reads mapped and paired:' *.stats | sed -e 's/:/'$'\t''/g' -e 's/\.stats//g' | cut -f 1,5) > mappedpaired.list
(echo -e "library\tproperly_paired_reads"; grep -H 'reads properly paired:' *.stats | sed -e 's/:/'$'\t''/g' -e 's/\.stats//g' | cut -f 1,5) > properlypaired.list
(echo -e "library\tpaired_reads"; grep -H 'reads paired:' *.stats | sed -e 's/:/'$'\t''/g' -e 's/\.stats//g' | cut -f 1,5) > paired.list
(echo -e "library\terror_rate"; grep -H 'error rate' *.stats | sed -e 's/:/'$'\t''/g' -e 's/\.stats//g' | cut -f 1,5) > error.list
(echo -ne "library\tins_total\tdel_total";grep -H '^ID' *.stats | sed -e 's/:/'$'\t''/g'  -e 's/\.stats//g' | awk 'BEGIN{library=""; OFS="\t"}{if ($1 != library){print library, INS, DEL; INS=0; DEL=0;library=$1} INS += $4; DEL+=$5}END{print library,INS, DEL}') >ins_del.list

join <(join total.list dup.list) <(join mapped.list mappedpaired.list) > tmp.list
join <(join tmp.list properlypaired.list) <(join paired.list error.list) > tmp2.list
join tmp2.list ins_del.list > output.list

Rscript --vanilla < plot.R
