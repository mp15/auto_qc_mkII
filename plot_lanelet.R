# plot.R
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
# Created by Martin Pollard on 21/10/2015.
#

#conf <- read.table('autoqc.conf', header=TRUE)
#conf <- dcast(conf, row~value, fill=0)



dat <- read.table('output.list', header=TRUE, comment.char="~")
dat$error_rate <-dat$error_rate *100
dat$dup_rate <- dat$dup_reads / dat$filtered_reads * 100
dat$map_rate <- dat$mapped_reads / dat$filtered_reads * 100
dat$properpair_rate <- dat$properly_paired_reads / dat$filtered_reads * 100
dat$indel_ratio <-dat$ins_total / dat$del_total
dat$map_minus_dup_coverage <- ( (dat$mapped_reads - dat$dup_reads) * 151 ) / 3000000000
dat$auto_qc_error_rate <- ifelse( dat$error_rate < 2,ifelse(dat$error_rate < 1,"PASS","WARNING"),"FAIL")
dat$auto_qc_dup_rate <- ifelse( dat$dup_rate < 20,ifelse(dat$dup_rate < 15,"PASS","WARNING"),"FAIL")
dat$auto_qc_map_rate <- ifelse( dat$map_rate > 90,ifelse(dat$map_rate > 95,"PASS","WARNING"),"FAIL")
dat$auto_qc_properpair_rate <- ifelse( dat$properpair_rate > 80,ifelse(dat$properpair_rate > 90,"PASS","WARNING"),"FAIL")
dat$auto_qc_indel_ratio <- ifelse( ((dat$indel_ratio > 1.105) + (dat$indel_ratio < 0.450)) > 0, "FAIL", ifelse(((dat$indel_ratio > 0.825) + (dat$indel_ratio < 0.675)) > 0, "WARNING", "PASS"))
dat$auto_qc <- ifelse(
(
(dat$auto_qc_error_rate == "FAIL") + 
(dat$auto_qc_dup_rate == "FAIL") +
(dat$auto_qc_map_rate == "FAIL") +
(dat$auto_qc_properpair_rate == "FAIL") +
(dat$auto_qc_indel_ratio == "FAIL") 
)
> 0, "FAIL", ifelse(
(
(dat$auto_qc_error_rate == "WARNING") +
(dat$auto_qc_dup_rate == "WARNING") +
(dat$auto_qc_map_rate == "WARNING") +
(dat$auto_qc_properpair_rate == "WARNING") +
(dat$auto_qc_indel_ratio == "WARNING") 
)
 > 0, "WARNING","PASS"))
summary(dat)

summary(as.factor(dat$auto_qc))
summary(as.factor(dat$auto_qc_error_rate))
summary(as.factor(dat$auto_qc_dup_rate))
summary(as.factor(dat$auto_qc_map_rate))
summary(as.factor(dat$auto_qc_properpair_rate))
summary(as.factor(dat$auto_qc_indel_ratio))
