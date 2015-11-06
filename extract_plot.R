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

library(bamcheckr)
library(reshape2)

load_stat <- function(lanelet_file)
{
    bamcheck <- read_bamcheck(lanelet_file)
    lanelet <- sub("[.][^.]*$", "", basename(lanelet_file), perl=TRUE)
    
    bamcheck$data$SN$lanelet <- rep(lanelet,nrow(bamcheck$data$SN))
    extract <-dcast(bamcheck$data$SN, lanelet~variable)
    extract[,c(
    'filtered sequences:',
    'reads duplicated:',
    'reads mapped:',
    'reads paired:',
    'reads mapped and paired:',
    'reads properly paired:',
    'error rate:',
    'fwd.percent.insertions.above.baseline:',
    'fwd.percent.insertions.below.baseline:',
    'fwd.percent.deletions.above.baseline:',
    'fwd.percent.deletions.below.baseline:',
    'rev.percent.insertions.above.baseline:',
    'rev.percent.insertions.below.baseline:',
    'rev.percent.deletions.above.baseline:',
    'rev.percent.deletions.below.baseline:',
    'quality.dropoff.fwd.high.iqr.start.read.cycle:',
    'quality.dropoff.fwd.high.iqr.end.read.cycle:',
    'quality.dropoff.fwd.high.iqr.max.contiguous.read.cycles:',
    'quality.dropoff.fwd.mean.runmed.decline.start.read.cycle:',
    'quality.dropoff.fwd.mean.runmed.decline.end.read.cycle:',
    'quality.dropoff.fwd.mean.runmed.decline.max.contiguous.read.cycles:',
    'quality.dropoff.fwd.mean.runmed.decline.high.value:',
    'quality.dropoff.fwd.mean.runmed.decline.low.value:',
    'quality.dropoff.rev.high.iqr.start.read.cycle:',
    'quality.dropoff.rev.high.iqr.end.read.cycle:',
    'quality.dropoff.rev.high.iqr.max.contiguous.read.cycles:',
    'quality.dropoff.rev.mean.runmed.decline.start.read.cycle:',
    'quality.dropoff.rev.mean.runmed.decline.end.read.cycle:',
    'quality.dropoff.rev.mean.runmed.decline.max.contiguous.read.cycles:',
    'quality.dropoff.rev.mean.runmed.decline.high.value:',
    'quality.dropoff.rev.mean.runmed.decline.low.value:',
    'quality.dropoff.high.iqr.threshold:',
    'quality.dropoff.runmed.k:',
    'quality.dropoff.ignore.edge.cycles:',
    'A.percent.mean.above.baseline:',
    'C.percent.mean.above.baseline:',
    'G.percent.mean.above.baseline:',
    'T.percent.mean.above.baseline:',
    'A.percent.mean.below.baseline:',
    'C.percent.mean.below.baseline:',
    'G.percent.mean.below.baseline:',
    'T.percent.mean.below.baseline:',
    'A.percent.max.above.baseline:',
    'C.percent.max.above.baseline:',
    'G.percent.max.above.baseline:',
    'T.percent.max.above.baseline:',
    'A.percent.max.below.baseline:',
    'C.percent.max.below.baseline:',
    'G.percent.max.below.baseline:',
    'T.percent.max.below.baseline:',
    'A.percent.max.baseline.deviation:',
    'C.percent.max.baseline.deviation:',
    'G.percent.max.baseline.deviation:',
    'T.percent.max.baseline.deviation:',
    'A.percent.total.mean.baseline.deviation:',
    'C.percent.total.mean.baseline.deviation:',
    'G.percent.total.mean.baseline.deviation:',
    'T.percent.total.mean.baseline.deviation:'
    )]
}

input_list <- readLines('../pomak.list')

dat <- mclapply(input_list, load_stat)


dat$error_rate <-dat$'error rate:' *100
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
