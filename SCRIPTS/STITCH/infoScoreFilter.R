args = commandArgs(trailingOnly=TRUE)
    if (! length(args)==3) {
        stop("supply arguments, 1-Input_RDataDir 2-info_score_filter 3-output_file_name", call.=FALSE)
} else {

input_EM_RData<-args[1]
infoScoreFilter<-args[2]
outputFileName<-args[3]

load(input_EM_RData)

filter<-(info<=infoScoreFilter)
excludeList<-paste0(pos[filter,"CHR"], ":", pos[filter,"POS"])

write.table(excludeList, file=outputFileName, quote=FALSE, row.names=FALSE, col.names=FALSE)

}
