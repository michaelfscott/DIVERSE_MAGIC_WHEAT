args = commandArgs(trailingOnly=TRUE)

sampleNames_path=args[1]
pos_1_path=args[2]
#pos_1_path="../../5STITCH/STITCH_input_old/chr1D_1:1-247726593/RData/pos.1D:1-247726593.1.247726593.RData"
pos_2_path=args[3]
#pos_2_path="../../5STITCH/STITCH_input_old/chr1D_2:1-247726593/RData/pos.1D:247726594-495453186.1.247726593.RData"
sampleReads_1_path=args[4]
#sampleReads_1_prefix="../../5STITCH/STITCH_input_old/chr1D_1:1-247726593/input/sample.1.input.1D:1-247726593.1.247726593.RData"
sampleReads_2_path=args[5]
#sampleReads_2_path="../../5STITCH/STITCH_input_old/chr1D_2:1-247726593/input/sample.1.input.1D:247726594-495453186.1.247726593.RData"
chromosome=args[6]
#chromosome="1D"
output_prefix=args[7]
#output_prefix="../../5STITCH/STITCH_input_old/chr1D"

source("../0Rfunctions/functions.r")

sampleReads_1_prefix=word(sampleReads_1_path, 1, sep="\\/input\\/")
sampleReads_1_suffix=word(sampleReads_1_path, 2, sep="\\.input\\.")
sampleReads_2_prefix=word(sampleReads_2_path, 1, sep="\\/input\\/")
sampleReads_2_suffix=word(sampleReads_2_path, 2, sep="\\.input\\.")

load(sampleNames_path)
sampleNumber=length(sampleNames)
save(sampleNames, file=paste0(output_prefix,"/RData/sampleNames.",chromosome,".RData"))

#load both pos files
load(pos_1_path)
pos_1<-pos
load(pos_2_path)
pos_2<-pos
#merge and save
pos<-rbind(update_pos(pos_1), update_pos(pos_2))
save(pos, file=paste0(output_prefix,"/RData/pos.1D.RData"))

for(sample in seq(1:sampleNumber)){
#load both sampleReads files
load(paste0(sampleReads_1_prefix, "/input/sample.", sample, ".input.", sampleReads_1_suffix))
sampleReads_1<-sampleReads
load(paste0(sampleReads_2_prefix, "/input/sample.", sample, ".input.", sampleReads_2_suffix))
sampleReads_2<-sampleReads

#add the length of pos_1 to positions in pos_2
add=nrow(pos_1)
new_sampleReads_2<-lapply(sampleReads_2, function(x)(list(x[[1]], x[[2]]+add, x[[3]], x[[4]]+add)))
#append and save
sampleReads<-append(sampleReads_1, new_sampleReads_2)
save(sampleReads, file=paste0(output_prefix,"/input/sample.",sample,".input.",chromosome,".RData"))

}
