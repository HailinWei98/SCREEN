library(scMAGeCK)

#get parameter
args<- commandArgs(T)
rds_dir<- args[1]
sg_dir<- args[2]
NTC<- args[3]
save_path<- args[4]

lr<- scmageck_lr(BARCODE = sg_dir,RDS = rds_dir, NEGCTRL = NTC,SAVEPATH = save_path)