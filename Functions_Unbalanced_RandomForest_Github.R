
################################ calculate variable importances

library(randomForest)
library(party)
library(varImp)
library(ranger)

FunMyVarImps <- function(InData){
	ImpsOut <- as.data.frame(matrix(nrow = ncol(InData)-1, ncol = 22, dimnames = list(names(InData)[-1], 
	c("RF_mean", "RF_sd", "Down_mean", "Down_sd", "DownAcc_mean", "DownAcc_sd", "WeightMin_mean", "WeightMin_sd", "WeightMax_mean", "WeightMax_sd", 
	"CFor_mean", "CFor_sd", "CForAUC_mean", "CForAUC_sd", 
	"Ranger_mean", "Ranger_sd", "RangWeightMin_mean", "RangWeightMin_sd", "RangWeightMax_mean", "RangWeightMax_sd", "Hell_mean", "Hell_sd"))))
	ImpsBase <- matrix(nrow = ncol(InData)-1, ncol = 10)
	InData$Y <- factor(InData$Y)
#################### calculate permutation importance for 10 repetitions of model with all variables
# Version RF
	Imps10 <- ImpsBase
	for (i in 1 : 10){																
		rFi <- randomForest(Y ~ ., data = InData, importance = T, ntree = 2000, mtry = round((ncol(InData)-1)/2))	
		Imps10[, i] <- randomForest::importance(rFi, scale = F, type = 1)
	}
# calculate mean permutation importance of 10 repetitions
	ImpsOut$RF_mean <- apply(Imps10, 1, mean)
# calculate sd perm imp of 10 rep							
	ImpsOut$RF_sd <- apply(Imps10, 1, sd)							
# Version RF-Down
	Imps10 <- ImpsBase
	NoMinClass <- length(InData$Y[InData$Y == 1])		
	for (i in 1 : 10){																
		rFi <- randomForest(Y~., data = InData, importance = T, ntree = 2000, mtry = round((ncol(InData)-1)/2), sampsize = c(NoMinClass, NoMinClass))	
		Imps10[, i] <- randomForest::importance(rFi, scale = F, type = 1)
	}
	ImpsOut$Down_mean <- apply(Imps10, 1, mean)						
	ImpsOut$Down_sd <- apply(Imps10,1, sd)
# Version RF-Down-Acc		
	Imps10 <- ImpsBase
	for (i in 1 : 10){																
		rFi <- randomForest(Y~., data = InData, importance = T, ntree = 2000, mtry = round((ncol(InData)-1)/2), sampsize = c(round(NoMinClass * 0.64), NoMinClass))	
		Imps10[, i] <- randomForest::importance(rFi, scale = F, type = 1)
	}
	ImpsOut$DownAcc_mean <- apply(Imps10, 1, mean)						
	ImpsOut$DownAcc_sd <- apply(Imps10,1, sd)
# Version RF-Priors-Min-Max
	Imps10 <- ImpsBase	
	NoMaxClass <- length(InData$Y[InData$Y == 0])
	MyWeightsMin <- c(NoMinClass/nrow(InData), NoMaxClass/nrow(InData))
	for (i in 1 : 10){																
		rFi <- randomForest(Y~., data = InData, importance = T, ntree = 2000, mtry = round((ncol(InData)-1)/2), classwt = MyWeightsMin)	
		Imps10[, i] <- randomForest::importance(rFi, scale = F, type = 1)
	}
	ImpsOut$WeightMin_mean <- apply(Imps10, 1, mean)						
	ImpsOut$WeightMin_sd <- apply(Imps10,1, sd)
# Version RF-Priors / RF-Priors-Min-Min
	Imps10 <- ImpsBase
	MyWeightsMax <- c(NoMaxClass/nrow(InData), NoMinClass/nrow(InData))
	for (i in 1 : 10){																
		rFi <- randomForest(Y~., data = InData, importance = T, ntree = 2000, mtry = round((ncol(InData)-1)/2), classwt = MyWeightsMax)	
		Imps10[, i] <- randomForest::importance(rFi, scale = F, type = 1)
	}
	ImpsOut$WeightMax_mean <- apply(Imps10, 1, mean)						
	ImpsOut$WeightMax_sd <- apply(Imps10,1 , sd)
# Versions CFor and CFor-AUC
	Imps10 <- ImpsBase
	ImpsAUC10 <- ImpsBase
	for (i in 1 : 10){											
		cFi <- cforest(Y~., data=InData, controls = cforest_unbiased(ntree=2000, mtry=round((ncol(InData)-1)/2)))	
 		Imps10[, i] <- party::varimp(cFi, conditional=F)
		ImpsAUC10[, i] <- varImp::varImpAUC(cFi, conditional=F)						
	}
	ImpsOut$CFor_mean <- apply(Imps10, 1, mean)						
	ImpsOut$CFor_sd <- apply(Imps10,1 , sd)
	ImpsOut$CForAUC_mean <- apply(ImpsAUC10, 1, mean)						
	ImpsOut$CForAUC_sd <- apply(ImpsAUC10,1, sd)
# Ranger without weights. Results not presented in the manuscript
	Imps10 <- ImpsBase
	for (i in 1 : 10){																
		rFi <- ranger(Y~., data = InData, importance = "permutation", num.trees = 2000, mtry = round((ncol(InData)-1)/2))	
		Imps10[, i] <- rFi$variable.importance
	}	
	ImpsOut$Ranger_mean <- apply(Imps10, 1, mean)						
	ImpsOut$Ranger_sd <- apply(Imps10,1, sd)
# Version: Rang-Weight-Min-Max
	Imps10 <- ImpsBase
	for (i in 1 : 10){																
		rFi <- ranger(Y~., data = InData, importance = "permutation", num.trees = 2000, mtry = round((ncol(InData)-1)/2), class.weights = MyWeightsMin)	
		Imps10[, i] <- rFi$variable.importance
	}	
	ImpsOut$RangWeightMin_mean <- apply(Imps10, 1, mean)						
	ImpsOut$RangWeightMin_sd <- apply(Imps10,1, sd)
#
# Version: Rang-Weight-Min-Min
	Imps10 <- ImpsBase
	for (i in 1 : 10){																
		rFi <- ranger(Y~., data = InData, importance = "permutation", num.trees = 2000, mtry = round((ncol(InData)-1)/2), class.weights = MyWeightsMax)	
		Imps10[, i] <- rFi$variable.importance
	}
	ImpsOut$RangWeightMax_mean <- apply(Imps10, 1, mean)						
	ImpsOut$RangWeightMax_sd <- apply(Imps10,1, sd)
# Version: Rang-Hellinger
	Imps10 <- ImpsBase
	for (i in 1 : 10){																
		rFi <- ranger(Y~., data = InData, importance = "permutation", num.trees = 2000, mtry = round((ncol(InData)-1)/2), splitrule = "hellinger")	
		Imps10[, i] <- rFi$variable.importance
	}
	ImpsOut$Hell_mean <- apply(Imps10, 1, mean)						
	ImpsOut$Hell_sd <- apply(Imps10,1, sd)
	return(ImpsOut)
}

############################################### variable selection

library(randomForest)
library(party)
library(varImp)
library(rpart)
library(ROCR)

# Versions CFor-AUC/OOB and CFor-AUC/AUC

VarImpcForest <- function(InData){									
	InData$Y <- factor(InData$Y)																
	Imps <- matrix(nrow = ncol(InData) - 1, ncol = 55, dimnames = list(names(InData)[-1], c(1 : 50, "mean", "sd", "VarNo", "CartPred", "ThreshPos")))
	for (i in 1 : 50){											
		cFi <- cforest(Y ~ ., data = InData, controls = cforest_unbiased(ntree = 2000, mtry = round((ncol(InData)-1)/2)))	
 		Imps[, i] <- varImp::varImpAUC(cFi, conditional = F)						
	}		
	Imps[, 51] <- apply(Imps[, 1 : 50], 1, mean)							
	Imps[, 52] <- apply(Imps[, 1 : 50], 1, sd)														
	sort <- order(Imps[, 51], decreasing = TRUE)									
	Imps <- Imps[sort, ]																			
	Imps[, 53] <- 1 : (length(InData)-1)																																				
	Imps[, 54] <- predict(rpart(Imps[, 52] ~ Imps[, 53]))					
	Imps[, 55] <- Imps[ , 51] - min(Imps[, 54])							
	Thresh <- ifelse(all(Imps[, 55] > 0), nrow(Imps), (which(Imps[, 55] < 0) - 1)[1]) 												 
	NestedData <- cbind(InData[, 1, drop = F], InData[ , -1][, sort])																										
	NestedData <- NestedData[, 1 : (Thresh + 1), drop = F]
	if(ncol(NestedData) <= 2){
		SelDataOOBs = NestedData 
		SelDataAUCs = NestedData 
		SelDataAUCPRs = NestedData
	}else{
# Version CFor-AUC/OOB 									
		OOBs <- matrix(ncol = 27, nrow = length(NestedData)-1, dimnames = list(paste0("NoOfVar", 1 : (length(NestedData)-1)), c(paste0("Rep", 1:25), "mean", "sd")))
# Version CFor-AUC/AUC
		AUCs <- OOBs					
# Additional version not presented in the manuscript
		AUCPRs <- OOBs
		for (j in 1 : (length(NestedData)-1)){								
			subs <- NestedData[, 1 : (j+1)]						
			for (k in 1:25){				
				rFk <- randomForest(Y ~ ., data = subs)
# Version CFor-AUC/OOB 	
				OOBs[j, k] <- rFk$err.rate[500, 1]	
# Version CFor-AUC/AUC
				MyPreds <- ROCR::prediction(rFk$votes[, 2], as.numeric(as.character(NestedData$Y)))
				AUCCalc <- ROCR::performance(MyPreds, "auc")
				AUCs[j, k] <- AUCCalc@y.values[[1]]
# Additional version not presented in the manuscript
				AUCPRCalc <- ROCR::performance(MyPreds, "aucpr")
				AUCPRs[j, k] <- AUCPRCalc@y.values[[1]]			
			}							
		}
# Version CFor-AUC/OOB	
		OOBs[, 26] <- apply(OOBs[, 1 : 25], 1, mean)							
		OOBs[, 27] <- apply(OOBs[, 1 : 25], 1, sd)			 
		OOBMinModel <- which.min(OOBs[, 26])							
		ThresModel <- which(OOBs[, 26] <= OOBs[OOBMinModel, 26] + OOBs[OOBMinModel, 27])  
		SelDataOOBs <- NestedData[, 1 : (ThresModel[1] + 1)]						
# Version CFor-AUC/AUC			
		AUCs[, 26] <- apply(AUCs[, 1 : 25], 1, mean)							
		AUCs[, 27] <- apply(AUCs[, 1 : 25], 1, sd)							 
		AUCMaxModel <- which.max(AUCs[, 26])							
		ThresModelAUC <- which(AUCs[, 26] >= AUCs[AUCMaxModel, 26] - AUCs[AUCMaxModel, 27])  
		SelDataAUCs <- NestedData[, 1 : (ThresModelAUC[1] + 1)]
# Additional version not presented in the manuscript				
		AUCPRs[, 26] <- apply(AUCPRs[, 1 : 25], 1, mean)							
		AUCPRs[, 27] <- apply(AUCPRs[, 1 : 25], 1, sd)							 
		AUCPRMaxModel <- which.max(AUCPRs[, 26])							
		ThresModelPR <- which(AUCPRs[, 26] >= AUCPRs[AUCPRMaxModel, 26] - AUCPRs[AUCPRMaxModel, 27])  
		SelDataAUCPRs <- NestedData[, 1 : (ThresModelPR[1] + 1)]
	}
	Out <- list(SelDataOOBs = SelDataOOBs, SelDataAUCs = SelDataAUCs, SelDataAUCPRs = SelDataAUCPRs)
	return(Out)					
}

# Version RF/OOB

VarImpRF <- function(InData){
	InData$Y <- factor(InData$Y)																
	Imps <- matrix(nrow = ncol(InData)-1, ncol = 55, dimnames = list(names(InData)[-1], c(1 : 50, "mean", "sd", "VarNo", "CartPred", "ThreshPos")))		
	for (i in 1 : 50){																
		rFi <- randomForest(Y ~ ., data = InData, importance = T, ntree = 2000, mtry = round((ncol(InData)-1)/2))	
		Imps[,i] <- randomForest::importance(rFi, scale = F, type = 1)
	}
	Imps[, 51] <- apply(Imps[, 1 : 50], 1, mean)						
	Imps[, 52] <- apply(Imps[, 1 : 50],1, sd)							
	sort <- order(Imps[, 51], decreasing = TRUE)									
	Imps <- Imps[sort,]																			
	Imps[, 53] <- 1 : (length(InData)-1)									
	Imps[, 54] <- predict(rpart(Imps[, 52] ~ Imps[, 53]))					
	Imps[, 55] <- Imps[,51] - min(Imps[, 54])							
	Thresh <- ifelse(all(Imps[, 55] > 0), nrow(Imps), (which(Imps[, 55] < 0)-1)[1]) 												
	NestedData <- cbind(InData[, 1, drop = F], InData[, -1][, sort])																										
	NestedData <- NestedData[, 1 : (Thresh + 1), drop = F]
	if(ncol(NestedData) <= 2){
		SelData <- NestedData
	}else{									
		OOBs <- matrix(ncol = 27, nrow = length(NestedData)-1, dimnames = list(paste0("NoOfVar", 1:(length(NestedData)-1)), c(paste0("Rep", 1:25), "mean", "sd")))					
		for (j in 1 : (length(NestedData)-1)){								
			subs <- NestedData[, 1 : (j+1)]						
			for (k in 1 : 25){				
				rFk <- randomForest(Y ~ ., data = subs)	
				OOBs[j, k] <- rFk$err.rate[500, 1]						
			}								
		}
		OOBs[, 26] <- apply(OOBs[, 1 : 25], 1, mean)	
		OOBs[, 27] <- apply(OOBs[, 1 : 25], 1, sd)									 
		OOBMinModel <- which.min(OOBs[, 26])							
		ThresModel <- which(OOBs[, 26] <= OOBs[OOBMinModel, 26] + OOBs[OOBMinModel, 27]) 
		SelData <- NestedData[, 1 : (ThresModel[1]+1)]	
	}					
	return(SelData)
}


############################## Extract ranks

FunRevRank <- function(x){rank(-x, ties.method = "average")}

FunExtractRanks <- function(MyPath, MyPattern, NoTrue){	
	Files <- list.files(MyPath, pattern = MyPattern)
	Files <- Files[order(Files)]
	RanksX1 <- as.data.frame(matrix(nrow = (length(Files) * 11), ncol = 5))
	colnames(RanksX1) <- c("RFVersion", "Rank", "X", "Beta", "Set")
	if(NoTrue > 1){RanksX2 <- RanksX1}
	if(NoTrue > 2){
		RanksX3 <- RanksX1
		RanksX4 <- RanksX1
	}
	for(i in 1 : length(Files)){
		CurFile <- read.csv(paste0(MyPath, Files[i]))
		CurFile <- CurFile[, c(1, grep("mean", colnames(CurFile)))]
		CurRanks <- apply(CurFile[, -1], 2, FunRevRank)
		CurSettings <- unlist(strsplit(Files[i], split = "_"))[3 : 4]
		IndX1 <- which(CurFile[, 1] == "X1")
		RanksX1$Rank[(i*11-10) : (i*11)] <- CurRanks[IndX1, ]
		RanksX1$RFVersion[(i*11-10) : (i*11)] <- colnames(CurRanks)
		RanksX1$Beta[(i*11-10) : (i*11)] <- CurSettings[1]
		RanksX1$Set[(i*11-10) : (i*11)] <- CurSettings[2]
		RanksX1$X[(i*11-10) : (i*11)] <- "X1"
		if(NoTrue > 1){
			IndX2 <- which(CurFile[, 1] == "X2")
			RanksX2$Rank[(i*11-10) : (i*11)] <- CurRanks[IndX2, ]
			RanksX2$RFVersion[(i*11-10) : (i*11)] <- colnames(CurRanks)
			RanksX2$Beta[(i*11-10) : (i*11)] <- CurSettings[1]
			RanksX2$Set[(i*11-10) : (i*11)] <- CurSettings[2]
			RanksX2$X[(i*11-10) : (i*11)] <- "X2"
		}
		if(NoTrue > 2){
			IndX3 <- which(CurFile[, 1] == "X3")
			RanksX3$Rank[(i*11-10) : (i*11)] <- CurRanks[IndX3, ]
			RanksX3$RFVersion[(i*11-10) : (i*11)] <- colnames(CurRanks)
			RanksX3$Beta[(i*11-10) : (i*11)] <- CurSettings[1]
			RanksX3$Set[(i*11-10) : (i*11)] <- CurSettings[2]
			RanksX3$X[(i*11-10) : (i*11)] <- "X3"
				
			IndX4 <- which(CurFile[, 1] == "X4")
			RanksX4$Rank[(i*11-10) : (i*11)] <- CurRanks[IndX4, ]
			RanksX4$RFVersion[(i*11-10) : (i*11)] <- colnames(CurRanks)
			RanksX4$Beta[(i*11-10) : (i*11)] <- CurSettings[1]
			RanksX4$Set[(i*11-10) : (i*11)] <- CurSettings[2]
			RanksX4$X[(i*11-10) : (i*11)] <- "X4"
		}
	}
	Ranks <- RanksX1
	if(NoTrue > 1){Ranks <- rbind(Ranks, RanksX2)}
	if(NoTrue > 2){Ranks <- rbind(Ranks, RanksX3, RanksX4)}
	PInd <- regexpr("P", Ranks$Set)
	Ranks$Size <- substr(Ranks$Set, 1, PInd-1)
	Ranks$Imb <- substr(Ranks$Set, PInd, nchar(Ranks$Set))
	Ranks$RFVersion <- gsub("_mean", "", Ranks$RFVersion)
	Ranks$ImbPerc <- 0
	Ranks[Ranks$Set == "A300P60", ]$ImbPerc <- 20
	Ranks[Ranks$Set == "A200P40", ]$ImbPerc <- 20
	Ranks[Ranks$Set == "A100P20", ]$ImbPerc <- 20
	Ranks[Ranks$Set == "A50P10", ]$ImbPerc <- 20
	Ranks[Ranks$Set == "A300P45", ]$ImbPerc <- 15
	Ranks[Ranks$Set == "A200P30", ]$ImbPerc <- 15
	Ranks[Ranks$Set == "A100P15", ]$ImbPerc <- 15
	Ranks[Ranks$Set == "A50P8", ]$ImbPerc <- 15
	Ranks[Ranks$Set == "A300P30", ]$ImbPerc <- 10
	Ranks[Ranks$Set == "A200P20", ]$ImbPerc <- 10
	Ranks[Ranks$Set == "A100P10", ]$ImbPerc <- 10
	Ranks[Ranks$Set == "A50P5", ]$ImbPerc <- 10
	Ranks[Ranks$Set == "A300P15", ]$ImbPerc <- 5
	Ranks[Ranks$Set == "A200P10", ]$ImbPerc <- 5
	Ranks[Ranks$Set == "A100P5", ]$ImbPerc <- 5
	Ranks[Ranks$Set == "A50P2", ]$ImbPerc <- 5
	return(Ranks)
}


#################### extract number of true and noise variables after variable selection

ExtractVarSel <- function(SetName, Ys){
	VersionNames <- c("Original", "AUCOOB", "AUCAUC", "AUCAUCPR")
	NAbs <- rep(c(300, 200, 100, 50), each = 4)
	NPres <- c(seq(60, 15, -15), seq(40, 10, -10), seq(20, 5, -5), round(seq(10, 2.5, -2.5)))
	OriginalOut <- as.data.frame(matrix(ncol = 4, nrow = length(NAbs) * 5))
	colnames(OriginalOut) <- c("NAbs", "NPres", "No_X",  "No_Noise")
	OriginalOut$NAbs <- rep(NAbs, each = 5)
	OriginalOut$NPres <- rep(NPres, each = 5)
	AUCOOBOut <- OriginalOut
	AUCAUCOut <- OriginalOut
	AUCAUCPROut <- OriginalOut
	for(k in 1 : length(NAbs)){
		for(m in 1 : 5){
			for(n in 1 : length(VersionNames)){
				CurFileName <- paste("VS", VersionNames[n], SetName, Ys, paste0("A", NAbs[k], "P", NPres[k]), paste0("Sel", m), sep = "_")
				CurFile <- read.csv(paste0(MyPath, CurFileName, ".csv"))
# Version RF/OOB
				if(n == 1){
					OriginalOut$No_X[k*5-5+m] <- length(grep("X", colnames(CurFile)))
					OriginalOut$No_Noise[k*5-5+m] <- length(grep("Noise", colnames(CurFile)))
				}
# Version CFor-AUC/OOB
				if(n == 2){
					AUCOOBOut$No_X[k*5-5+m] <- length(grep("X", colnames(CurFile)))
					AUCOOBOut$No_Noise[k*5-5+m] <- length(grep("Noise", colnames(CurFile)))
				}
# Version CFor-AUC/AUC
				if(n == 3){
					AUCAUCOut$No_X[k*5-5+m] <- length(grep("X", colnames(CurFile)))
					AUCAUCOut$No_Noise[k*5-5+m] <- length(grep("Noise", colnames(CurFile)))
				}
# Additional version not presented in the manuscript
				if(n == 4){
					AUCAUCPROut$No_X[k*5-5+m] <- length(grep("X", colnames(CurFile)))
					AUCAUCPROut$No_Noise[k*5-5+m] <- length(grep("Noise", colnames(CurFile)))
				}
			}
		}
	}
	OriginalOut$Version <- "Original"
	AUCOOBOut$Version <- "AUCOOB"
	AUCAUCOut$Version <- "AUCAUC"
	AUCAUCPROut$Version <- "AUCAUCPR"
	Out <- rbind(OriginalOut, AUCOOBOut, AUCAUCOut, AUCAUCPROut)
	Out$Beta <- Ys
	Out$Set <- paste0("A", Out$NAbs, "P", Out$NPres)
	Out$ImbPerc <- Out$NPres * 100 / Out$NAbs
# only if very small dataset with 50 absences considered:
	Out[Out$ImbPerc == 4, ]$ImbPerc <- 5
	Out[Out$ImbPerc == 16, ]$ImbPerc <- 15
# calculate sensitivity
	NoTrue <- as.numeric(gsub("[A-Za-z]", "", SetName))
	Out$Sens <- Out$No_X / NoTrue
# calculate specificity
	NoNoise <- 121 - NoTrue
	Out$Spec <- (NoNoise - Out$No_Noise) / NoNoise
	return(Out)
}


	