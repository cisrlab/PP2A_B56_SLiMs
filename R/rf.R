


rForestVarImpPlots<-function(x, pdf.width=7,pdf.height=7,..., filestem) {

  if (!missing(filestem)) {
    pdf(makeExtFilePath(filestem, "rf.varimp", "pdf"), width=pdf.width, height=pdf.height);
    varImpPlot(x=x,...);
    dev.off();
    png(makeExtFilePath(filestem, "rf.varimp", "png"));
    varImpPlot(x=x,...);
    dev.off();

  }
  ans = varImpPlot(x);
  if (!missing(filestem)) {
    write.csv(ans, makeExtFilePath(filestem, "rf.varimp", "csv"));
  }
  invisible(ans);
}


rForestPlotTuneVars<-function(
  tune.vars,
  metric="AUC.ROC",
  filestem
) {

  if (!missing(filestem)) {
    pdf(makeExtFilePath(filestem, "rforest.tune.vars","pdf"))
  }

  par(mfrow=c(1,2))
  plot(tune.vars$ntree, tune.vars[,metric]);
  plot(tune.vars$mtry, tune.vars[,metric]);


  if (!missing(filestem)) {
    dev.off();
  }
}


rForestTune<-function(
  train_data,
  class_name,
  feature_names,
  nfold=-1,
  nrepeat = 10,
  fold.vec=getFoldVector(train_data, class_name, nfold),
  metric = "AUC.ROC",
  ntrees = c(1:5,10, 25, 50, 100, 200, 500, 1000),
  mtries = 1:length(feature_names),
  sort.it = TRUE,
  debug=FALSE
  ) {

  sort.order = TRUE;
  if (metric == "AUC.ROC.greater.random.pvalue") {
    sort.order=FALSE;
  }


  summary.sum.df = data.frame();
  summary.full.df = data.frame();

  for (idx in 1:nrepeat) {
    fold.vec = getFoldVector(train_data, class_name, nfold);
    summary.df = data.frame();
    for (ntree in ntrees) {
      for (mtry in mtries) {
        if (debug) {
            cat("ntree:",ntree," mtry:",mtry,"\n");
        }
        accuracys = getCVAccuracyFolds(train_data, class_name, feature_names, rForestTrainTest, fold.vec,
          tune=FALSE,
          ntree=ntree,
          mtry=mtry,
          importance=FALSE
        );



        sub.df = data.frame(
            ntree=ntree,
            mtry=mtry,
            stringsAsFactors=FALSE
        );
        sub.df = cbind(sub.df, getCVAccuracySummary(accuracys));
        if (debug) { print(sub.df);}
        #print(sub.df);
        if (nrow(summary.df)) {
          summary.df = rbind(summary.df, sub.df);
        } else {
          summary.df = sub.df;
        }

        if (debug) {
          rForestPlotTuneVars(summary.df, metric=metric);
        }

      }
    }
    if (nrow(summary.sum.df)) {
      summary.sum.df[,3:ncol(summary.sum.df)] = summary.sum.df[,3:ncol(summary.sum.df)] + summary.df[,3:ncol(summary.sum.df)];
    } else {
      summary.sum.df = summary.df;

    }
    summary.df$Repeat = idx;
    if (nrow(summary.full.df) == 0) {
      summary.full.df = summary.df;
    } else {
      summary.full.df = rbind(summary.full.df, summary.df);
    }
    cat("Repeat ",idx, " finished\n");
  }

  summary.sum.df[,3:ncol(summary.sum.df)] = summary.sum.df[,3:ncol(summary.sum.df)] / nrepeat;

  rForestPlotTuneVars(summary.sum.df, metric=metric);

  ans = summary.sum.df;
  if (sort.it) {
    ans = ans[order(ans[,metric], decreasing=sort.order),];
  }

  attr(ans, "summary.full.df") = summary.full.df;

  return(ans);
}



#' Trains and validates a random forest from the randomForest class
#'
#' @param train_data
#' @param test_data
#' @param class_name
#' @param feature_names
#' @param ntree
#' @param mtry
#' @param tune
#' @param importance
#' @param na.action
#' @param debug
#' @param pans
#'
#' @return
#' @export
#'
#' @examples
rForestTrainTest <- function(
  train_data,
  test_data,
  class_name,
  feature_names,
  ntree = 100,
  mtry=sqrt(length(feature_names)),
  tune=FALSE,
  importance=FALSE,
  na.action = randomForest::na.roughfix,
  debug = FALSE,
  pans = list()
  ) {
   #cat("rForestTrainTest\n");
   ans = rForestTrain(train_data, class_name, feature_names, ntree = ntree, mtry=mtry, na.action = na.action, tune=tune, importance=importance, pans = pans);
   #cat("rForestTest\n");
   ans = rForestTest(ans$model, test_data, feature_names, pans = ans);
   return(ans);

}

rForestAssignPValues<-function(imp, debug=FALSE) {
  if (debug) {print(colnames(imp));}
  if ("MeanDecreaseAccuracy" %in% colnames(imp)) {
    imp$MDA.pvalue = pnorm(-imp$MeanDecreaseAccuracy);
    imp$MDA.padj = p.adjust(imp$MDA.pvalue, method="BH");
    imp = imp[order(imp$MDA.pvalue, decreasing=FALSE),];
  }
  if ("%IncMSE" %in% colnames(imp)) {
    imp$PIMSE.pvalue = pnorm(-imp[,"%IncMSE"]);
    imp$PIMSE.padj = p.adjust(imp$PIMSE.pvalue, method="BH");
    imp = imp[order(imp$PIMSE.pvalue, decreasing=FALSE),];
  }
  return(imp);
}

################
# Trains a random forest model
# Returns in a list
#   model - learned model
#   importance - importance dataframe
#   train.prob - train probabilities
#   train.predict - training predictions
################
rForestTrain<-function(train_data, class_name, feature_names,
  ntree=500,
  mtry=max(1,as.integer(sqrt(length(feature_names)))),
  tune=FALSE,
  tune.nfold=10,
  importance=FALSE,
  na.action = randomForest::na.roughfix,
  debug=FALSE,
  seed = 1003,
  pans = list()
  ) {
   ans = pans;
   regression = is.numeric(train_data[,class_name])

   #cat("class_name:", class_name,"\n");
   #cat("feature_names:\n");
   #print(feature_names);

   #cat("rforestTrain\n");
   library(randomForest);
   #cat("Training model\n");

   if (length(feature_names) == 1) {
     train_data$temp_ = rnorm(nrow(train_data));
     feature_names = c(feature_names, "temp_");
   }

   class_values = train_data[,c(class_name)];
   if (!regression) {
     class_values = factor(class_values)
   }
  if (tune) {
    ans$tune.vars = rForestTune(train_data, class_name, feature_names, tune.nfold);
    print(head(ans$tune.vars))
    set.seed(seed)
    ans$model = randomForest(class_values ~.,
      data=train_data[,feature_names],
      ntree=ans$tune.vars$ntree[1],
      mtry=ans$tune.vars$mtry[1],
      na.action = na.action,
      importance=importance #This will calculate the importance dataframe

      );
  } else {
    set.seed(seed)
    ans$model = randomForest(class_values ~.,
      data=train_data[,feature_names],
      ntree=ntree,
      mtry=mtry,
      na.action = na.action,
      importance=importance #This will calculate the importance dataframe
      );
   }
   #cat("Getting predictions\n");
   model = ans$model;

   if (importance) {
     ans$importance = as.data.frame(importance(model,scaled=TRUE));
     ans$importance = rForestAssignPValues(ans$importance);  #TODO fix assignPValues to handle regression
     ans$importanceSD = model$importanceSD;
     ans$importance.unscaled = as.data.frame(model$importance);
     ans$importance.unscaled = rForestAssignPValues(ans$importance.unscaled);

   }
   ans$train.actual = class_values;
   ans$train.predict = predict(model, train_data[,feature_names]);

   if (!regression) {
     ans$train.prob = predict(model, train_data[,feature_names], type="prob");
   } else {
     temp = ans$train.actual - ans$train.predict;
     ans$train.mse = sum(temp * temp) / length(temp);
     ans$train.mae = sum(abs(temp)) / length(temp);
     ans$train.sstot = ans$train.actual - mean(ans$train.actual);
     ans$train.sstot = sum(ans$train.sstot * ans$train.sstot);
     ans$train.ssreg = sum(temp * temp);
     ans$train.rsq = 1 - (ans$train.ssreg/ans$train.sstot);


   }
   return(ans);
}

rForestTest<-function(model, test_data, feature_names, pans = list()) {
  library(randomForest)
  #cat("rForestTest: start\n");
  if (length(feature_names) == 1) {
    test_data$temp_ = rnorm(nrow(test_data));
    feature_names = c(feature_names, "temp_");
  }
  ans = pans;

   regression = model$type == "regression";
   ans$test.predict = predict(model, test_data[,feature_names]);
   if (!regression) {
     ans$test.prob = try(predict(model, test_data[,feature_names], type="prob"));
     ans$test.prob <- as.data.frame(ans$test.prob)
     rownames(ans$test.prob) = rownames(test_data);
   } else {
     ans$test.ci = predict(model, test_data[,feature_names], interval="confidence");

   }
   return(ans);

}



