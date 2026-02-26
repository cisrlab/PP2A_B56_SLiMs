

#' Title
#'
#' @param fs
#' @param fn
#' @param rname
#' @param transform_alg
#' @param feature_select
#' @param regression_algorithm
#' @param cv
#'
#' @return
#' @export
#'
#' @examples
runAnalysis<-function(
        fs,
        fn,
        rname = "ITC",
        transform_alg = "none",
        feature_select="all",
        regression_algorithm = "rforest",
        cv = FALSE
) {
    fs2 = fs;
    fs2$Transform = NA;
    #fs2$Transform_SD = fs2$ITC_SD;

    fs2$Transform = doTransform(fs[,rname], transform_alg)
    #print(fs2);
    #print(summary(fs2))
    tt.res = runFSTraining(fs2, fn, "Transform", feature_select, regression_algorithm, cv=cv, transform_alg = transform_alg);

    p = NA;
    if (cv) {
        p = tt.res$predicts
    } else {
        p = tt.res$test.predict;
    }

    p = undoTransform(p, transform_alg);

    #Enforce a cap;
    top_cap = max(fs$ITC, na.rm=TRUE) * 2;
    bottom_cap = min(fs$ITC, na.rm=TRUE) / 2;

    p[p > top_cap] = top_cap;
    p[p < bottom_cap] = bottom_cap;

    tt_mse = mse(fs$ITC, p)
    tt_mae = mae(fs$ITC, p)
    tt_rmse = rmse(fs$ITC, p)
    tt_nmae = nmae(fs$ITC, p)
    tt_mperr = mperr(fs$ITC, p)
    tt_med_mperr = med.perr(fs$ITC,p)

    tt_mse_hv = mse(fs$ITC[fs$HasValue],p[fs$HasValue])
    tt_mae_hv = mae(fs$ITC[fs$HasValue],p[fs$HasValue])
    tt_rmse_hv = rmse(fs$ITC[fs$HasValue], p[fs$HasValue])
    tt_nmae_hv = nmae(fs$ITC[fs$HasValue], p[fs$HasValue])
    tt_mperr_hv = mperr(fs$ITC[fs$HasValue], p[fs$HasValue])
    tt_med_mperr_hv = med.perr(fs$ITC[fs$HasValue],p[fs$HasValue])


    ans = data.frame(
        Transform = transform_alg,
        FeatureSelect = feature_select,
        Regressor = regression_algorithm,
        MSE = tt_mse,
        MAE = tt_mae,
        RMSE = tt_rmse,
        NMAE = tt_nmae,
        MPERR = tt_mperr,
        MED_PERR = tt_med_mperr,
        MSE_HV = tt_mse_hv,
        MAE_HV = tt_mae_hv,
        RMSE_HV = tt_rmse_hv,
        NMAE_HV = tt_nmae_hv,
        MPERR_HV = tt_mperr_hv,
        MED_PERR_HV = tt_med_mperr_hv,
        stringsAsFactors = FALSE
    );
    attr(ans, "rname") = rname;
    attr(ans, "transform_alg") = transform_alg;
    attr(ans, "feature_select") = feature_select;
    attr(ans, "regression_algorithm") = regression_algorithm;
    attr(ans, "fs") = fs;
    attr(ans, "fs2") = fs2;
    attr(ans, "tt.res") = tt.res;
    attr(ans, "actual") = fs$ITC;
    attr(ans, "predict") = p;
    attr(ans, "top_cap") = top_cap;
    attr(ans, "bottom_cap") = bottom_cap;


    return(ans);
}


getTrainTestFxn<-function(regression_algorithm, transform_alg="none") {
    if (regression_algorithm == "rforest") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(rForestTrainTest(train_data,test_data,class_name, feature_names, na.action = na.omit, importance=TRUE))
            }
        )
    }
    if (regression_algorithm == "rfsrc")
    {
        return (rfsrcTrainTestReg);
    }
    if (regression_algorithm == "rfsrc_w1")
    {
        return(rfsrcTrainTestW1Reg);
    }
    if (regression_algorithm == "rfsrc_w2")
    {
        return(rfsrcTrainTestW2Reg);
    }

    if (regression_algorithm == "lm") {
        return(lmTrainTestReg);
    }

    if (regression_algorithm == "lmi") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(lmTrainTestReg(train_data, test_data, class_name, feature_names, transform_fxn = transform_alg, var.interact = TRUE, ...))
            })

    }

    if (regression_algorithm == "lm_wp") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(lmTrainTestWPReg(train_data, test_data, class_name, feature_names, transform_fxn = transform_alg, var.interact = FALSE, ...))
            })
    }

    if (regression_algorithm == "lm_wfc") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(lmTrainTestWEReg(train_data, test_data, class_name, feature_names, transform_fxn = transform_alg, var.interact = FALSE, err_fxn = fold.err, ...))
            })
    }


    if (regression_algorithm == "lmi_wp") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(lmTrainTestWPReg(train_data, test_data, class_name, feature_names, transform_fxn = transform_alg, var.interact = TRUE, ...))
            })
    }

    if (regression_algorithm == "lmi_wfc") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(lmTrainTestWEReg(train_data, test_data, class_name, feature_names, transform_fxn = transform_alg, var.interact = TRUE, err_fxn = fold.err, ...))
            })
    }


    if (regression_algorithm == "rfsrc_wz") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(rfsrcTrainTestWZReg(train_data, test_data, class_name, feature_names,
                                           transform_alg = transform_alg, ...))
            })
    }
    if (regression_algorithm == "rfsrc_wz2")
    {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(rfsrcTrainTestWMPReg(train_data, test_data, class_name, feature_names,
                                            transform_alg = transform_alg, square = TRUE, ...))
            }
        )
    }
    if (regression_algorithm == "rfsrc_wp")
    {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(rfsrcTrainTestWMPReg(train_data, test_data, class_name, feature_names,
                                            transform_alg = transform_alg, ...))
            }
        )
    }

    if (regression_algorithm == "rfsrc_wfc")
    {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(
                    rfsrcTrainTestWEReg(
                        train_data, test_data, class_name, feature_names,
                        transform_alg = transform_alg, err_fxn = fold.err, ...
                    )
                )
            }
        )
    }

    if (regression_algorithm == "rfsrc_wp2") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(rfsrcTrainTestWMPReg(train_data, test_data, class_name, feature_names,
                                            transform_alg = transform_alg, square = TRUE, ...))
            }
        )
    }

    if (regression_algorithm == "glmnet_v2") {
        return(PP2A.B56.SLiMs::glmNetTrainTestReg2)
    }
    if (regression_algorithm == "glmnet_wp") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(glmNetTrainTestWPReg(
                    train_data,
                    test_data,
                    class_name,
                    feature_names,
                    transform_fxn = transform_alg,
                    ...))
            }
        )
    }

    if (regression_algorithm == "glmnet_wfc") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                return(glmNetTrainTestWEReg(
                    train_data,
                    test_data,
                    class_name,
                    feature_names,
                    transform_fxn = transform_alg,
                    err_fxn = fold.err,
                    ...))
            }
        )
    }



    if (regression_algorithm == "glmnet_wz") {
        if (transform_alg == "log") {
            return(glmNetTrainTestWZRegLog)
        }
        if (transform_alg == "sqrt")
        {
            return(glmNetTrainTestWZRegSqrt)
        }
        return(glmNetTrainTestWZReg)
    }

    if (regression_algorithm == "glmnet_wz2") {
        if (transform_alg == "log") {
            return(glmNetTrainTestWZ2RegLog)
        }
        if (transform_alg == "sqrt")
        {
            return(glmNetTrainTestWZ2RegSqrt)
        }
        return(glmNetTrainTestWZ2Reg)
    }

    if (regression_algorithm == "M5P") {
        return(M5PTrainTest);
    }
    if (regression_algorithm == "cforest") {
        return(function(train_data, test_data, class_name, feature_names, ...) {
            return(cForestTrainTest(train_data, test_data, class_name, feature_names, ntree=500, regression=TRUE,...))
        });
    }
    if (regression_algorithm == "glmnet") {
        return(glmNetTrainTestReg);
    }
    if (regression_algorithm == "ridge") {
        return(glmNetTrainTestRegRidge)
    }
    if (regression_algorithm == "glmnet_w1") {
        return(glmNetTrainTestW1Reg);
    }

    if (regression_algorithm == "glmnet_w2") {
        return(glmNetTrainTestW2Reg);
    }

    if (regression_algorithm == "glm") {
        return(glmTrainTestReg);
    }

    stop("Unknown regression:", regression_algorithm);
}

getTrainTestFxn2<-function(train_test_fxn, feature_select = "all") {
    if (feature_select == "all") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                fn <- feature_names[feature_names != "Substrate"]
                #print(rule_fn)
                res = train_test_fxn(train_data, test_data, class_name, fn, ...)
                return(res);
            }
        )
    }
    if (feature_select == "all_no_meme") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                fn <- feature_names[feature_names != "Substrate"]
                fn <- fn[fn != "Clip_meme_Ba"]
                #print(rule_fn)
                res = train_test_fxn(train_data, test_data, class_name, fn, ...)
                return(res);
            }
        )
    }
    if (feature_select == "all_s") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                fn <- unique(c(feature_names, "Substrate"))
                #print(rule_fn)
                res = train_test_fxn(train_data, test_data, class_name, fn, ...)
                return(res);
            }
        )
    }
    if (feature_select == "rules") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                rule_fn = c(feature_names[grep("Rule", feature_names)])
                #print(rule_fn)
                res = train_test_fxn(train_data, test_data, class_name, rule_fn, ...)
                return(res);
            }
        )
    }

    if (feature_select == "rules_s") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                rule_fn = unique(c(feature_names[grep("Rule", feature_names)], "Substrate"))
                #print(rule_fn)
                res = train_test_fxn(train_data, test_data, class_name, rule_fn, ...)
                return(res);
            }
        )
    }


    if (feature_select == "rules_cc") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                rule_fn = c(feature_names[grep("Rule", feature_names)], "Clip_charge");
                res = train_test_fxn(train_data, test_data, class_name, rule_fn, ...)
                return(res);
            }
        )
    }
    if (feature_select == "rules_cc_s") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                rule_fn = unique(c(feature_names[grep("Rule", feature_names)], "Clip_charge", "Substrate"))
                res = train_test_fxn(train_data, test_data, class_name, rule_fn, ...)
                return(res);
            }
        )
    }

    if (feature_select == "norules") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                rule_fn = c(feature_names[grep("Rule", feature_names, invert=TRUE)]);
                rule_fn <- rule_fn[rule_fn != "Substrate"]
                res = train_test_fxn(train_data, test_data, class_name, rule_fn, ...)
                return(res);
            }
        )
    }
    if (feature_select == "norules_s") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                rule_fn = unique(c(feature_names[grep("Rule", feature_names, invert=TRUE)], "Substrate"))
                res = train_test_fxn(train_data, test_data, class_name, rule_fn, ...)
                return(res);
            }
        )
    }
    if (feature_select == "esm2") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                fn = unique(c(feature_names[grep("esm2", feature_names, invert=FALSE)]))
                res = train_test_fxn(train_data, test_data, class_name, fn, ...)
                return(res)
            }
        )
    }
    if (feature_select == "all_no_esm2") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                fn = unique(c(feature_names[grep("esm2", feature_names, invert=TRUE)]))
                fn <- fn[fn != "Substrate"]
                res = train_test_fxn(train_data, test_data, class_name, fn, ...)
                return(res)
            }
        )

    }
    if (feature_select == "all_no_meme_no_esm2") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
              fn <- feature_names[grep("esm2", feature_names, invert=TRUE)]
              fn <- fn[fn != "Substrate"]
              fn <- fn[fn != "Clip_meme_Ba"]
              res <- train_test_fxn(train_data, test_data, class_name, unique(fn), ...)
              return(res)
            }
        )
    }
    if (feature_select == "zpep") {
        return(
            function(train_data, test_data, class_name, feature_names, ...) {
                zpep_fn = c(feature_names[grep("zpep", feature_names, invert=FALSE)]);
                res = train_test_fxn(train_data, test_data, class_name, zpep_fn, ...)
                return(res);
            }
        )
    }


    if (feature_select == "forward_mse") {
        return (
            function(train_data, test_data, class_name, feature_names, ...) {
                res = forwardTrainTest(
                    train_data,
                    test_data,
                    class_name,
                    feature_names,
                    train.nfold = -1,
                    metric = "MSE",
                    smaller_is_better = TRUE,
                    forward_train_test_fxn = train_test_fxn,
                    debug=TRUE,
                    ...
                )
                return(res);
            }
        );
    }
    if (feature_select == "forward_mperr") {
        return (
            function(train_data, test_data, class_name, feature_names, ...) {
                res = forwardTrainTest(
                    train_data,
                    test_data,
                    class_name,
                    feature_names,
                    train.nfold = -1,
                    metric = "MPERR",
                    smaller_is_better = TRUE,
                    forward_train_test_fxn = train_test_fxn,
                    debug=TRUE,
                    ...
                )
                return(res);
            }
        );
    }

    if (feature_select == "rfs_mperr") {
        return (
            function(train_data, test_data, class_name, feature_names, ...) {
                res = rForestScoredTrainTest(
                    train_data,
                    test_data,
                    class_name,
                    feature_names,
                    train.nfold = -1,
                    metric = "MPERR",
                    rf.metric = "%IncMSE",
                    smaller_is_better = FALSE,
                    metric_smaller_is_better = TRUE,
                    scored_train_test_fxn = train_test_fxn,
                    debug=FALSE,
                    filestem = "rfs_mperr",
                    ...
                )
                return(res);
            }
        );
    }

    stop("Unknown:", feature_select);
}


runFSTraining<-function(
        fs,
        fn,
        rname,
        feature_select= "all",
        regression_algorithm = "rforest",
        transform_alg = "none",
        cv = FALSE
) {

    tt_fxn = getTrainTestFxn(regression_algorithm, transform_alg = transform_alg);
    tt_fs_fxn = getTrainTestFxn2(tt_fxn, feature_select);

    if (cv) {
        res = getCVAccuracy(fs, rname, fn, tt_fs_fxn, nfold=-1)
    } else {
        res = tt_fs_fxn(fs, fs, rname, fn)
    }
    res$fn <- fn
    return(res);
}

#' Title
#'
#' @param fs
#' @param remove_uninformative
#'
#' @return
#' @export
#'
#' @examples
extractFN<-function(fs, remove_uninformative=FALSE) {

    fn = colnames(fs)

    fn = fn[fn != "Clip_Sequence"]
    fn = fn[fn != "Full_Sequence"]

    #fn = fn[fn != "sequence_length"]
    fn = fn[fn != "Clip_meme_Ba_idx"]
    fn = fn[fn != "Full_meme_Ba_idx"]
    fn = fn[fn != "ITC"]
    fn = fn[fn != "ITC_SD"]
    fn = fn[fn != "Transform"]
    fn = fn[fn != "Transform_SD"]
    fn = fn[fn != "lITC"]
    fn = fn[fn != "sITC"]
    fn = fn[fn != "zITC"]
    fn = fn[fn != "HasValue"]
    fn = fn[fn != "Orig_Sequence"]
    fn = fn[fn != "SAlpha"]
    fn <- fn[fn != "ITC_Missing_Status"]
    fn <- fn[fn != "Substrate"]
    fn <- fn[fn != "Researcher"]
    #Remove uninformative-features
    to_remove = c();
    to_keep = c();
    for (feature in fn) {
        if (is.logical(fs[,feature])) {
            if (length(table(fs[,feature]))<2) {
                to_remove <- c(to_remove, feature)
            } else {
                to_keep <- c(to_keep, feature)
            }
        }
        if (is.factor(fs[,feature])) {
            if (length(table(fs[,feature]))<nlevels(fs[,feature])) {
                to_remove <- c(to_remove, feature)
            } else {
                to_keep <- c(to_keep, feature)
            }
        }
        if (is.numeric(fs[,feature])) {
            if (length(unique(fs[,feature]))==1) {
                to_remove <- c(to_remove, feature)
            } else {
                to_keep <- c(to_keep, feature)
            }
        }
    }


    #message("Marked ", length(to_remove), " uninformative features");
    if (remove_uninformative) {
        fn = to_keep
    }
    attr(fn, "to_remove") = to_remove
    attr(fn, "to_keep") <- to_keep
    return(fn)
}


getITCFeatures<-function(res, fn) {
    if (res$FeatureSelect[1] == "all") {
        itc_features = fn[fn != "Substrate"]
    } else if (res$FeatureSelect[1] == "all_no_meme") {
        itc_features <- fn[fn != "Substrate"]
        itc_features <- itc_features[itc_features != "Clip_meme_Ba"]
    } else if(res$FeatureSelect[1] == "all_s") {
        itc_features <- unique(c(fn, "Substrate"))
    }
    else if (res$FeatureSelect[1] == "rules") {
        itc_features = fn[grep("Rule", fn)]
    } else if (res$FeatureSelect[1] == "rules_s") {
        itc_features <- unique(c(fn[grep("Rule", fn)], "Substrate"))
    } else if (res$FeatureSelect[1] == "norules") {
        itc_features = fn[grep("Rule", fn, invert=TRUE)]
        itc_features <- fn[fn != "Substrate"]
    } else if (res$FeatureSelect[1] == "norules_s") {
        itc_features = unique(c(fn[grep("Rule", fn, invert=TRUE)], "Substrate"))
    } else if (res$FeatureSelect[1] == "rules_cc") {
        itc_features = c(fn[grep("Rule", fn)], "Clip_charge");
    } else if (res$FeatureSelect[1] == "rules_cc_s") {
        itc_features <- unique(c(fn[grep("Rule", fn)], "Clip_charge", "Substrate"))
    } else if (res$FeatureSelect[1] == "rfs_mperr"){
        tt.res = attr(res, "tt.res")
        itc_features = tt.res$forward.model$selected.features;
    } else if (res$FeatureSelect[1] == "esm2") {
        itc_features = unique(c(fn[grep("esm2", fn, invert=FALSE)]))
    } else if (res$FeatureSelect[1] == "all_no_esm2") {
        itc_features = unique(c(fn[grep("esm2", fn, invert=TRUE)]))
        itc_features <- itc_features[itc_features != "Substrate"]
    } else if (res$FeatureSelect[1] == "all_no_meme_no_esm2") {
        itc_features <- fn[fn != "Substrate"]
        itc_features <- itc_features[grep("esm2", itc_features, invert=TRUE)]
        itc_features <- itc_features[itc_features != "Clip_meme_Ba"]
        itc_features <- unique(itc_features)
    } else {
        stop("Unknown feature set:", res$FeatureSelect[1]);
    }
    return(itc_features);
}


