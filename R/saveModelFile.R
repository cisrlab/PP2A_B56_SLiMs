#' Title
#'
#' @param res
#' @param model_file_path
#'
#' @return
#' @export
#'
#' @examples
saveModelFile<-function(res, model_file_path) {

    tt.res = attr(res, "tt.res")
    itc_model = tt.res$model;
    itc_transform = res$Transform[1];
    itc_regressor = res$Regressor[1];
    itc_label = paste0(res$Label[1],"_",date_tag);
    itc_get_fs_fxn = getFS2b

    itc_top_cap = attr(res,"top_cap")
    itc_bottom_cap = attr(res, "bottom_cap")
    itc_features <- res$FeatureSelect
    itc_tt.res <- tt.res
    fn <- attr(tt.res, "fn")

    itc_features <- getITCFeatures(res,fn);

    itc_train_predict <- attr(res, "predict")

    #print(itc_label);
    save(list = c(
        "itc_model",
        "itc_transform",
        "itc_regressor",
        "itc_features",
        "itc_label",
        "itc_get_fs_fxn",
        "itc_top_cap",
        "itc_bottom_cap",
        "itc_tt.res",
        "itc_train_predict"
        ), file = model_file_path)


}
