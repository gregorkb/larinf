#' Face temperature data
#'
#' The data set contains for each of 933 subjects an oral temperature reading,
#' which we use as the response, as well as several measurements derived from infrared
#' thermographs of the subject's face, the values of some demographic variables,
#' the subject's distance to the thermal camera, the relative humidity, the ambient temperature,
#' and an indicator of whether the subject was wearing cosmetics, for a total of 23 covariates.
#'
#' @docType data
#'
#' @usage data(facetemp)
#'
#' @keywords datasets
#'
#' @references
#'
#' Wang, Q., Zhou, Y., Ghassemi, P., Chenna, D., Chen, M., Casamento, J., Pfefer, J., & Mcbride, D. (2023).
#' Facial and oral temperature data from a large set of human subject volunteers (version 1.0.0). *PhysioNet*.
#' [https://doi.org/10.13026/3bhc-9065](https://doi.org/10.13026/3bhc-9065).
#'
#' Wang, Q., Zhou, Y., Ghassemi, P., McBride, D., Casamento, J. P., & Pfefer, T. J. (2022).
#' Infrared Thermography for Measuring Elevated Body Temperature: Clinical Accuracy, Calibration, and Evaluation.
#' Sensors, 22, 215.
#'
#' Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000).
#' PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals.
#' Circulation (Online). 101 (23), pp. e215–e220.
#'
#' @source PhysioNet: [https://physionet.org/content/face-oral-temp-data/1.0.0/](https://physionet.org/content/face-oral-temp-data/1.0.0/)
#'
#'
#' @examples
#' data(facetemp)
#' X <- facetemp$X
#' y <- facetemp$y
#'
#' lar_out <- lar(X,y)
#' larinf_out <- larinf(X,y)
#'
#' plot(lar_out,which_lab = which(larinf_out$rej))
#' plot(larinf_out,omaadd=c(0,0,1,1))
#'
#' \dontrun{
#' # This code shows how the facetemp data are constructed from the file 'FLIR_groups1and2.csv',
#' # downloadable at https://physionet.org/content/face-oral-temp-data/1.0.0/
#' data0 <- read.csv("facetemp_data/FLIR_groups1and2.csv",skip = 2, header = TRUE)
#'
#' # remove some columns
#' colrm <- c(1,2,30,58,86,114,115,124,125)
#' data1 <- data0[,-colrm]
#'
#' # remove rows with missing data
#' data2 <- data1[apply(!is.na(data1),1,all),]
#'
#' # average the measurements from the four images
#' X1 <- data2[,1:27]
#' X2 <- data2[,28:54]
#' X3 <- data2[,55:81]
#' X4 <- data2[,82:108]
#' X <- as.matrix((X1 + X2 + X3 + X4) / 4)
#'
#' # keep only these measurements
#' Tkeep <- c("T_LC_Dry1","T_LC_Wet1","T_RC_Dry1","T_RC_Wet1","T_FHCC1",
#'            "T_FHRC1","T_FHLC1","T_FHBC1","T_FHTC1","T_OR1")
#'
#' X <- X[,Tkeep]
#'
#' # update names of columns of X
#' colnames(X) <- c("T_LC_Dry","T_LC_Wet","T_RC_Dry","T_RC_Wet","T_FHCC",
#'                  "T_FHRC","T_FHLC","T_FHBC","T_FHTC","T_OR")
#'
#' # collect some more covariates
#' W <- cbind(data2$T_atm,
#'            data2$Humidity,
#'            data2$Distance,
#'            as.numeric(data2$Gender == "Male"),
#'            as.numeric(data2$Ethnicity == "Black or African-American"),
#'            as.numeric(data2$Ethnicity == "White"),
#'            as.numeric(data2$Ethnicity == "Asian"),
#'            as.numeric(data2$Ethnicity == "Hispanic/Latino"),
#'            as.numeric(data2$Cosmetics),
#'            as.numeric(data2$Age == "18-20"),
#'            as.numeric(data2$Age == "21-25"),
#'            as.numeric(data2$Age == "26-30"),
#'            as.numeric(data2$Age == "31-40"))
#'
#' colnames(W) <- c("T_atm","Humidity","Distance","Male","Black","White","Asian",
#'                  "Hisp","Cosmetics","Age18_20","Age21_25","Age26_30","Age31_40")
#'
#' define centered design matrix
#' n <- nrow(X)
#' Xn <- (diag(n) - matrix(1/n,n,n)) %*% cbind(X,W)
#'
#' # define centered response vector
#' yn <- data2$aveOralM - mean(data2$aveOralM)
#'
#' # export the data set as a list
#' facetemp <- list(Xn = Xn, yn = yn)
#' }
"facetemp"
