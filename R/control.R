
boost_control <- function(mstop = 100, nu = 0.1, constraint = FALSE, 
                          risk = c("inbag", "oobag", "none"), savedata = TRUE) {

   risk <- match.arg(risk)
   RET <- list(mstop = mstop, nu = nu, constraint = constraint,
               risk = risk, savedata = savedata)
   class(RET) <- c("boost_control")
   RET
}
