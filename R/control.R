
boost_control <- function(mstop = 100, nu = 0.1, constraint = FALSE, 
                          risk = TRUE, savedata = TRUE) {

   RET <- list(mstop = mstop, nu = nu, constraint = constraint,
               risk = risk, savedata = savedata)
   class(RET) <- c("boost_control")
   RET
}
