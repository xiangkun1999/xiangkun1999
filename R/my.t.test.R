#' My Student's t-Test
#'
#' Performs one and two sample t-tests on vectors of data. Then draws corresponding image.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param conf.level confidence level of the interval.
#' @param normal.te a logical indicating whether you want a normality test.
#'
#' @details The formula interface is only applicable for the 2-sample tests.
#'
#'          alternative = "greater" is the alternative that x has a larger mean than y.
#'
#           If paired is TRUE then both x and y must be specified and they must be the same length. Missing values are silently removed (in pairs if paired is TRUE). If var.equal is TRUE then the pooled estimate of the variance is used. By default, if var.equal is FALSE then the variance is estimated separately for both groups and the Welch modification to the degrees of freedom is used.
#'
#'          If the input data are effectively constant (compared to the larger of the two means) an error is generated.
#' @return  A list with class "htest" containing the following components:
#'
#'          <statistic>	     the value of the t-statistic.
#'
#'          <parameter>	     the degrees of freedom for the t-statistic.
#'
#'          <p.value>	       the p-value for the test.
#'
#'          <conf.int>	     a confidence interval for the mean appropriate to the specified alternative hypothesis.
#'
#'          <estimate>	     the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.
#'
#'          <null.value>	   the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample test or a two-sample test.
#'
#'          <stderr>	       the standard error of the mean (difference), used as denominator in the t-statistic formula.
#'
#'          <alternative>	   a character string describing the alternative hypothesis.
#'
#'          <method>	       a character string indicating what type of t-test was performed.
#'
#'          <data.name>	     a character string giving the name(s) of the data.
#'
#' @export
#' @examples my.t.test(x, y, alternative = "le", paired = TRUE, var.equal = TRUE, conf.level = 0.99)
#' my.t.test(x, mu = 0.5, normal.te = TRUE)
my.t.test<-function (x, y = NULL, alternative = c("two.sided", "less",
                                                  "greater"), mu = 0, paired = FALSE, var.equal = FALSE,
                     conf.level = 0.95,normal.te = FALSE, ...)
{   #Test for variable input
  alternative <- match.arg(alternative) #Allow default matching
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  #If it is a two-parameter test
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    # dname is x and y
    if (paired)
      xok <- yok <- complete.cases(x, y) #If it matches, remove the NA value in x & y and assign it to xok, yok
    #in x & y and assign it to xok, yok
    else {
      yok <- !is.na(y) #yok is the NA value in y
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else { #If it is a single parameter test
    dname <- deparse(substitute(x)) #dname is x
    if (paired)
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if(paired){
    par(mfrow=c(1,2))
    boxplot(x,y,names=c('x','y'),col=c('maroon1','mediumpurple1'),
            main="the plot of two sample t-test",varwidth = T,notch = T)
    if(normal.te) {qqplot(x,y,main = "the plot of normality tests")
      qqline(c(x,y),col="red")
      legend("bottomright",legend=c("theory","sample"),
             pch=list("",1),col=c("red","black"),lwd=c(1,-1))
    }
  } #Double-sample box plot
  if (paired) {
    x <- x - y #Pairwise test if the means are equal
    y <- NULL
  }
  # At this point, only the x variable is left
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  if (is.null(y)) {
    if (nx < 2)
      stop("not enough 'x' observations")
    df <- nx - 1 #Degrees of freedom n-1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx))
      stop("data are essentially constant") #Data is essentially constant
    tstat <- (mx - mu)/stderr #t value
    method <- if (paired) #method is Single sample or paired test
      "Paired t-test"
    else "One Sample t-test"
    estimate <- setNames(mx, if (paired)  #estimate is mean difference or mean
      "mean of the differences"
      else "mean of x")
  }
  else { #y Not empty and unpaired
    ny <- length(y)
    if (nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if (ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if (var.equal && nx + ny < 3)
      stop("not enough observations")
    my <- mean(y)
    vy <- var(y)
    method <- paste(if (!var.equal)
      "Welch", "Two Sample t-test") #welch or two-sample test (equal variance)
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    if (var.equal) {
      df <- nx + ny - 2 #Degrees of freedom
      v <- 0
      if (nx > 1)
        v <- v + (nx - 1) * vx
      if (ny > 1)
        v <- v + (ny - 1) * vy
      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
    }
    else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny -
                                                        1))
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx),
                                                abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df))
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (paired || !is.null(y))
    "difference in means"
  else "mean"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               stderr = stderr, alternative = alternative, method = method,
               data.name = dname)
  #return t,df,p.value
  class(rval) <- "htest"
  if(!is.null(y)&&(!paired)){
    par(mfrow=c(1,2))
    boxplot(x,y,names=c('x','y'),col=c('maroon1','mediumpurple1'),
            main="the plot of two sample t-test",varwidth = T,notch = T)
    if(normal.te){
      qqplot(x,y,main = "the plot of normality tests")
      qqline(c(x,y),col="red")
      legend("bottomright",legend=c("Normality","Sample"),
             pch=list("",1),col=c("red","black"),lwd=c(1,-1))
    }
  } #Double-sample box plot
  if(!paired&&is.null(y)){
    par(mfrow=c(1,2))
    hist(x,main = "the plot of one sample t-test",breaks = "Scott",freq = F)
    abline(v=mean(x),lwd=2,col="red")
    lines(sort(x),dnorm(sort(x),mean=mean(x),sd=sd(x)),col="purple")
    lines(density(x),col="red",lwd=2)#Density curve
    legend("topright",legend=c("Normality","Sample"),col=c("purple","red"),lwd=c(1,1))
    if(normal.te){
      qqnorm(x,main = "the plot of normality tests")
      qqline(x, col=2, lwd=2)
      legend("bottomright",legend=c("theory","sample"),pch=list("",1),
             col=c("red","black"),lwd=c(1,-1))
    } #qq plot for one-sample normality test
  } #Single sample draw histogram
  rval
}
