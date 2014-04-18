\docType{package}
\name{quantspec-package}
\alias{quantspec}
\alias{quantspec-package}
\title{Quantile-based Spectral Analysis of Time Series}
\description{
Methods to determine, smooth and plot quantile (i. e.,
Laplace or copula) periodograms for univariate time series.
}
\details{
\tabular{ll}{ \cr Package: \tab quantspec \cr Type: \tab
Package \cr Version: \tab 1.0-0 \cr Date: \tab 2014-04-18
\cr License: \tab GPL (>= 2) }
}
\section{Contents}{
  The \pkg{quantspec} package contains an hierachy of S4
  classes with corresponding methods and functions serving
  as constructors. The following class diagrams provide an
  overview on the structure of the package. In the first
  class diagram the classes implementing the estimators are
  implemented.

  \if{html}{\figure{main.png}{options: width=960}}
  \if{latex}{\figure{main.pdf}{options: width=15cm}}

  In the second class diagram the classes implementing
  model quantities are displayed. A relation to the
  ``empirical classes'' is given via the fact that the
  quantile spectral densities are computed by simulation of
  quantile periodograms and a common abstract superclass
  \code{QSpecQuantity} which is used to provide a common
  interface to quantile spectral quantities.

  \if{html}{\figure{csd.png}{options: width=768}}
  \if{latex}{\figure{csd.pdf}{options: width=12cm}}

  Besides the object-oriented design a few auxiliary
  functions exists. They serve as parameters or are mostly
  for internal use. A more detailed description of the
  framework can be found in the paper on the package (Kley,
  2014b).
}

\section{Organization of the source code / files in the \code{/R} folder}{
  All of the source code related to the specification of a
  certain class is contained in a file named
  \code{Class-[Name_of_the_class].R}. This includes, in the
  following order, \enumerate{ \item all roxygen
  \code{@include} to insure the correctly generated collate
  for the DESCRIPTION file.  \item \code{\\setClass}
  preceded by a meaningful roxygen documentation.  \item
  specification of an \code{initialize} method, where
  appropriate.  \item all accessor and mutator method (i.
  e., getter and setter); first the ones returning
  attributes of the object, then the ones returning
  associated objects.  \item constructors; use generics if
  there is more than one of them.  \item \code{show} and
  \code{plot} methods. }
}

\section{Coding Conventions}{
  To improve readability of the software and documentation
  this package was written in the spirit of the ``Coding
  conventions of the Java Programming Language'' (Oracle,
  2014). In particular, the naming conventions for classes
  and methods have been adopted, where ``Class names should
  be nouns, in mixed case with the first letter of each
  internal word capitalized.'' and ``Methods should be
  verbs, in mixed case with the first letter lowercase,
  with the first letter of each internal word
  capitalized.''
}

\section{Naming Conventions for the Documentation}{
  To reflect the structure of the contents of the package
  in the documentation file, the following system for
  naming of the sections is adopted: \itemize{ \item
  Documentation of an S4 class is named as the name of the
  class followed by ``-class''. [cf.
  \code{\link{QuantilePG-class}}] \item Documentation of a
  constructor for an S4-class is named as the name of the
  class followed by ``-constructor''.  [cf.
  \code{\link{QuantilePG-constructor}}] \item Documentation
  of a method dispaching to an object of a certain S4 class
  is named by the name of the method, followed by ``-'',
  followed by the name of the Class.  [cf.
  \code{\link{getValues-QuantilePG}}] }
}
\author{
Tobias Kley
}
\references{
Kley, T. (2014a). Quantile-Based Spectral Analysis of
Stationary Time-Series. unpublished Ph.D. Dissertation,
Ruhr University Bochum.

Kley, T. (2014b). An Object-oriented Framework for
Quantile-based Spectral Analysis and a Reference
Implementation in R: The quantspec Package. in preparation.

Dette, H., Hallin, M., Kley, T. & Volgushev, S. (2014+). Of
Copulas, Quantiles, Ranks and Spectra: an
\eqn{L_1}{L1}-approach to spectral analysis.
\emph{Bernoulli}, \bold{forthcoming}.

Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2014).
Quantile Spectral Processes: Asymptotic Analysis and
Inference. \url{http://arxiv.org/abs/1401.8104}.

Oracle (2014). Coding conventions of the Java Programming
Language.
\url{http://www.oracle.com/technetwork/java/codeconv-138413.html}.
}

