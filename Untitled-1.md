---
title: Lab Two -- Basic Tasks in R # Title of the document

# Bibliography Details
bibliography: Lab2References.bib  # Name of bib file in the same folder 
bibliographystyle: apa                   # Bibliography style

execute:
  warning: false
  message: false

# Below will do html by default. To render the entire project
# Use the following code:
# quarto::quarto_render("Lab2WriteUp.qmd")

format:
  html: 
    number-sections: TRUE
    code-fold: TRUE
    code-line-numbers: TRUE
    cap-location: top 
  pdf:
    documentclass: article
    number-sections: TRUE
    reference-section-title: "References" 
    code-line-numbers: TRUE
    cap-location: top 
    include-in-header:
      - text: |
          \usepackage{amsmath}
          \usepackage[dvipsnames]{xcolor}
          \usepackage[shortlabels]{enumitem}
          \usepackage{parskip}
          \usepackage{hyperref}
          \usepackage{geometry}
          \usepackage{float}
    geometry:
      - top=0.75in
      - left=0.75in
      - bottom=0.75in
      - right=0.75in
format-links: [pdf]
---
```{r}
#| echo: false

default_chunk_hook  <- knitr::knit_hooks$get("chunk")

latex_font_size <- c("Huge", "huge", "LARGE", "Large", 
                     "large", "normalsize", "small", 
                     "footnotesize", "scriptsize", "tiny")

knitr::knit_hooks$set(chunk = function(x, options) {
  x <- default_chunk_hook(x, options)
  if(options$size %in% latex_font_size) {
    paste0("\n \\", options$size, "\n\n", 
      x, 
      "\n\n \\normalsize"
    )
  } else {
    x
  }
})

```

```{r}
#| label: load-packages
#| include: false
# Load Libraries here
```

* Use `R` [@R] to complete the tasks below. Make sure to start your solutions in on a new line that starts with "**Solution**:".
* Make sure to use the Quarto Cheatsheet. This will make completing and writing up the lab *much* easier.

# Libraries and Packages

## Part a
Install the `esquisse` package for `R`. Report the code below, ensuring to set the `eval` option to `false`, using "\#| eval: false", so that we don't install the package every time you compile your Quarto document.

**Solution**

```{r}
#| eval: false
install.packages ("esquisse")
```

## Part b
Load the `esquisse` package for `R`. Report the code below. It does not matter whether you set "\#| eval: false" because loading the package does not add a lot of time. Note you will not evaluate any code that uses the `esquisse` package, so there is no need to load it as you compile your Quarto document.

**Solution**
```{r}
library ("esquisse")
```

## Part c
How can you ask for help about the  `esquisse` package for `R`? Report the code below, ensuring to set the `eval` option to `false`, using "\#| eval: false", so that the help document isn't loaded every time you compile your Quarto document.

**Solution**
```{r}
#| eval: false
help ("esquisse")
```

## Part d
Is a demo or vignette available for the  `esquisse` package for `R`? Report the code below, ensuring to set the `eval` option to `false`, using "\#| eval: false", so that the you don't get errors when you compile your Quarto document.

**Solution**
```{r}
#| eval: false
demo ("esquisse")
```
```{r}
#| eval: false
vignette ("esquisse")
```

## Part e
Add the BibTeX citation for `esquisse` package for `R` to your `.bib` file and cite it in a sentence describing what the package does below. Note you can reference a citation named `esquisse` using "[@\hspace{-0em}esquisse]".
<!-- Note the \hspace{-0em} above is for printing 
the command instead of actually citing. -->

**Solution**

We use `esquisse` [@esquisse] to explore and visualize data interactively.

## Part f
Run `esquisser()` in your console -- not in a chunk of `R` code in your Quarto document. 

i. Select `palmerpenguins` from the "Select an environment in which to search:" dropdown. This should automatically select `penguins` in the "Select a data.frame:" dropdown. Click "Import Data".
ii. Drag `body_mass_g` to the X box and `species` to the fill box. 
iii. Describe what you see in words. What can you conclude about Adelie, Chinstrap, and Gentoo penguins based on the resulting plot? 

**Solution**

Adelie penguins seem to be largest in quantity. Adelie and Chinstrap penguins, on average, have lower body mass than Gentoo penguines. 

# Objects and Vectors
Create the following vectors in `R`. In some cases, you may be able to use `seq()` or `rep()` while in others you cannot. Use these functions when possible, otherwise manually create the vector and explain why that was necessary.

**Some Snowday Notes**

There are three ways we will create a vector. Most simply, we can create one from scratch. For example, I create a vector of odds, and evens less than 10 below.

```{r}
#| size: "scriptsize"
#| code-fold: FALSE
(odds <- c(1, 3, 5, 7, 9))
(evens <- c(2, 4, 6, 8))
```

There are also built-in functions like `seq(...)` for doing this:
```{r}
#| size: "scriptsize"
#| code-fold: FALSE
(odds <- seq(from=1, to=9, by=2))
(evens <- seq(from=2, to=8, by=2))
```
Either approach is easy enough with a small number of elements, but what if I wanted odds less than 100? 1000? 1 million? The `rep(...)` and `seq()` approaches would be far more efficient.

We can also create repeating sequences by hand
```{r}
#| size: "scriptsize"
#| code-fold: FALSE
(repeating.seq1 <- c(1, 2, 3, 1, 2, 3, 1, 2, 3))
(repeating.seq2 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3))
```
or using a built-in function `rep(...)`
```{r}
#| size: "scriptsize"
#| code-fold: FALSE
(repeating.seq1 <- rep(c(1,2,3), times=3))
(repeating.seq2 <- rep(c(1,2,3), each=3))
```

There are some vectors for which we can't use a `seq()` or `rep()`. For example, consider the prime numbers less than 10. The primes are not sequential (e.g., jump by a fixed amount), nor are they repeating.
```{r}
#| size: "scriptsize"
#| code-fold: FALSE
primes <- c(2, 3, 5, 7)
```
**Note**: There is a `primes` package for `R` [@primes] that contains a `generate_primes(min, max)` function for generating a vector of primes from `min` to `max`.


## Part a
The Fibonacci Sequence is a recursive formula:
$$F_n = F_{n-1} + F_{n-2}$$
where $F_0=0$ and $F_1=1$. 

Create a vector of the first 10 Fibonacci numbers.

**Solutions**

I can't use `seq(...)` because the distance between 2 observations in this vector is not fixed.

```{r}
#| size: "scriptsize"
#| code-fold: FALSE
(fibonacci <- c(0, 1, 1, 2, 3, 5, 8, 13, 21, 34))
```

## Part b
Triangular Numbers are the cumulative sums of natural numbers:
$$F_n = \frac{n(n+1)}{2}.$$

Create a vector of the first 10 Triangular Numbers. 

**Solutions**

```{r}
#| size: "scriptsize"
#| code-fold: FALSE
(triangular <- c(1, 3, 6, 10, 15, 21, 28, 36, 45, 55))
```

## Part c
Suppose I were designing a repeated meausres experiment with three treatment conditions. Each of $n=10$ participants (with `ID` 1 to 10) will receive *all* experimental conditions, call them "Control", "Treatment A", and "Treatment B". 

Consider setting up data entry for this experiment.

i. Create a vector containing each `ID` repeated three times, once for each treatment.
ii. Create a vector containing each `Treatment` repeated for each participant ID.

**Solutions**

```{r}
#| size: "scriptsize"
#| code-fold: FALSE
(idrepeated <- rep(1:10, each=3))
(txrepeated <- rep(c("Control", "Treatment A", "Treatment B"), each=10))
```

## Part d
Create a vector containing the `character` "MATH" and the `numeric` 240. What is the resulting class? Explain why in a sentence.

**Solutions**

The resulting class is character because that class is to the right of the numeric class.