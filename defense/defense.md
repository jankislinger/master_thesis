<style>
.reveal h1, .reveal h2, .reveal h3 {
  word-wrap: normal;
  -moz-hyphens: none;
}
.small-code pre code {
  font-size: 1em;
}
.midcenter {
    position: fixed;
    top: 50%;
    left: 50%;
}
.footer {
    color: black; background: #E8E8E8;
    position: fixed; top: 90%;
    text-align:center; width:100%;
}
.pinky .reveal .state-background {
  background: #FF69B4;
} 
.pinky .reveal h1,
.pinky .reveal h2,
.pinky .reveal p {
  color: black;
}
</style>


Dynamic Fare Model
========================================================
author: Jan Kislinger
date: 1. 2. 2017
autosize: true

Thesis Content
========================================================

Theoretical part:

- Inhomogenous Markov process
- Estimating intensity of MP
- Simulated optimization

Practical applications:

-

Slide With Code
========================================================


```r
summary(cars)
```

```
     speed           dist       
 Min.   : 4.0   Min.   :  2.00  
 1st Qu.:12.0   1st Qu.: 26.00  
 Median :15.0   Median : 36.00  
 Mean   :15.4   Mean   : 42.98  
 3rd Qu.:19.0   3rd Qu.: 56.00  
 Max.   :25.0   Max.   :120.00  
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-2](defense-figure/unnamed-chunk-2-1.png)
