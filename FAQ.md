<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [How do I calculate p values from Z scores?](#how-do-i-calculate-p-values-from-z-scores)

<!-- markdown-toc end -->

# How do I calculate p values from Z scores?

In R, do
```
p = pnorm(-abs(z)) * 2
```
