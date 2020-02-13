# Missing Data Analysis

## 1. Intro
불완전 자료 분석에 대해서 공부한 repository입니다.
불완전 자료에 대한 분석으로 많은 사람들이 적절한 값을 채워 넣는 것을 생각합니다. 이는 어느정도 맞는 말입니다. 결측치를 적절한 값으로 채워 넣는 imputation도 불완전 자료 분석의 한 분야죠. 하지만 이것외에도 할 수 있는 것이 생각보다 많습니다.

우리는 눈에 보이는 것을 주로 결측치로 간주하는데, 관측되지는 않았지만 내재하는 무엇인가를 가정하고 이를 missing data로 간주하여 분석을 진행할 수 있습니다. 예를 들어 어떤 정치인의 정치 성향은 눈으로 보이지는 않지만 충분히 존재하는 것이고, 이를 missing data로 보고 분석하는 것입니다. 이와 같이 응용할 수 있는 점이 많이 있어서, 이를 정리해두려고 합니다.

## 2. Contents
순서는 아래와 같습니다.
- [x] EM 알고리즘
 - Gaussian Mixture Ver. 1
 - Gaussian Mixture Ver. 2
 - Ideal Point Estimation (정치인들의 정치 성향 분석)
- [ ] Imputation
 - MICE
