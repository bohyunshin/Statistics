## Longitudinal Data Anlaysis

### 1. Intro

경시적 자료 분석에 대해서 공부한 repository입니다. 우선 경시적 자료에 대해서 알아보겠습니다. 위키피디아에서는 경시적 자료를 아래와 같이 정의합니다.

*A longitudinal study is a rsearch design that involves **repeated observations** of the same variables over short or long periods of time*

여기서 핵심은 **repeated observations**이라는 것입니다. 경시적 자료는 한 subject에 대해서 여러번 관측한 자료를 말합니다. 예를 들어, 약의 효능을 알아보기 위해서 시험군과 대조군을 총 3일에 걸쳐 다섯 번 관찰하는 경우, 한 subject에 대해서 다섯 번의 관측치가 생기기 때문에 경시적 자료라고 볼 수 있습니다.

통계학에서 경시적 자료를 어떻게 바라보는지 알아보기 이전에, 가장 기본적인 자료에 대해서 선형 회귀를 하는 경우를 생각해보겠습니다. 선형 회귀는 통계학에서 자료를 분석할 때, 시작점이라고 해도 과언이 아니기 때문입니다. 선형 회귀에서 관심이 있는 종속변수 $y_i$ 는 서로 독립이라고 가정합니다. 이는 오차항에 대한 $iid$ 가정 때문입니다. 하지만 경시적 자료에서 더이상 $y_i$는 독립이 아닙니다. 왜냐하면 어떤 subject A의 반복 측정 자료는 서로 독립이 아니라 correlate 될 수밖에 없기 때문입니다. 따라서 경시적 자료에 대해서 선형 회귀를 수행한다면 가정 사항이 위배되므로 regression fit이 안 좋게 나온다든가, 변수가 유의하게 나오지 않는 등의 문제가 발생할 것입니다.

경시적 자료를 분석하는 통계적 방법으로 통계학에서 mixed model을 사용하며 단순하게 random effect을 포함하는 model, covariance pattern model과 random coefficient model이 있습니다.

* random effect을 포함하는 model
  * 단순하게 반복 측정된 subject 단위를 random effect로 두는 방법이다.
  * 이는 곧, 한 subject 내에서 측정된 데이터들 간에 동일한 correlation, 즉 compound symmetry을 가정하는 방법이다.
  * 아래의 covariance pattern model에서 pattern을 cs로 둔 것과 결과가 동일하다.
  * cs말고 다른 pattern을 가정하고 싶으면 covariance pattern model을 사용한다. 실제로 cs가 적절하지 않을 수 있기 때문이다.
  * 사용 R 패키지 (자료의 형태를 막론하고 모두 lme4 package에 구현되어 있다)
    * continuous: lme4 package
    * binary: lme4 package
    * count: lme4 package
  
* Covariance pattern model
  * random effect을 사용하지 않는다.
  * correlation 구조를 직접 지정하기 때문에 이해하기 쉽다.
  * 측정 시점이 모든 subject에게서 동일해야만 한다.
  * 사용 R 패키지
    * continous: R의 gls function
    * binary: R에서는 아직 못 찾음. SAS에서는 glimmix을 통해 가능
    * count: R에서는 아직 못 찾음. SAS에서는 glimmix을 통해 가능
    
* Random coefficient model
  * random effect을 사용한다.
  * 모형이 직관적으로 reasonable하고 측정 시점이 모든 subject에게서 동일할 필요가 없다.
  * correlation 구조가 복잡하며 이해하기 어렵다.
  * 사용 R 패키지
    * 찾는 중...

### 2. Contents

순서는 아래와 같습니다.

- [ ] Random effec Model
- [X] Covariance Pattern Model (2.27 완료)
- [ ] Random Coefficient Model
