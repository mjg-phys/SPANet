# SPANet for 4Tops

Forked library with additional tools to create and use SPANet models to reconstruct four top quark (tttt) decays in the 2LSS, 3L and 4L channels with ATLAS data.

**Please first read [README2.md](README2.md) for instructions on how to install and use this library.** Here are additional instructions specific to use with ATLAS 4Tops data.

Created by Dean Reiter, Cornell University '25.


## Introduction

This library provides helper scripts needed to create and use `SPANet` models to reconstruct tttt decays in ATLAS data. With this library, tttt Monte Carlo data stored in the n-tuple format in `.root` files can be used to train `SPANEt` models. These models can then be evaluated on real tttt n-tuple data, and the models will attempt to reconstruct the parent top quark and W boson particles in the decays. Lastly, the kinematics of the parent particles, jet and lepton assignments, and reconstruction probablities can be stored as new variables in the n-tuple which was evaluated.

Separate `SPANet` models must be created and trained for each tttt decay channel (2LSS, 3L, 4L). As such, tttt data used for training and evaluation will be separated by channel. The tttt n-tuples used for training and evaluation do not need to be modified prior to using this library. However, after evaulating these `SPANet` models to reconstruct tttt decays, the n-tuples produced (with the reconstruction output variables added) *will* be separated by decay channel. If you wish to use this library to reconstruct all events in an n-tuple regardless of decay channel, you will need to create and train three `SPANet` models, evaluate them individually on the same n-tuple, and finally combine the three output n-tuples into one.

Data formatting scripts are provided for each decay channel in [data/multileptonic_tttt](data/multileptonic_tttt). These data formatting scripts take N number of n-tuples, combine them, select data from one decay channel, format the data, and create the `.h5` files needed for either training or evaluation. Event specification `.yaml` files are provided for each channel in [event_files](event_files). Training settings `.json` files are provided for each channel in [options_files/multileptonic_tttt](options_files/multileptonic_tttt), but these settings are not optimized. Output formatting scripts are provided for each channel in [predict](predict).

## Event Topology

Here is an example of the event topology for the 2LSS decay channel defined in [event_files/tttt_2lss_multi.yaml](event_files/tttt_2lss_multi.yaml): 

```yaml
INPUTS:
    SEQUENTIAL:
      Source:
        e: log_normalize
        pt: log_normalize
        eta: normalize
        phi: normalize
        btag: none
        mtag: none
        etag: none
        q: none 

EVENT:
  t1:
    - q1
    - q2
    - b
  t2:
    - q1
    - q2
    - b
  t3:
    - l
    - b
  t4:
    - l
    - b

PERMUTATIONS:
    EVENT:
      - [ t1, t2 ]
      - [ t3, t4 ]
    t1:
      - [ q1, q2 ]
    t2:
      - [ q1, q2 ]
```

The kinematic variables used as inputs are particle energy `e`, transverse momentum `pt`, angle `phi`, pseudorapidity `eta`, b-tagging variable at an 80 percent working point `btag`, muon or electron tags `mtag` and `etag`, and electric charge `q`.

The event particles we define are four top quarks `t1` through `t4`. The decay products of the hadronic top quarks (`t1` and `t2` in this example) are two jets `q1`, `q2`, and a b-jet `b`. The decay products of the leptonic quarks (`t3` and `t4` in this example) are a lepton `l` and a b-jet `b`. Note that neutrinos are not defined as decay products in the event topology since we do not have neutrino kinematics, and the model does not attempt to reconstruct them. 

The symmetries defined are between the hadronic top quarks, the leptonic top quarks, and the two hadronic top quark decay products `q1` and `q2`.

## Creating Training Data

## Configuring Training Options

## Formatting Evaluation Data

## Evaluation

## Output Formatting






