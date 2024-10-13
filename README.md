# SPANet for 4Tops

Forked library with additional tools to create and use SPANet models to reconstruct four top quark (tttt) decays in the 2LSS, 3L and 4L channels with ATLAS data.

**Read [README2.md](README2.md) for instructions and tutorials on how to install and use the SPANet library before reading further.** Here are additional instructions specific to usung SPANet with ATLAS 4Tops data.

Created by Dean Reiter, Cornell University '25. 
Contact: dar333@cornell.edu


## Introduction

This library provides helper scripts needed to create and use `SPANet` models to reconstruct tttt decays in ATLAS data. With this library, tttt Monte Carlo data stored in the n-tuple format in `.root` files can be used to train `SPANEt` models. These models can then be evaluated on real tttt n-tuple data, and the models will attempt to reconstruct the parent top quark and W boson particles in the decays. Lastly, the kinematics of the parent particles, jet and lepton assignments, and reconstruction probablities can be stored as new variables back in the n-tuple which was evaluated.

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

To create a tttt training dataset, you must run one of the scripts `format_**_train.py` in the folder [data/multileptonic_tttt](data/multileptonic_tttt). The script takes as input as many n-tuples as you would like to chain together. The n-tuples must contain the reco-level truth matching variables from TopCPToolkit because the SPANet models need as much matching information as possible for training. The script outputs a `.h5` file with a chosen file name.

```
# 2LSS Channel
python data/multileptonic_tttt/format_2lss_train.py PATH/OUTPUT.h5 NTUPLE_1.root NTUPLE_2.root ...

# 3L Channel
python data/multileptonic_tttt/format_3l_train.py PATH/OUTPUT.h5 NTUPLE_1.root NTUPLE_2.root ...

# 4L Channel
python data/multileptonic_tttt/format_4l_train.py PATH/OUTPUT.h5 NTUPLE_1.root NTUPLE_2.root ...
```

## Configuring Training Options

Before training, you should modify the following files, depending on which decay channel you are creating a model to learn:
- 2LSS: [options_files/2lss.json](options_files/2lss.json)
- 3L: [options_files/3l.json](options_files/3l.json)
- 4L: [options_files/4l.json](options_files/4l.json)

## Formatting Evaluation Data

To format an n-tuple for evaluation with a SPANet model which you have trained, you must use one of the scripts `format_**_test.py` in the folder [data/multileptonic_tttt](data/multileptonic_tttt). The script takes as input one n-tuple and outputs a `.h5` file containing the kinematic inputs needed for `SPANet`. 

Note that the script only saves events with the correct number of leptons. For example, `format_2lss_test.py` selects events with two same-sign leptons and only saves their information in the `.h5` file.

```
# 2LSS Channel
python data/multileptonic_tttt/format_2lss_test.py OUTPUT.h5 NTUPLE.root

# 3L Channel
python data/multileptonic_tttt/format_3l_test.py OUTPUT.h5 NTUPLE.root 

# 4L Channel
python data/multileptonic_tttt/format_4l_test.py OUTPUT.h5 NTUPLE.root 
```

## Output Formatting

After creating a `SPANet` model, training it, and evaluating it on a tttt dataset, you will have a `.h5` file which was output by `spanet.predict`. This file contains the reconstruction information generated by SPANet, as well as the same kinematic information which was input. 

To store this reconstruction information in new branches in the n-tuple you evaluated, you can use the scripts `format_**_output.py` in the folder [predict/](predict/). This script creates a copy of the original n-tuple which was evaluated, selecting only 2LSS, 3L or 4L events. The script then reads the `.h5` file output by the SPANet evalation and stores the reconstruction information in this copied n-tuple.

```
# 2LSS Channel
python predict/format_2lss_output.py ORIGINAL_NTUPLE.root SPANET_OUTPUT.h5 NEW_NTUPLE.root

# 3L Channel
python predict/format_3l_output.py ORIGINAL_NTUPLE.root SPANET_OUTPUT.h5 NEW_NTUPLE.root

# 4L Channel
python predict/format_4l_output.py ORIGINAL_NTUPLE.root SPANET_OUTPUT.h5 NEW_NTUPLE.root
```

The script creates the following new variables and stores them in the final n-tuple:
- `b_index_SPANET`
    - Jet indices* of b-jets assigned to targets
    - Type: std::vector[int], Size: 4
    - Example: (1, 8, 3, 7)
- `q1_index_SPANET`
    - Jet indices of hadronic decay jet 1 assigned to targets
    - Type: std::vector[int], Size: 4
    - Example: (4, 6, -1, -1)**
- `q2_index_SPANET`
    - Jet indices of hadronic decay jet 2 assigned to targets
    - Type: std::vector[int], Size: 4
    - Example: (5, 2, -1, -1)
- `el_index_SPANET`
    - Electron indices* of leptonic decay electron assigned to targets
    - Type: std::vector[int], Size: 4
    - Example: (-1, -1, 0, -1)
- `mu_index_SPANET`
    - Muon indices* of leptonic decay muon assigned to targets
    - Type: std::vector[int], Size: 4
    - Example: (-1, -1, -1, 0)
- `top_isHadronic_SPANET`
    - 1 for hadronic target, 0 for leptonic target
    - Type: std::vector[int], Size: 4
    - Example: (1, 1, 0, 0)
- `top_m_SPANET`
    - Invariant mass of reconstructed hadronic top quarks, computed from kinematics of b, q1, q2 assigned to targets
    - Type: std::vector[float], Size: 4
    - Example: (142470, 184340, -1, -1)**
- `W_m_SPANET`
    - Invariant mass of reconstructed hadronic W bosons, computed from kinematics of q1, q2 assigned to targets
    - Type: std::vector[float], Size: 4
    - Example: (90450, 75610, -1, -1)
- `mbl_SPANET`
    - Invariant mass of b-jet and lepton system from leptonic decays, computed from kinematics of b, l assigned to targets
    - Type: std::vector[float], Size: 4
    - Example: (-1, -1, 93670, 114580)
- `top_assign_prob_SPANET`
    - Assignment probability for each target (roughly, probability that all assignments per target are correct - see SPANet documentation)
    - Type: std::vector[float], Size: 4
    - Example: (0.8423, 0.5682, 0.9534, 0.7588)
- `top_detect_prob_SPANET`
    - Detection probability for each target (see SPANet documentation)
    - Type: std::vector[float], Size: 4
    - Example: (0.4753, 0.7342, 0.6124, 0.8352)
- `top_margin_prob_SPANET`
    - Marginal probability for each target (see SPANet documentation)
    - Type: std::vector[float], Size: 4
    - Example: (0.1245, 0.3681, 0.0245, 0.2985)

*Jet index refers to index in TopCPToolkit variables `jet_pt`, `jet_eta`, etc. Electron index refers to index in TopCPToolkit variables `el_pt`, etc. Muon index refers to index in TopCPToolkit variables `mu_pt`, etc.

**Leptonic targets are given a `q1_index`, `q2_index`, `top_m`, and `W_m` of -1. Likewise, hadronic targets are given a `mu_index`, `el_index`, and `mbl` of -1. 


