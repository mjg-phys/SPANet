import h5py
import sys
import numpy as np


# Set up h5 file strutcture 

og = h5py.File(str(sys.argv[1]), 'r')
hf = h5py.File(str(sys.argv[3]), 'w')

max_events = int(sys.argv[2])

inputs = hf.create_group('INPUTS')
targets = hf.create_group('TARGETS')
# classes = hf.create_group('CLASSIFICATIONS')

t1 = targets.create_group('t1')
t2 = targets.create_group('t2')


mask = np.array(og.get('INPUTS/Source/MASK'))
mask = mask[:max_events]

eta = np.array(og.get('INPUTS/Source/eta'))
eta = eta[:max_events]

btag = np.array(og.get('INPUTS/Source/btag'))
btag = btag[:max_events]

phi = np.array(og.get('INPUTS/Source/phi'))
phi = phi[:max_events]

pt = np.array(og.get('INPUTS/Source/pt'))
pt = pt[:max_events]

mass = np.array(og.get('INPUTS/Source/mass'))
mass = mass[:max_events]

source = inputs.create_group('Source')
mask_data = source.create_dataset('MASK', data=mask)
eta_data = source.create_dataset('eta', data=eta)
btag_data = source.create_dataset('btag', data=btag)
phi_data = source.create_dataset('phi', data=phi)
pt_data = source.create_dataset('pt', data=pt)
mass_data = source.create_dataset('mass', data=mass)


b1 = np.array(og.get('TARGETS/t1/b'))
b1 = b1[:max_events]

q11 = np.array(og.get('TARGETS/t1/q1'))
q11 = q11[:max_events]

q12 = np.array(og.get('TARGETS/t1/q2'))
q12 = q12[:max_events]

b2 = np.array(og.get('TARGETS/t2/b'))
b2 = b2[:max_events]

q21 = np.array(og.get('TARGETS/t2/q1'))
q21 = q21[:max_events]

q22 = np.array(og.get('TARGETS/t2/q2'))
q22 = q22[:max_events]

b1_data = t1.create_dataset('b', data=b1)
q11_data = t1.create_dataset('q1', data=q11)
q12_data = t1.create_dataset('q2', data=q12)
b2_data = t2.create_dataset('b', data=b2)
q21_data = t2.create_dataset('q1', data=q21)
q22_data = t2.create_dataset('q2', data=q22)

print("done copying")

hf.close()
og.close()

#user.dreiter.38180875._000001.output.root



#





