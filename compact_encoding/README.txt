- doc/ contains the description of the alt (new) encoding.

- psvn/ contains psvn code for generating the whole state space and features
  that make up the input for the encoder; only blocksworld and gripper has been 
  encoded in psvn.

- encoder/ contains the encoder/decoders. The former produces a sat theory
  while the latter reads a model (output of minisat) and decodes it.

