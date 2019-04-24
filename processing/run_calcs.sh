pw.x < relax.in > relax.out
python3 ScfCreator.py relax.in relax.out
pw.x < scf.in > scf.out
pw.x < nscf.in > nscf.out
dos.x < dos.in > dos.out
projwfc.x < pdos.in > pdos.out
